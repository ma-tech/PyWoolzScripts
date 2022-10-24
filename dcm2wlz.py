#!/usr/bin/python3
##
# \file         dcm2wlz.py
# \author       Bill Hill
# \date         December 2019
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               MRC Institute of Genetics and Molecular Medicine,
#               University of Edinburgh,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C), [2019],
# The University Court of the University of Edinburgh,
# Old College, Edinburgh, UK.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be
# useful but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
# \brief        Creates Woolz spatial domain objects with values from
#               a DICOM image directory.
##

import os
import sys
import glob
import argparse
import math as m
import numpy as np
import pydicom as dcm
import ctypes as c
import Wlz as w

libc = c.CDLL('libc.so.6')

fopen = libc.fopen
fopen.argtypes = c.c_char_p, c.c_char_p,
fopen.restype = c.c_void_p

fwrite = libc.fwrite
fwrite.argtypes = c.c_void_p, c.c_size_t, c.c_size_t, c.c_void_p
fwrite.restype = c.c_size_t

fclose = libc.fclose
fclose.argtypes = c.c_void_p,
fclose.restype = c.c_int


verbose = False

class WlzError(Exception):
  pass

def vrbMsg(msg):
  global verbose
  if verbose:
    prog = sys.argv[0]
    print(prog + ': ' + msg, file=sys.stderr)

def wrnMsg(msg):
  prog = sys.argv[0]
  print(prog + ': Warning - ' + msg, file=sys.stderr)

# Create single Woolz image from DICOM slices
def makeWlzImageObj(slices, rescale): 
  vrbMsg('creating Woolz object')
  obj = None
  gvw = None
  errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
  s0 = slices[0]
  nx = s0.Columns
  ny = s0.Rows
  nz = len(slices)
  sx = float(s0.PixelSpacing[0])
  sy = float(s0.PixelSpacing[1])
  r_intercept = 0.0
  r_slope = 1.0
  bgd_v = 0
  if rescale:
    if ('RescaleIntercept' in s0) and ('RescaleSlope' in s0):
      r_intercept = float(s0.RescaleIntercept)
      r_slope = float(s0.RescaleSlope)
      bgd_v = int(r_intercept)
    else:
      wrnMsg('Unable to rescale image grey values to Hounsfield units as ' +
             'rescale parameters not in DICOM metadata.')
  sz = m.fabs(slices[1].ImagePositionPatient[2] - s0.ImagePositionPatient[2])
  if (sz < sys.float_info.epsilon) and ('SpacingBetweenSlices' in s0):
    sz = m.fabs(s0.SpacingBetweenSlices)
  if sz < sys.float_info.epsilon:
    if nz > 1:
      wrnMsg('Multiple slices at same position, slice thickness set to 1.0.')
    sz = 1.0
  x1 = int(np.floor(s0.ImagePositionPatient[0] / sx))
  y1 = int(np.floor(s0.ImagePositionPatient[1] / sy))
  z1 = int(np.floor(s0.ImagePositionPatient[2] / sz))
  g_type = w.WlzGreyType(w.WLZ_GREY_ERROR)
  if s0.BitsAllocated == 8:
    g_type = int(w.WLZ_GREY_UBYTE)
  elif s0.BitsAllocated == 16:
    g_type = int(w.WLZ_GREY_SHORT)
  else:
    raise Exception('Unsupported voxel grey type.')
  err_num = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
  vrbMsg('Creating cuboid object (' + 
          str(z1) + ', ' + str(z1 + nz - 1) + ', ' +
          str(y1) + ', ' + str(y1 + ny - 1) + ', ' +
          str(x1) + ', ' + str(x1 + nx - 1) +')')
  obj = w.WlzMakeCuboidI(z1, z1 + nz - 1, y1, y1 + ny - 1, x1, x1 + nx - 1,
             g_type, bgd_v, None, None, c.byref(err_num))
  if not bool(err_num):
    err_num = w.enum__WlzErrorNum(
              w.WlzSetVoxelSize(obj, sx, sy, sz))
  if not bool(err_num):
    vvp = obj.contents.values.vox.contents.values
    for iz in range(0, nz):
      vrbMsg('Setting values for slice ' + str(iz)+ ' of ' + str(nz))
      si = slices[iz]
      vp = vvp[iz].r.contents.values
      for iy in range(0, ny):
        offset = nx * iy
        for ix in range(0, nx):
          if g_type == int(w.WLZ_GREY_UBYTE):
            ubp = vp.ubp[offset + ix]
            ubp.contents = c.c_char(si.pixel_array[iy,ix])
          elif g_type == int(w.WLZ_GREY_SHORT):
            vp.shp[offset + ix] = c.c_short(
                int((si.pixel_array[iy,ix] * r_slope) + r_intercept))
          else:
            raise Exception('Unsupported voxel grey type.')
    vrbMsg('Object complete.')
  return obj
    
# If the given path is a regular file add it otherwise if it's a directory
# find DICOM image files in the directory
def getFiles(img_path):
  vrbMsg('getting image files from path ' + img_path)
  img_files = []
  if os.path.isdir(img_path):
    for f in glob.glob(img_path + '/*', recursive=False):
      vrbMsg('Reading DICOM file ' + f)
      img_files.append(dcm.read_file(f))
  else:
    img_files.append(dcm.read_file(img_path))
  return img_files

# Collect the files for each DICOM image series
def collectSeries(img_files):
  vrbMsg('collecting image slices from files ' + str(img_files))
  img_series = {}
  for f in img_files:
    if hasattr(f, 'SliceLocation') and hasattr(f, 'SeriesNumber'):
      n = f.SeriesNumber
      if not (n in img_series):
        img_series[n] = []
      img_series[n].append(f)
  return img_series

# Ensure that the slices of each series are in the correct order
def sortSlices(img_series):
  vrbMsg('sorting slices')
  for i in img_series:
    img_series[i] = sorted(img_series[i], key=lambda s: s.SliceLocation)

# Write Woolz object and text files
def outputFiles(out_dir, img_series, rescale, nowrite):
  for i in img_series:
    vrbMsg('Encoding ' + img_series[i][0].ProtocolName + '_' + str(i))
    obj = makeWlzImageObj(img_series[i], rescale)
    file_base = out_dir + '/' + img_series[i][0].ProtocolName + '_' + str(i)
    wlz_file_nm = file_base + '.wlz'
    txt_file_nm = file_base + '.txt'
    if not nowrite:
      fp = c.cast(fopen(wlz_file_nm.encode('utf-8'), b'wb'), c.POINTER(w.FILE))
      err_num = w.WlzWriteObj(fp, obj)
      if bool(err_num):
        raise WlzError()
      fclose(fp)
      vrbMsg('Woolz object written to ' + wlz_file_nm)
    w.WlzFreeObj(obj)
    obj = None
    if not nowrite:
      f = open(txt_file_nm, 'wt')
      f.write(str(img_series[i][0]))
      f.close()
      vrbMsg('Text written to ' + txt_file_nm)

# Parse the command line arguments
def parseArgs():
  parser = argparse.ArgumentParser(description = 
  'Creates Woolz spatial domain objects with values from a DICOM image or ' +
  'DICOM image directory. The file names are set from within the DICOM file' +
  'using the protocol name and image series index (<protocol>_<index>.wlz).' +
  'The DICOM metadata is also saved to a text file (<protocol>_<index>.txt).')
  parser.add_argument('-o', '--outdir',
      type=str, required=True,
      help='Output directory for Woolz object files.')
  parser.add_argument('-n', '--nowrite',
      action='store_true', default=False,
      help='Don\'t write any files (mainly useful for debugging).')
  parser.add_argument('-r', '--rescale',
      action='store_true', default=False,
      help='Rescale the image grey values using its grey value rescale ' +
      'parameters (use to get values as Hounsfield units where appropriate.')
  parser.add_argument('-v', '--verbose',
      action='store_true', default=False,
      help='Verbose output (mainly useful for debugging).')
  parser.add_argument('inpath',
      help='Input DICOM directory or image.')
  args = parser.parse_args()
  return(args)

def main():
  global verbose
  args = parseArgs()
  verbose = args.verbose
  img_files = getFiles(args.inpath)
  img_series = collectSeries(img_files)
  sortSlices(img_series)
  outputFiles(args.outdir, img_series, args.rescale, args.nowrite)

if __name__ == '__main__':
  main()

