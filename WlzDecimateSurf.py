#!/usr/bin/python
##
# \file         WlzDecimateSurf.py
# \author       Bill Hill
# \date         July 2015
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               MRC Institute of Genetics and Molecular Medicine,
#               University of Edinburgh,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C), [2015],
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
# \brief        Decimates, smooths and flips normals for surfaces.
#               This script makes use of Woolz binaries, PyWoolz
#               and MeshLab.
##

from __future__ import print_function
import os
import sys
import argparse
import subprocess
import tempfile
import ctypes
import glob
import Wlz

libc = ctypes.CDLL("libc.so.6")

libc.fopen.restype = ctypes.POINTER(Wlz.FILE)

ma_bin_dir        = '/opt/MouseAtlas/bin'
WlzExtFFConvert   = ma_bin_dir + '/WlzExtFFConvert'
MeshLabServer     = '/opt/vis/bin/meshlabserver'

def ParseArgs():
  parser = argparse.ArgumentParser(description= \
  'Decimates, smooths and flips normals for surfaces.')
  parser.add_argument('-f', '--flip', \
      action='store_true', default=False, \
      help='Flip face normals')
  parser.add_argument('-o', '--outfile', \
      type=str, required=True,\
      help='Output surface file.')
  parser.add_argument('-m', '--maxface', \
      type=int, default=20000, \
      help='Maximum number of faces in the surface triangulation.')
  parser.add_argument('-s', '--smooth', \
      action='store_true', default=False, \
      help='Use Taubin smoothing.')
  parser.add_argument('-t', '--tmpdir', \
      type=str, default='/tmp', \
      help='Temporary directory for working files.')
  parser.add_argument('-v', '--verbose', \
      action='store_true', default=False, \
      help='Verbose output (mainly useful for debugging).')
  parser.add_argument('infile',
      help='Input Woolz domain.')
  args = parser.parse_args()
  return(args)


def ReadWoolzFile(filename):
  obj = None
  errNum = Wlz.enum__WlzErrorNum(Wlz.WLZ_ERR_FILE_OPEN)
  fp = libc.fopen(filename, 'rb')
  if(bool(fp)):
    obj = Wlz.WlzReadObj(fp, ctypes.byref(errNum))
    libc.fclose(fp)
  return(obj)


def WriteWoolzObject(filename, obj):
  errNum = Wlz.enum__WlzErrorNum(Wlz.WLZ_ERR_FILE_OPEN)
  fp = libc.fopen(filename, 'wb')
  if(bool(fp)):
    errNum = Wlz.WlzWriteObj(fp, obj)
    libc.fclose(fp)
  return(errNum)


def StartMeshLabFilter(flt):
  flt.append('<!DOCTYPE FilterScript>')
  flt.append('<FilterScript>')
  flt.append('<filter name="Merge Close Vertices">')
  flt.append('<Param type="RichAbsPerc" value="0.1" min="0" ' + \
             'name="Threshold" max="16"/>')
  flt.append('</filter>')


def EndMeshLabFilter(flt):
  flt.append('</FilterScript>')


def AppendMeshLabFilter(flt, target_face_count, taubin_steps):
  if(target_face_count > 0):
    flt.append('<filter name="Quadric Edge Collapse Decimation">')
    flt.append('<Param type = "RichInt"   value = "' + \
               str(target_face_count) + '" ' + \
               'name = "TargetFaceNum"/>')
    flt.append('<Param type = "RichFloat" value = "0"     ' + \
               'name = "TargetPerc"/>')
    flt.append('<Param type = "RichFloat" value = "0.3"   ' + \
               'name = "QualityThr"/>')
    flt.append('<Param type = "RichBool"  value = "false" ' + \
               'name = "PreserveBoundary"/>')
    flt.append('<Param type = "RichFloat" value = "1"     ' + \
               'name = "BoundaryWeight"/>')
    flt.append('<Param type = "RichBool"  value = "false" ' + \
               'name = "PreserveNormal"/>')
    flt.append('<Param type = "RichBool"  value = "false" ' + \
               'name = "PreserveTopology"/>')
    flt.append('<Param type = "RichBool"  value = "true"  ' + \
               'name = "OptimalPlacement"/>')
    flt.append('<Param type = "RichBool"  value = "false" ' + \
               'name = "PlanarQuadric"/>')
    flt.append('<Param type = "RichBool"  value = "false" ' + \
               'name = "QualityWeight"/>')
    flt.append('<Param type = "RichBool"  value = "true"  ' + \
               'name = "AutoClean"/>')
    flt.append('<Param type = "RichBool"  value = "false" ' + \
               'name = "Selected"/>')
    flt.append('</filter>')
  if(taubin_steps > 0):
    flt.append('<filter name="Taubin Smooth">')
    flt.append('<Param type = "RichFloat" value = "0.5"   ' + \
               'name = "lambda"/>')
    flt.append('<Param type = "RichFloat" value = "-0.53" ' + \
               'name = "mu"/>')
    flt.append('<Param type = "RichInt"   value = "' + \
               str(taubin_steps) + '"    ' + \
               'name = "stepSmoothNum"/>')
    flt.append('<Param type = "RichBool"  value = "false" ' + \
               'name = "Selected"/>')
    flt.append('</filter>')


def WriteMeshLabFilter(filename, flt):
  with open(filename, "w") as ff:
    for fr in flt:
      print(fr, file=ff)
  ff.close()

def CleanExit(stat):
  for f in glob.glob(workfile + '[0-1].*'):
    os.remove(f)
  exit(stat)

if __name__ == '__main__':
  args = ParseArgs()
  prog = sys.argv[0];
  errNum = Wlz.enum__WlzErrorNum(Wlz.WLZ_ERR_NONE)

  if(args.verbose):
    print('Args = ' + str(args))

  workfile = tempfile.mktemp(dir=args.tmpdir, prefix='wd2vs')

  # Convert the input file to a Woolz contour object, read it in and then
  # find the number of faces in the model.
  cmdline = [WlzExtFFConvert, '-o' + workfile + '0.wlz', \
             args.infile]
  if(args.verbose):
    print(cmdline)
  rtn = subprocess.call(cmdline)
  if(bool(rtn)):
    print(prog + ': WlzExtFFConvert failed to convert input file to wlz.')
    CleanExit(1)
  obj = Wlz.WlzAssignObject(ReadWoolzFile(workfile + '0.wlz'), None)
  if(not bool(obj)):
    print(prog + ': Failed to read Woolz file.')
    CleanExit(1)
  model = obj.contents.domain.ctr.contents.model
  n_faces = model.contents.res.face.numElm
  if(args.verbose):
    print('Initial number of faces = ' + str(n_faces))
  if(n_faces < 4):
    print(prog + ': Invalid surface model (n_faces = ' + str(n_faces) + ')')
    CleanExit(1)

  # Flip face orientation if required.
  if(args.flip):
    if(args.verbose):
      print('Flipping face orientation.')
    errNum = Wlz.WlzGMFilterFlipOrient(model);
    if(bool(errNum)):
      print(prog + ': Failed to flip face orientation (' +
            Wlz.WlzStringFromErrorNum(errNum, None) + ').')
      CleanExit(1)

  # Write the Woolz object file using the stl file format.
  if(args.verbose):
    print('Writing working Woolz file.')
  if(bool(WriteWoolzObject(workfile + '0.wlz', obj))):
    print(prog + ': Failed to write working Woolz file.')
    CleanExit(1)
  Wlz.WlzFreeObj(obj)
  cmdline = [WlzExtFFConvert, '-o' + workfile + '0.stl', \
             workfile + '0.wlz']
  if(args.verbose):
    print(cmdline)
  rtn = subprocess.call(cmdline)
  if(bool(rtn)):
    print(prog + ': WlzExtFFConvert failed to convert working wlz file to stl.')
    CleanExit(1)

  # Create meshlab filter script to reduce the number of faces and smooth
  # the surface
  if(args.verbose):
    print('Creating MeshLab filter file.')
  flt = []
  if(args.smooth):
    taubin_steps = 10
  else:
    taubin_steps = 0
  StartMeshLabFilter(flt)
  if(n_faces > args.maxface):
    face_count = n_faces
    while(face_count > args.maxface):
      if(face_count > args.maxface * 4):
        face_count = face_count / 4
      else:
        face_count = args.maxface
      AppendMeshLabFilter(flt, face_count, taubin_steps)
      if(args.smooth):
        taubin_steps = 4
  else:
    AppendMeshLabFilter(flt, 0, taubin_steps)
  EndMeshLabFilter(flt)
  WriteMeshLabFilter(workfile + '0.mlx', flt)
  
  # Execute the MeshLab filter
  cmdline = [MeshLabServer,
             '-i ' + workfile + '0.stl', \
             '-o ' + workfile + '1.stl', \
             '-s ' + workfile + '0.mlx']
  if(args.verbose):
    print(cmdline)
  rtn = subprocess.call(cmdline)
  if(bool(rtn)):
    print(prog + ':meshlabserver failed to run filter script.')
    CleanExit(1)

  # Convert the STL file to required output format
  cmdline = [WlzExtFFConvert, '-o' + args.outfile, \
             workfile + '1.stl']
  if(args.verbose):
    print(cmdline)
  rtn = subprocess.call(cmdline)
  if(bool(rtn)):
    print(prog + ': WlzExtFFConvert failed to convert working stl file format.')
    CleanExit(1)

  CleanExit(0)
