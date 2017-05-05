#!/usr/bin/python
##
# \file         WlzITKAffineTransform.py
# \author       Bill Hill
# \date         April 2017
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
# \brief        Converts affine transform representations between ANTs/ITK and
#               Woolz.
##

from __future__ import print_function
import os
import re
import sys
import ctypes as c
import numpy as np
import argparse
from xml.etree import ElementTree as et
import Wlz as w

libc = c.CDLL('libc.so.6')

class WlzError(Exception):
  pass


def ErrorMsg(msg): #{
  print(prog + ': ' + msg, file=sys.stderr)
  exit(1)
#}

def VerbMsg(msg): #{
  if(args.verbose): #{
    print(prog + ': ' + msg)
  #}
#}

def ParseArgs(): #{
  parser = argparse.ArgumentParser(description = \
      'Converts affine transform representations between ANTs/ITK and ' +
      'Woolz. ' +
      'Reading / writting ITK transforms does not yet work. ' +
      'Note ANTs expects an inverse transform for images and a ' +
      'forward transform for coordinates.')
  parser.add_argument('-i', '--invert', \
      action='store_true', default=False, \
      help='invert the affine transform before writting it out.')
  parser.add_argument('-I', '--intype', \
      type=str, default='w', \
      help='Valid input types: ' +
           '  l - Landmarks to an ITK affine transform ' +
           '  k - ITK affine transform to a ascii homogeneous matrix ' +
           '  m - 4x4 ascii homogeneous matrix ' +
           '  w - 3D oolz affine transform to an ITK affine transform ')
  parser.add_argument('-O', '--outtype', \
      type=str, default='w', \
      help='Valid input types: ' +
           '  k - ITK affine transform to a ascii homogeneous matrix ' +
           '  m - 4x4 ascii homogeneous matrix ')
  parser.add_argument('-c', '--cor', \
      type=str, default='0,0,0', \
      help='Centre of rotation to be used when outputting an ITK affine ' +
           'transform with the form x,y,z. ')
  parser.add_argument('-o', '--output', \
      type=str, default='-', \
      help='output transform file.')
  parser.add_argument('-v', '--verbose', \
      action='store_true', default=False, \
      help='verbose output (mainly useful for debugging).')
  parser.add_argument('input', help = \
      'input transform file.')
  args = parser.parse_args()
  return(args,parser)
#}

def Invert(mat): #{
  VerbMsg('Inverting matrix ' + str(mat))
  mat = np.linalg.inv(mat)
  VerbMsg('Inverted matrix ' + str(mat))
  return(mat)
#}

def ReadITKTr(): #{
  try: #{
    if args.input == '-': #{
      f = sys.stdin
    else: #}{
      f = open(args.input, 'r')
    #}
    for rec in f: #{
      if rec.startswith('Parameters:'): #{
        prm = rec.split(':')[1].lstrip().rstrip().split(' ')
      #}
      if rec.startswith('FixedParameters:'): #{
        fpm = rec.split(':')[1].lstrip().rstrip().split(' ')
      #}
    #}
    if not args.input == '-': #{
      f.close()
    #}
    VerbMsg('Parameters ' + str(prm))
    VerbMsg('FixedParameters ' + str(fpm))
    if ((not len(prm) == 12) or (not len(fpm) == 3)): #{
      raise IOError()
    #}
    mat = np.array( \
          [[float(prm[0]), float(prm[1]), float(prm[2]), float(prm[9])], \
           [float(prm[3]), float(prm[4]), float(prm[5]), float(prm[10])], \
           [float(prm[6]), float(prm[7]), float(prm[8]), float(prm[11])], \
           [0.0,           0.0,           0.0,           1.0]])
    cen = np.array([float(fpm[0]), float(fpm[1]), float(fpm[2])])
    #Assume that the ITK transform uses 
    for i in range(0, 3): #{
      mat[i][3] = cen[i] + mat[i][3] - \
                  (mat[i][0] * cen[0] + \
                   mat[i][1] * cen[1] + \
                   mat[i][2] * cen[2])
    #}
    VerbMsg('mat ' + str(mat))
  except IOError: #}{
    ErrorMsg('Failed to read ITK affine transform from file ' + args.input)
  #}
  return(mat)
#}

def ReadMatTr(): #{
  try: #{
    if args.input == '-': #{
      f = sys.stdin
    else: #}{
      f = open(args.input, 'r')
    #}
    ln = 0
    mat = np.zeros((4,4))
    for rec in f: #{
      prm = re.sub(' +',' ', rec).lstrip().rstrip().split(' ')
      VerbMsg(str(prm))
      if not len(prm) == 4: #{
        raise IOError()
      #}
      for i in range(0, 4): #{
        mat[ln][i] = float(prm[i])
      #}
      ln = ln + 1
    #}
    VerbMsg(str(ln))
    if not args.input == '-': #{
      f.close()
    #}
    if not ln == 4: #{
      raise IOError()
    #}
  except IOError: #}{
    ErrorMsg('Failed to read affine transform homogeneous matrix from file ' \
             + args.input)
  #}
  return(mat)
#}

def GetLmkVertex(g): #{
  vtx = w.WlzDVertex3()
  vtx.vtX = float(g.find('X').text)
  vtx.vtY = float(g.find('Y').text)
  vtx.vtZ = float(g.find('Z').text)
  return(vtx)
#}

def WlzTrToMat(tr): #{
  m = tr.contents.mat
  mat = np.zeros((4,4))
  for i in range(0, 4): #{
    for j in range(0, 4): #{
      mat[i][j] = m[i][j]
    #}
  #}
  return(mat)
#}

def ReadLmkTr(): #{
  sv = []
  tv = []
  try: #{
    if args.input == '-': #{
      f = sys.stdin
    else: #}{
      f = open(args.input, 'r')
    #}
    doc = et.parse(f)
    idx = 0
    for lm in doc.findall('Landmark'): #{
      s = lm.find('Source')
      t = lm.find('Target')
      sv.append(GetLmkVertex(s))
      tv.append(GetLmkVertex(t))
      VerbMsg('{:d} {:g} {:g} {:g} {:g} {:g} {:g}'. \
              format( \
                idx, \
                sv[idx].vtX, sv[idx].vtY, sv[idx].vtZ, \
                tv[idx].vtX, tv[idx].vtY, tv[idx].vtZ))
      idx = idx + 1
    #}
    if not args.input == '-': #{
      f.close()
    #}
    if ((not len(sv) == len(tv)) or (len(sv) < 3)): #{
      raise IOError()
    #}
  except IOError: #}{
    ErrorMsg('Failed to read ITK affine transform from file ' + args.input)
  #}
  try: #{
    n = len(sv)
    svtx = (w.WlzDVertex3 * n)()
    tvtx = (w.WlzDVertex3 * n)()
    errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
    for i in range(0, n): #{
      svtx[i] = sv[i]
      tvtx[i] = tv[i]
    #}
    tr = w.WlzAffineTransformLSq3D(n, tvtx, n, svtx, 0, None, \
                         w.enum__WlzTransformType(w.WLZ_TRANSFORM_3D_AFFINE), \
                         c.byref(errNum))
    if bool(errNum) or \
       (not int(tr.contents.type) == int(w.WLZ_TRANSFORM_3D_AFFINE)): #{
      raise WlzError()
    #}
  except WlzError: #}{
    ErrorMsg('Failed to compute LSq affine transform from landmarks (' +
             w.WlzStringFromErrorNum(errNum, None) + ')')
  #}
  mat = WlzTrToMat(tr)
  return(mat)
#}

def ReadWlzTr(): #{
  try: #{
    errNum = w.enum__WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
    if args.input== '-': #{
      f = libc.stdin
    else: #}{
      f = libc.fopen(args.input, 'rb')
    #}
    if not bool(f): #{
      raise IOError()
    #}
    obj = w.WlzAssignObject( \
          w.WlzReadObj(f, c.byref(errNum)), None)
    if not args.input == '-': #{
      libc.fclose(f)
    #}
    if not bool(obj): #{
      errNum = w.enum__WlzErrorNum(w.WLZ_ERR_OBJECT_NULL)
      raise WlzError()
    #}
    if ((not int(obj.contents.type) == int(w.WLZ_AFFINE_TRANS)) or \
        (not int(obj.contents.domain.core.contents.type) ==  \
             int(w.WLZ_TRANSFORM_3D_AFFINE))): #{
      errNum = w.enum__WlzErrorNum(w.WLZ_ERR_OBJECT_TYPE)
      raise WlzError()
    #}
  except IOError: #}{
    ErrorMsg('Failed to open ' + args.input)
  except WlzError: #}{
    ErrorMsg('Failed to read affine transform from ' + args.input +
             '(' + w.WlzStringFromErrorNum(errNum, None) + ')')
  #}
  mat = WlzTrToMat(obj.contents.domain.t)
  return(mat)
#}

def WriteITKTr(mat): #{
  try: #{
    if args.output == '-': #{
      f = sys.stdout
    else: #}{
      f = open(args.output, 'w')
    #}
    VerbMsg(prog + ': ' + str(args.cor))
    print('#Insight Transform File V1.0', file=f)
    print('#Transform 0', file=f)
    print('Transform: MatrixOffsetTransformBase_double_3_3', file=f)
    print(('Parameters: ' +
           '{:g} {:g} {:g} ' +
           '{:g} {:g} {:g} ' +
           '{:g} {:g} {:g} ' +
           '{:g} {:g} {:g} '). \
          format(mat[0][0], mat[0][1], mat[0][2], \
                 mat[1][0], mat[1][1], mat[1][2], \
                 mat[2][0], mat[2][1], mat[2][2], \
                 args.cor[0] + mat[0][3], \
                 args.cor[1] + mat[1][3], \
                 args.cor[2] + mat[2][3]), \
                 file=f)
    #}
    print(('FixedParameters: {:g} {:g} {:g}'). \
          format(args.cor[0], args.cor[1], args.cor[2]),
          file=f)
    if not args.output == '-': #{
      f.close()
    #}
  except IOError: #}{
    ErrorMsg('Failed to write transform matrix to file.' +  args.output)
#}

def WriteMatTr(mat): #{
  try: #{
    if args.output == '-': #{
      f = sys.stdout
    else: #}{
      f = open(args.output, 'w')
    #}
    for i in range(0, 4): #{
      print('{:g} {:g} {:g} {:g}'.\
            format(mat[i][0], mat[i][1], mat[i][2], mat[i][3]), file=f)
    #}
    if not args.output == '-': #{
      f.close()
    #}
  except IOError: #}{
    ErrorMsg('Failed to write transform matrix to file.' +  args.output)
#}

if __name__ == '__main__': #{
  # Process the command line
  usage = False
  prog = sys.argv[0];
  args,parser = ParseArgs()
  if args.outtype == 'k': #{
    try: #{
      cor = re.sub(' +','', args.cor).split(',')
      if not len(cor) == 3: #{
        raise IOError()
      else: #}{
        args.cor = [float(cor[0]), float(cor[1]), float(cor[2])]
      #}
    except: #}{
      usage = True
    #}
  #}
  VerbMsg(prog + ': args = ' + str(args))
  if not usage: #{
    errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
    if args.intype == 'l': #{
      mat = ReadLmkTr()
    #}
    elif args.intype == 'k': #{
      mat = ReadITKTr()
    #}
    elif args.intype == 'm': #{
      mat = ReadMatTr()
    #}
    elif args.intype == 'w': #{
      mat = ReadWlzTr()
    #}
    else: #{
      usage = True
    #}
  #}
  if (not usage) and (not args.invert): #{
    mat = Invert(mat)
  #}
  if not usage: #{
    if args.outtype == 'k': #{
      WriteITKTr(mat)
    #}
    elif args.outtype == 'm': #{
      WriteMatTr(mat)
    #}
    else: #{
      usage = True
  #}
  if usage: #{
    parser.print_help()
    exit(1)
  #}
#}
