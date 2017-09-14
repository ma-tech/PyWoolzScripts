#!/usr/bin/python
##
# \file         WlzNumpyArrayDemo.py
# \author       Bill Hill
# \date         September 2017
# \version      $Id$
# \par
# Address:
#               MRC Human Genetics Unit,
#               MRC Institute of Genetics and Molecular Medicine,
#               University of Edinburgh,
#               Western General Hospital,
#               Edinburgh, EH4 2XU, UK.
# \par
# Copyright (C), [2017],
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
# \brief        Demo of linking Woolz to Numpy through PyWoolz which
#               converts a 2 or 3D Woolz domain object either with or
#               without values to a Numpy array and then back to a
#               Woolz object again.
##
from __future__ import print_function
import os
import sys
import argparse as ap
import numpy as np
import ctypes as c
import Wlz as w
import skimage.morphology as skm
import skimage.filters as skf

libc = c.CDLL("libc.so.6")

class WlzError(Exception):
  pass

def ParseArgs(): #{
  parser = ap.ArgumentParser(description=
          'Demo of linking Woolz to Numpy through PyWoolz which ' +
          'converts a 2 or 3D Woolz domain object either with or ' +
          'without values to a Numpy array and then back to a ' +
          'Woolz object again. This is a very slow approach, better ' +
          'would be to write the Wlz-Numpy-Wlz conversions in c.')
  parser.add_argument('-o', '--outfile',
      type=str, required=True,
      help='Output object file.')
  parser.add_argument('-s ', '--skeletonize',
      action='store_true', default=False,
      help='Skeletonize the numpy array before converting back to Woolz ' +
           '(to demonstrate the use of skimage).')
  parser.add_argument('-v ', '--verbose',
      action='store_true', default=False,
      help='Verbose output.')
  parser.add_argument('infile',
      help='Input Woolz domain object.')
  args = parser.parse_args()
  return(args)
#}

def ErrorExit(msg): #{
  print(prog + ': ' + msg, file=sys.stderr)
  exit(1)
#}

def PrintWlzFacts(obj, f=sys.stderr): #{
  s = c.c_char_p(0)
  w.WlzObjectFacts(obj, None, c.byref(s), 0)
  print(str(s.value), file=f)
#}

def ReadWlzObj(filename): #{
  obj = None
  errNum = w.enum__WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
  try: #{
    fp = libc.fopen(filename, 'rb')
    if (bool(fp)): #{
      obj = w.WlzReadObj(fp, c.byref(errNum))
      libc.fclose(fp)
    #}
    if (bool(errNum)): #{
      raise WlzError()
    #}
  except: #}{
    ErrorExit('Failed to read Woolz object from file ' + filename + ' (' +
              w.WlzStringFromErrorNum(errNum, None) + ').')
  #}
  return(obj)
#}

def WriteWlzObj(filename, obj): #{
  errNum = w.enum__WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
  try: #{
    fp = libc.fopen(filename, 'wb')
    if (bool(fp)): #{
      errNum = w.WlzWriteObj(fp, obj)
      libc.fclose(fp)
    #}
    if (bool(errNum)): #{
      raise WlzError()
    #}
  except: #}{
    ErrorExit('Failed to write Woolz object to file ' + filename + ' (' +
              w.WlzStringFromErrorNum(errNum, None) + ').')
  #}
#}

def ArrayFromWlz(obj): #{
  ary = None
  org = [0]
  shape = (0)
  otype = None
  bbox = [0, 0, 0]
  errNum = w.WLZ_ERR_NONE
  if (not bool(obj)): #{
    errNum = w.WLZ_ERR_OBJECT_NULL
  elif ((obj.contents.type == w.WLZ_2D_DOMAINOBJ) or
        (obj.contents.type == w.WLZ_3D_DOMAINOBJ)): #}{
    otype = obj.contents.type
    if (not bool(obj.contents.domain.core)): #{
      errNum = w.WLZ_ERR_DOMAIN_NULL
    else: #}{
      if (bool(obj.contents.values.core)): #{
        errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
        gtype = w.WlzGreyTypeFromObj(obj, c.byref(errNum))
      else: #}{
        gtype = w.WLZ_GREY_BIT
      #}
    #}
  else: #}{
    errNum = w.WLZ_ERR_OBJECT_TYPE
  #}
  if (not bool(errNum)): #{
    errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
    bbox = w.WlzBoundingBox3I(obj, c.byref(errNum))
    if (not bool(errNum)): #{
      if (otype == w.WLZ_2D_DOMAINOBJ): #{
        org = [bbox.xMin, bbox.yMin]
        shape = (bbox.xMax - bbox.xMin + 1, bbox.yMax - bbox.yMin + 1)
      else: #}{
        org = [bbox.xMin, bbox.yMin, bbox.zMin]
        shape = (bbox.xMax - bbox.xMin + 1, bbox.yMax - bbox.yMin + 1, \
                 bbox.zMax - bbox.zMin + 1)
      #}
    #}
  #}
  if (not bool(errNum)): #{
    if (gtype == w.WLZ_GREY_BIT): #{
      atype = np.bool_
    elif (gtype == w.WLZ_GREY_UBYTE): #}{
      atype = np.uint8
    elif (gtype == w.WLZ_GREY_SHORT): #}{
      atype = np.int16
    elif (gtype == w.WLZ_GREY_INT): #}{
      atype = np.int32
    elif (gtype == w.WLZ_GREY_FLOAT): #}{
      atype = np.float32
    elif (gtype == w.WLZ_GREY_DOUBLE): #}{
      atype = np.float64
    else: #}{
      errNum = w.WLZ_ERR_GREY_TYPE
    #}
    if (not bool(errNum)): #{
      ary = np.zeros(shape, dtype=atype)
    #}
  #}
  if (not bool(errNum)): #{
    if (gtype == w.WLZ_GREY_BIT): #{
      if (otype == w.WLZ_2D_DOMAINOBJ): #{
        for ln in range(bbox.yMin, bbox.yMax + 1): #{
          y = ln - org[1]
          for kl in range(bbox.xMin, bbox.xMax + 1): #{
            x = kl - org[0]
            ary[x, y] = bool(w.WlzInsideDomain(obj, 0, ln, kl, None))
          #}
        #}
      else: #}{
        for pl in range(bbox.zMin, bbox.zMax + 1): #{
          z = pl - org[2]
          for ln in range(bbox.yMin, bbox.yMax + 1): #{
            y = ln - org[1]
            for kl in range(bbox.xMin, bbox.xMax + 1): #{
              x = kl - org[0]
              ary[x, y, z] = bool(w.WlzInsideDomain(obj, pl, ln, kl, None))
            #}
          #}
        #}
      #}
    else: #}{
      errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
      gvwsp = w.WlzGreyValueMakeWSp(obj, c.byref(errNum))
      if (not bool(errNum)): #{   
        if (otype == w.WLZ_2D_DOMAINOBJ): #{
          for ln in range(bbox.yMin, bbox.yMax + 1): #{
            y = ln - org[1]
            for kl in range(bbox.xMin, bbox.xMax + 1): #{
              x = kl - org[0]
              w.WlzGreyValueGet(gvwsp, 0, ln, kl)
              gv = gvwsp.contents.gVal[0]
              if (gtype == w.WLZ_GREY_UBYTE): #{
                ary[x, y] = gv.ubv
              elif (gtype == w.WLZ_GREY_SHORT): #}{
                ary[x, y] = gv.shv
              elif (gtype == w.WLZ_GREY_INT): #}{
                ary[x, y] = gv.inv
              elif (gtype == w.WLZ_GREY_FLOAT): #}{
                ary[x, y] = gv.flv
              elif (gtype == w.WLZ_GREY_DOUBLE): #}{
                ary[x, y] = gv.dbv
              #}
            #}
          #}
        else: #}{
          for pl in range(bbox.zMin, bbox.zMax + 1): #{
            z = pl - org[1]
            for ln in range(bbox.yMin, bbox.yMax + 1): #{
              y = ln - org[1]
              for kl in range(bbox.xMin, bbox.xMax + 1): #{
                x = kl - org[0]
                w.WlzGreyValueGet(gvwsp, pl, ln, kl)
                if (gtype == w.WLZ_GREY_UBYTE): #{
                  ary[x, y] = gvwsp.contents.gVal[0].ubv
                elif (gtype == w.WLZ_GREY_SHORT): #}{
                  ary[x, y] = gvwsp.contents.gVal[0].shv
                elif (gtype == w.WLZ_GREY_INT): #}{
                  ary[x, y] = gvwsp.contents.gVal[0].inv
                elif (gtype == w.WLZ_GREY_FLOAT): #}{
                  ary[x, y] = gvwsp.contents.gVal[0].flv
                elif (gtype == w.WLZ_GREY_DOUBLE): #}{
                  ary[x, y] = gvwsp.contents.gVal[0].dbv
                #}
              #}
            #}
          #}
      #}
      w.WlzGreyValueFreeWSp(gvwsp)
    #}
  #}
  return(ary, org, errNum)
#}

def ArrayToWlz(ary, org): #{
  obj = None
  otype = None
  shape = ary.shape
  dim = len(shape)
  atype = ary.dtype
  gsz = 0
  gtype = None
  val = w.WlzGreyP(0)
  errNum = w.WLZ_ERR_NONE
  if ((atype == np.bool) or (atype == np.bool_) or (atype == np.uint8)): #{
    gtype = w.WLZ_GREY_UBYTE
  elif ((atype == np.int8) or (atype == np.int16)): #}{
    gtype = w.WLZ_GREY_SHORT
  elif ((atype == np.uint16) or (atype == np.uint32) or
        (atype == np.uint64) or (atype == np.int32) or
        (atype == np.int64) or (atype == np.int_) or 
        (atype == np.intc) or (atype == np.intp)): #}{
    gtype = w.WLZ_GREY_INT
  elif ((atype == np.float16) or (atype == np.float32)): #}{
    gtype = w.WLZ_GREY_FLOAT
  elif (atype == np.float64) or (atype == np.float_): #}{
    gtype = w.WLZ_GREY_DOUBLE
  else: #}{
    errNum = w.WLZ_ERR_GREY_TYPE
  #}
  if not bool(errNum): #{
    gsz = w.WlzGreySize(gtype)
    errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
    if (dim == 2): #{
      val.v = w.AlcCalloc(shape[0] * shape[1], gsz)
      obj = w.WlzMakeRectI(org[1], org[1] + shape[1] - 1,
                           org[0], org[0] + shape[0] - 1,
                           gtype, val.inp, 0, None, None, c.byref(errNum))
    elif (dim == 3): #}{
      obj = w.WlzMakeCuboidI(org[2], org[2] + shape[2] - 1,
                             org[1], org[1] + shape[1] - 1,
                             org[0], org[0] + shape[0] - 1,
                             gtype, 0, None, None, c.byref(errNum))
    else: #}{
      errNum = w.WLZ_ERR_OBJECT_TYPE
    #}
  #}
  if not bool(errNum): #{
    errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
    if not bool(errNum): #{
      if (dim == 2): #{
        for y in range(0, shape[1]): #{
          offy = y * shape[0]
          for x in range(0, shape[0]): #{
            off = offy + x
            if (gtype == w.WLZ_GREY_UBYTE): #{
              val.ubp[off] = c.c_ubyte(ary[x, y])
            elif (gtype == w.WLZ_GREY_SHORT): #}{
              val.shp[off] = c.c_short(ary[x, y])
            elif (gtype == w.WLZ_GREY_INT): #}{
              val.inp[off] = c.c_int(ary[x, y])
            elif (gtype == w.WLZ_GREY_FLOAT): #}{
              val.flp[off] = c.c_float(ary[x, y])
            elif (gtype == w.WLZ_GREY_DOUBLE): #}{
              val.dbp[off] = c.c_double(ary[x, y])
            #}
          #}
        #}
      else: #}{
        for z in range(0, shape[2]): #{
          rval = obj.contents.values.vox.contents.values[z].r
          val.v =rval.contents.values.v
          for y in range(0, shape[1]): #{
            offy = y * shape[0]
            for x in range(0, shape[0]): #{
              off = offy + x
              if (gtype == w.WLZ_GREY_UBYTE): #{
                val.ubp[off] = c.c_ubyte(ary[x, y, z])
              elif (gtype == w.WLZ_GREY_SHORT): #}{
                val.shp[off] = c.c_short(ary[x, y, z])
              elif (gtype == w.WLZ_GREY_INT): #}{
                val.inp[off] = c.c_int(ary[x, y, z])
              elif (gtype == w.WLZ_GREY_FLOAT): #}{
                val.flp[off] = c.c_float(ary[x, y, z])
              elif (gtype == w.WLZ_GREY_DOUBLE): #}{
                val.dbp[off] = c.c_double(ary[x, y, z])
              #}
            #}
          #}
        #}
      #}
    #}
  #}
  if ((not bool(errNum)) and
      ((atype == np.bool) or (atype == np.bool_))): #{
    errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
    obj1 = w.WlzThresholdI(obj, w.WLZ_THRESH_HIGH, 1, c.byref(errNum))
    if bool(errNum): #{
      w.WlzFreeObj(obj)
      obj = None
    else: #}{
      obj1.contents.values.core = None
      obj.contents.values.core.linkcount = 1
      w.WlzFreeObj(obj)
      obj = obj1
    #}
  #}
  return(obj, errNum)
#}

if __name__ == '__main__': #{
  # Process the command line
  usage = False
  prog = sys.argv[0];
  args = ParseArgs()
  if (args.verbose): #{
    print(prog + ': Reading Woolz object from file ' + args.infile)
  #}
  inobj = ReadWlzObj(args.infile)
  w.WlzAssignObject(inobj, None)
  if (args.verbose): #{
    print(prog + ': WlzFacts inobj')
    PrintWlzFacts(inobj)
  #}
  if (args.verbose): #{
    print(prog + ': Creating a numpy array from the input object.')
  #}
  inary, org, err = ArrayFromWlz(inobj)
  if bool(err): #{
    ErrorExit('Failed to create numpy array from Woolz object (' +
              w.WlzStringFromErrorNum(err, None) + ').')
  #}
  if (args.verbose): #{
    print(prog + ': dtype = ' + str(inary.dtype) + ' shape = ' +
          str(inary.shape)) 
  #}
  if (args.skeletonize): #{
    if (not (inary.dtype == np.bool)) and (not (inary.dtype == np.bool_)): #{
      if (args.verbose): #{
        print(prog + ': Binarizing image using threshold above mean.')
      #}
      thresh = skf.threshold_mean(inary)
      binary = inary > thresh
      inary = binary
    #}
    if (args.verbose): #{
      print(prog + ': Skeletonizing binary array.')
    #}
    if (len(inary.shape) == 2): #{
      outary = skm.skeletonize(inary)
    else: #}{#
      outary = skm.skeletonize_3d(inary)
    #}
  else: #}{
    if (args.verbose): #{
      print(prog + ': No action (output array = input array).')
    #}
    outary = inary
  #}
  if (args.verbose): #{
    print(prog + ': Creating Woolz object from numpy array.')
  #}
  outobj, err = ArrayToWlz(outary, org)
  if bool(err): #{
    ErrorExit('Failed to create Woolz object from numpy array (' +
              w.WlzStringFromErrorNum(err, None) + ').')
  #}
  w.WlzAssignObject(outobj, None)
  if (args.verbose): #{
    print(prog + ': WlzFacts outobj')
    PrintWlzFacts(outobj)
  #}
  if (args.verbose): #{
    print(prog + ': Writing Woolz object to file ' + args.outfile)
  #}
  WriteWlzObj(args.outfile, outobj)
  if (args.verbose): #{
    print(prog + ': Clean exit.')
  #}
  exit(0)
#}


