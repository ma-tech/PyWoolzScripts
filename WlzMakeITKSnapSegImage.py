#!/usr/bin/python
##
# \file         WlzMakeITKSnapSegImage.py
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
# \brief        Exports Woolz domains as an ITK-SnAP (indexed) segmentation
#               image corresponding to the given ITK-SnAP label description
#               file.
##

from __future__ import print_function
import os
import re
import sys
import ctypes as c
import argparse
import subprocess
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
      'Exports Woolz domains as an ITK-SnAP (indexed) segmentation image ' + \
      'corresponding to the given ITK-SnAP label description file. All ' + \
      'files are both read and written using the Woolz file format, use ' + \
      'WlzExtFFConvert to convert them to NIfTI as required.')
  parser.add_argument('-c', '--compound', \
      action='store_true', default=False, \
      help='output compound object instead of index object.')
  parser.add_argument('-d', '--domaindir', \
      type=str, default='.', \
      help='directory containing the Woolz domains.')
  parser.add_argument('-o', '--output', \
      type=str, default='-', \
      help='output segmentation image.')
  parser.add_argument('-r', '--reference', \
      type=str, default='', \
      help='reference image (used to ensure index object is the same size' + \
           'the reference object).')
  parser.add_argument('-v', '--verbose', \
      action='store_true', default=False, \
      help='verbose output (mainly useful for debugging).')
  parser.add_argument('ldf', help = \
      'input ITK-SnAP label description file in which label text must ' + \
      'correspond to domain file names with a .wlz extension.')
  args = parser.parse_args()
  return(args)
#}

def ExportDomainsToIndexObj(): #{
  # Get domain indices and names from the label description file
  # adding the domain at each step.
  cpd_obj = None
  idx_obj = None
  max_idx = 0
  obj_type = w.WLZ_NULL
  dom_objs = []
  dom_idxs = []
  with open(args.ldf) as ldf: #{
    for rec in ldf: #{
      VerbMsg('rec = ' + rec)
      if not re.match(r'^[ \t]*#', rec): #{
        rec = rec.lstrip()
        dom_idx = int(re.split('\s+', rec)[0])
        dom_name = re.split('"', rec)[1]
        VerbMsg('idx = ' + str(dom_idx) + ' dom = ' + dom_name)
        if(dom_idx > 0): #{
          try: #{
            if dom_idx > max_idx: #{
              max_idx = dom_idx
            #}
            dom_idxs.append(dom_idx)
            errNum = w.enum__WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
            dom_file = args.domaindir + '/' + dom_name + '.wlz'
            fp = libc.fopen(dom_file, 'rb')
            if not bool(fp): #{
              raise IOError()
            #}
            obj = w.WlzAssignObject( \
                  w.WlzReadObj(fp, c.byref(errNum)), None)
            if(bool(errNum)): #{
              raise WlzError()
            #}
            p_lst = w.WlzMakePropertyList(None)
            p_nam = w.WlzMakeNameProperty(dom_name, c.byref(errNum))
            if(bool(errNum)): #{
              raise WlzError()
            #}
            ft = c.CFUNCTYPE(c.c_void_p, c.c_void_p)\
                            (w.WlzFreePropertyListEntry)
            alcerr = w.AlcDLPListEntryAppend(p_lst.contents.list, None, \
                           p_nam, ft)
            if(bool(alcerr)): #{
              raise WlzError()
            #}
            obj.contents.plist = w.WlzAssignPropertyList(p_lst, None)
            dom_objs.append(obj)
            libc.fclose(fp)
            if(bool(errNum)): #{
              raise WlzError()
            #}
            VerbMsg('domain object type (idx == ' + str(dom_idx) + ') ' + \
                    str(w.WlzStringFromObjTypeValue( \
                    obj.contents.type, None)))
            if(obj_type == w.WLZ_NULL): #{
              if((obj.contents.type == w.WLZ_2D_DOMAINOBJ) or \
                 (obj.contents.type == w.WLZ_3D_DOMAINOBJ)): #{
                obj_type = obj.contents.type
              #}
            #}
          except IOError: #}{
            ErrorMsg('Failed to read domain from file ' + dom_file + \
                     ' (' + w.WlzStringFromErrorNum(errNum, None) + ')')
          except WlzError: #}{
            ErrorMsg('Failed to process domain from file ' + dom_file + \
                     ' (' + w.WlzStringFromErrorNum(errNum, None) + ')')
          #}
        #}
      #}
    #}
  #}
  VerbMsg('Object type ' + str(w.WlzStringFromObjTypeValue(obj_type, None)))
  if((max_idx < 1) or
     ((not obj_type == w.WLZ_2D_DOMAINOBJ) and
      (not obj_type == w.WLZ_3D_DOMAINOBJ))): #{
    ErrorMsg('Require at least one 2 or 3D domain with a positive index.')
  #}
  try: #{
    VerbMsg('Creating compound array object of domains.')
    cpd_obj = w.WlzMakeCompoundArray( \
                w.enum__WlzObjectType(w.WLZ_COMPOUND_ARR_2), \
                1, max_idx + 1, None, obj_type, c.byref(errNum))
    if(bool(errNum)): #{
      raise WlzError()
    #}
    for i in range(0, len(dom_objs)): #{
      obj = dom_objs[i]
      idx = dom_idxs[i]
      cpd_obj.contents.o[idx] = obj
    #}
  except: #}{
    ErrorMsg('Failed to create compound object from domains (' + \
             w.WlzStringFromErrorNum(errNum, None) + ')')
  #}
  if not args.compound: #{
    try: #{
      VerbMsg('Creating index object from compound of domains.')
      idx_obj = w.WlzAssignObject( \
                w.WlzIndexObjFromCompound(cpd_obj, c.byref(errNum)), None)
      if(bool(errNum)): #{
        raise WlzError()
      #}
    except: #}{
      ErrorMsg('Failed to create index object from compound of domains (' + \
               w.WlzStringFromErrorNum(errNum, None) + ')')
    #}
  #}
  if bool(args.reference) and (not args.compound): #{
    try: #{
      VerbMsg('Cutting index object to the same bounding box as the\n' + \
              'reference object (read from ' + args.reference + ').')
      fp = None
      errNum = w.enum__WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
      if args.reference == '-': #{
        fp = libc.stdin
      else: #}{
        fp = libc.fopen(args.reference, 'rb')
      #}
      if not bool(fp): #{
        raise IOError()
      #}
      obj = w.WlzAssignObject( \
            w.WlzReadObj(fp, c.byref(errNum)), None)
      if not args.reference == '-': #{
        libc.fclose(fp)
      #}
      if(bool(errNum)): #{
        raise WlzError()
      #}
      box = w.WlzBoundingBox3I(obj, c.byref(errNum))
      if(not bool(errNum)): #{
        gtype = w.WlzGreyTypeFromObj(idx_obj, c.byref(errNum))
      #}
      w.WlzFreeObj(obj)
      if(bool(errNum)): #{
        raise WlzError()
      #}
      VerbMsg('Using bounding box ' + \
              '(' + str(box.xMin) + ',' + str(box.xMax) + '),' + \
              '(' + str(box.yMin) + ',' + str(box.yMax) + '),' + \
              '(' + str(box.zMin) + ',' + str(box.zMax) + ')')
      VerbMsg('Preserving grey type ' + \
              w.WlzStringFromGreyType(gtype, None) + \
              '.')
      obj = w.WlzCutObjToBox3D(idx_obj, box, gtype, 0, 0.0, 0.0, \
                                 c.byref(errNum))
      if(bool(errNum)): #{
        raise WlzError()
      #}
      w.WlzFreeObj(idx_obj)
      idx_obj = obj
    except IOError: #}{
      ErrorMsg('Failed to read the reference object read from ' + \
               args.reference + '.')
    except WlzError: #}{
      ErrorMsg('Failed to cut index object to the bounding box of the\n' + \
               'reference object read from (' + args.reference + ') (' + \
               w.WlzStringFromErrorNum(errNum, None) + ')')
    #}
  #}
  try: #{
    ft = 'index'
    if args.compound: #{
      ft = 'compound'
    #}
    VerbMsg('Writing ' + ft + ' object to output file (' + args.output + ').')
    errNum = w.enum__WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
    if args.output == '-': #{
      fp = libc.stdout
    else: #}{
      fp = libc.fopen(args.output, 'wb')
    #}
    if not bool(fp): #{
      raise IOError()
    #}
    if args.compound: #{
      errNum = w.WlzWriteObj(fp, cpd_obj)
    else: #}{
      errNum = w.WlzWriteObj(fp, idx_obj)
    #}
    if not args.output == '-': #{
      libc.fclose(fp)
    #}
    if(bool(errNum)): #{
      raise WlzError()
    #}
  except: #}{
    ErrorMsg('Failed to write ' + ft + ' object to file (' + \
             w.WlzStringFromErrorNum(errNum, None) + ')')
  #}
#}

if __name__ == '__main__': #{
  # Process the command line
  args = ParseArgs()
  prog = sys.argv[0];
  errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
  if(args.verbose): #{
    print(prog + ': args = ' + str(args))
  #}
  ExportDomainsToIndexObj()
#}

