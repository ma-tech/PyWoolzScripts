#!/usr/bin/python
##
# \file         WlzWarpPrjWarpValues.py
# \author       Bill Hill
# \date         October 2015
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
# \brief        Extracts landmarks from a WlzWarp file
##

from __future__ import print_function

import os
import sys
import argparse
from xml.etree import ElementTree as et

def ParseArgs(): #{
  parser = argparse.ArgumentParser(description= \
      'Reads a WlzWarp project file and outputs the warp values to the\n' + \
      'standard output. Default action s to retrieve all values, but if\n' + \
      'any are specified then only those will be retrieve.')
  parser.add_argument('-b', '--basis', \
      action='store_true', default=False, \
      help='Retrieve warp basis function type.')
  parser.add_argument('-d', '--delta', \
      action='store_true', default=False, \
      help='Retrieve warp delta value.')
  parser.add_argument('-s', '--snap', \
      action='store_true', default=False, \
      help='Retrieve warp snap fit distance.')
  parser.add_argument('-v', '--verbose', \
      action='store_true', default=False, \
      help='Verbose output.')
  parser.add_argument('infile', \
      help='WlzWarp project input file.')
  args = parser.parse_args()
  return(args)
#}

if __name__ == '__main__': #{
  args = ParseArgs()
  if (not args.basis) and (not args.delta) and (not args.snap): #{
    args.basis = True
    args.delta = True
    args.snap = True
  #}
  if not args.infile == '-': #{
    f = open(args.infile)
  else: #}{
    f = sys.stdin
  #}
  doc = et.parse(f)
  wrp = doc.find('Warping')
  if(args.basis): #{
    basis = wrp.find('BasisFunctionType').text
    if(args.verbose): #{
      print('basis ', end='')
    #}
    print(basis)
  #}
  if(args.delta): #{
    delta = wrp.find('Delta').text
    if(args.verbose): #{
      print('delta ', end='')
    #}
    print(delta)
  #}
  if(args.snap): #{
    snap = wrp.find('SnapToFitDist').text
    if(args.verbose): #{
      print('snap ', end='')
    #}
    print(snap)
  #}
  if f is not sys.stdin: #{
    f.close()
  #}
#}
exit(0)

