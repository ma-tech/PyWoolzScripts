#!/usr/bin/python
##
# \file         LandmarksNum2XML.py
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
# \brief        Puts simple white space seperated tie points into a
#               WlzWarp XML landmarks file
##

from __future__ import print_function

import os
import re
import sys
import argparse

def ParseArgs(): #{
  parser = argparse.ArgumentParser(description= \
      'Reads a white space seperated (.num style) tie points file and\n' + \
      'then writes it as a WlzWarp XLM landmark to the standard output.')
  parser.add_argument('-a', '--absolute', \
      action='store_true', default=False, \
      help='Absolute rather than relative displacements in input')
  parser.add_argument('infile', \
      help='Input tie points file.')
  args = parser.parse_args()
  return(args)
#}

if __name__ == '__main__': #{
  cnt = 0
  args = ParseArgs()
  prog = sys.argv[0];
  rex = re.compile('[ \t]')
  print('<?xml version="1.0" encoding="UTF-8"?>')
  print('<LandmarkSet>')
  print('  <Type>3D</Type>')
  if not args.infile == '-': #{
    f = open(args.infile)
  else: #}{
    f = sys.stdin
  #}
  for ln in f: #{
    cnt = cnt + 1
    rec = rex.split(ln)
    if not len(rec) == 6: #{
      print(prog + ': Invalid input at line ' + str(cnt))
      exit(1)
    #}
    lmk = [0.0] * 6
    for idx in range(0,6): #{
      lmk[idx] = float(rec[idx])
    #}
    if not args.absolute: #{
      lmk[3] = lmk[0] + lmk[3];
      lmk[4] = lmk[1] + lmk[4];
      lmk[5] = lmk[2] + lmk[5];
    #}
    print('<Landmark>')
    print('  <Source>')
    print('    <X>' + str(lmk[0]) + '</X>')
    print('    <Y>' + str(lmk[1]) + '</Y>')
    print('    <Z>' + str(lmk[2]) + '</Z>')
    print('  </Source>')
    print('  <Target>')
    print('    <X>' + str(lmk[3]) + '</X>')
    print('    <Y>' + str(lmk[4]) + '</Y>')
    print('    <Z>' + str(lmk[5]) + '</Z>')
    print('  </Target>')
    print('</Landmark>')
  #}
  print('</LandmarkSet>')
  if f is not sys.stdin: #{
    f.close()
  #}
#}
exit(0)

