#!/usr/bin/python
##
# \file         LandmarksXML2Num.py
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

def GetVertex(g): #{
  x = float(g.find('X').text)
  y = float(g.find('Y').text)
  z = float(g.find('Z').text)
  return list((x,y,z))
#}

def ParseArgs(): #{
  parser = argparse.ArgumentParser(description= \
      'Read an XLM landmark or project file from WlzWarp then output the\n' + \
      'landmarks as a .num style landmarks file to the standard output.')
  parser.add_argument('-a', '--absolute', \
      action='store_true', default=False, \
      help='Absolute rather than relative displacements in output')
  parser.add_argument('infile', \
      help='Input WlzWarp file.')
  args = parser.parse_args()
  return(args)
#}

if __name__ == '__main__': #{
  args = ParseArgs()
  doc = et.parse(args.infile)
  lms = doc.find('LandmarkSet')
  for lm in lms.findall('Landmark'):
    src = lm.find('Source')
    tgt = lm.find('Target')
    s = GetVertex(src)
    t = GetVertex(tgt)
    print(str(s[0]) + ' ' + str(s[1]) + ' ' + str(s[2]) + ' ', end='')
    if not args.absolute: #{
      t[0] = t[0] - s[0]
      t[1] = t[1] - s[1]
      t[2] = t[2] - s[2]
    #}
    print(str(t[0]) + ' ' + str(t[1]) + ' ' + str(t[2]))
#}
exit(0)

