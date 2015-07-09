#!/usr/bin/python
##
# \file         MAExportAnatomyToNIfTI.py
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
# \brief        Creates NIfTI format files containing EMAP anatomy
#               as indexed image files along with ITK-SnAP label
#               description files.
##

from __future__ import print_function
import os
import re
import sys
import argparse
import subprocess
import tempfile
import ctypes
import urllib2
import json

class Label(object):
  def __init__(self, name, accession, color, index):
    self.name = name
    self.accession = accession
    self.color = color
    self.index = index
  def __repr__(self):
    return "<name: %s, accession: %s, color: %d,%d,%d>, index %d" % \
           (self.name, self.accession, \
            self.color[0], self.color[1], self.color[2], \
            self.index)
  
singleStage        = None
emapSrvDescURL     = \
    'http://www.emouseatlas.org/jsonservice/json_service?service=range'
emapAnatomyDescURL = \
    'http://www.emouseatlas.org/jsonservice/json_service?service=domain'
WlzExtFFConvert    = '/opt/MouseAtlas/bin/WlzExtFFConvert'

def ParseArgs():
  parser = argparse.ArgumentParser(description= \
  'Creates NIfTI format files containing EMAP anatomy as indexed image\n' +
  'files along with ITK-SnAP label description files.')
  parser.add_argument('-d', '--destdir', \
      type=str, default='.', \
      help='Destination directory for the output files.')
  parser.add_argument('-n', '--nofiles', \
      action='store_true', default=False, \
      help='Dont write any files (label files written stdout, mainly\n' + \
           'useful for debugging).')
  parser.add_argument('-i', '--noindex', \
      action='store_true', default=False, \
      help='Dont make a NIfTI index file.')
  parser.add_argument('-s', '--stage', \
      type=str, default=None, \
      help='Single stage (eg TS14).')
  parser.add_argument('-v', '--verbose', \
      action='store_true', default=False, \
      help='Verbose output (mainly useful for debugging).')
  args = parser.parse_args()
  return(args)

def CleanExit(stat):
  exit(stat)

def MakeFlatLabelTable(table, entry): #{
  name = entry['name']
  if('children' in entry): #{
    for ent in entry['children']:
      MakeFlatLabelTable(table, ent)
  else: #}{
    row = Label(name, entry['domain']['accession'], entry['domain']['color'],
                entry['domain']['index'])
    table.append(row)
  #}
  return(table)
#}

def CommandLineString(args): #{
  rtn = ''
  sep = ''
  for a in args: #{
    rtn = rtn + sep + str(a)
    sep = ' '
  #}
  return(rtn)
#}


if __name__ == '__main__':
  # Process the command line
  args = ParseArgs()
  prog = sys.argv[0];
  if(args.stage): #{
    match = re.match(r'^TS\d\d?$', args.stage)
    if(match): #{
      singleStage = match.group(0)
    else: #}{
      print(prog + ': Bad stage, specify using TS[0-9][0-9]? (eg TS07\n' +
            'or TS17).')
      CleanExit(1)
    #}
  #}

  if(args.verbose): #{
    print(prog + ': args = ' + str(args))
  #}

  # Fetch the EMAP service description file
  if(args.verbose): #{
    print(prog + ': Fetching service description file (using url =\n' + \
          emapSrvDescURL + ').')
  #}
  try: #{
    url = urllib2.urlopen(emapSrvDescURL)
  except urllib2.HTTPError as e: #}{
    print(prog + ': Failed to open service description file (url =\n' + \
          emapSrvDescURL + ').') 
    CleanExit(1)
  #}
  srvDesc = json.load(url)
  #print(json.dumps(srvDesc))
  # For each entry of the EMAP service description file
  for entry in srvDesc: #{
    idx = 0
    stage = entry['stage']
    if(stage and ((not singleStage) or (singleStage == stage))): #{
      idxFile = None
      if(args.verbose): #{
        print(prog + ': Processing stage ' + stage + '.')
      #}
      if 'voxelModel' in entry: #{
        voxMod = entry['voxelModel']
        # There may be many models ie voxMod[i] with i > 0
        # use  model = "indexEMA118" stripped of the index prefix, eg:
        # http://www.emouseatlas.org/jsonservice/json_service? \
        # service=domain&stage=TS17&model=EMA118
        for mod in voxMod: #{
          for mfm in mod: #{
            model = None
            if(args.verbose): #{
              print(prog + ': Checking for model, name = ' + \
                    mfm['name'] + '.')
            #}
            # Sometimes there is not EMA model just default, so allow
            # default here
            match = re.match(r'[a-z]*((EMA\d+)|(default))', mfm['name'])
            if(match): #{
              model = match.group(1)
            #}
            if(model and (mfm['type'] == 'index')): #{
              idxFile = mfm['fileName']
              if(args.verbose): #{
                print(prog + ': At stage ' + stage + \
                      ' with model ' + model + ', idxFile = ' + \
                      str(idxFile) + '.')
              #}
              # If the index file exists  convert it to NIfTI then create the
              # label description file
              if(idxFile): #{
                filebase = stage + '_' + model + '_anatomy'
                # Convert the Woolz index file to NIfTI
                if((not args.nofiles) and (not args.noindex)): #{
                  cmdline = [WlzExtFFConvert, '-o' + filebase + '.nii', \
                             idxFile]
                  if(args.verbose): #{
                    print(prog + \
                          ': Converting index object to NIfTI format\n' + \
                          'using: ' + CommandLineString(cmdline))
                  #}
                  rtn = subprocess.call(cmdline)
                  if(bool(rtn)): #{
                    print(prog + \
                          ': WlzExtFFConvert failed to convert index\n' + \
                          'from Woolz to NIfTI format.')
                  #}
                #}
                # Fetch the anatomy description file and then make a ITK SnAP
                # label file
                adu = emapAnatomyDescURL + '&stage=' + stage + \
                      '&model=' + model
                if(args.verbose):
                  print(prog + ': Fetching anatomy description file (using\n' + \
                        adu + ').') 
                try:
                  url = urllib2.urlopen(adu)
                except urllib2.HTTPError as e:
                  print(prog + ': Failed to open anatomy description file\n' +
                        '(url' + adu + ').') 
                  CleanExit(1)
                anatomyDesc = json.load(url)
                #print(json.dumps(anatomyDesc))
                if(args.nofiles): #{
                  f = sys.stdout
                else: #}{
                  outfile = filebase + '.txt'
                  f = open(outfile, 'w')
                #}
                table = MakeFlatLabelTable([], anatomyDesc)
                # Sort the label table by index
                table = sorted(table, key=lambda label: label.index)
                # Print label table to file
                print('################################################\n' + \
                      '# ITK-SnAP Label Description File\n' + \
                      '# File format: \n' + \
                      '# IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL\n' + \
                      '# Fields: \n' + \
                      '#    IDX:   Zero-based index \n' + \
                      '#    -R-:   Red color component (0..255)\n' + \
                      '#    -G-:   Green color component (0..255)\n' + \
                      '#    -B-:   Blue color component (0..255)\n' + \
                      '#    -A-:   Label transparency (0.00 .. 1.00)\n' + \
                      '#    VIS:   Label visibility (0 or 1)\n' + \
                      '#    IDX:   Label mesh visibility (0 or 1)\n' + \
                      '#  LABEL:   Label description \n' + \
                      '################################################\n' + \
                      '     0    0    0    0 0.0      0 0 "Clear Label"',
                      file=f)
                for row in table: #{
                  print('% 6d' % row.index + ' ', end='', file=f)
                  print('% 4d' % row.color[0] + ' ', end='', file=f)
                  print('% 4d' % row.color[1] + ' ', end='', file=f)
                  print('% 4d' % row.color[2] + ' ', end='', file=f)
                  print('1.0      1 1 ' + \
                        '"' + row.accession + ' ' + row.name + '"', file=f)
                #}
                if(not args.nofiles):
                  f.close()
              #}
            #}
          #}
        #}
      #}
    #}
  #}
        

  if(args.verbose):
    print(prog + ': Cleaning up.')
  CleanExit(0)
