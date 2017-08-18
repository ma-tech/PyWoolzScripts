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
import json
import glob
import errno
import urllib2
import tempfile
import argparse
import subprocess
import ctypes as c
import Wlz as w

libc = c.CDLL('libc.so.6')

class WlzError(Exception):
  pass

class Label(object): #{
  def __init__(self, name, accession, color, index):
    self.name = name
    self.accession = accession
    self.color = color
    self.index = index
  def __repr__(self): #{
    return '<name: %s, accession: %d, color: %d,%d,%d>, index %d' % \
           (self.emapa, self.name, self.accession, \
            self.color[0], self.color[1], self.color[2], \
            self.index)
  #}
#}

prog = ''
args = None

origin           = 'http://www.emouseatlas.org'

db_image_area    = '/opt/emageDBLocation/dbImageArea/atlas/voxel3D/'

wlzbin           = '/opt/MouseAtlas/bin/'
wlzextffconvert  = wlzbin + 'WlzExtFFConvert'

single_stage       = None
emap_srv_desc_url  = \
    'http://www.emouseatlas.org/jsonservice/json_service?service=range'
emap_anat_desc_url = \
    'http://www.emouseatlas.org/jsonservice/json_service?service=domain'
emap_atlasviewer_url = \
    'http://www.emouseatlas.org/eAtlasViewer_ema/application/ema/anatomy/'
emap_atlasviewer_tree = \
    '/tree/treeData.jso'
emap_atlasviewer_tiledimagemodel = 'tiledImageModelData.jso'

def ChkWoolzError(errNum, status): #{
  if(bool(errNum)): #{
    msg = 'Woolz error ' + w.WlzStringFromErrorNum(errNum, None)
    if(status < 1): #{
      MsgWarn(msg)
    else: #}{
      MsgError(msg, status)
    #}
    raise WlzError()
  #}
#}

def CleanExit(stat): #{
  MsgVerbose('Exit with status ' + str(stat))
  exit(stat)
#}

def MsgVerbose(msg): #{
  if(args.verbose): #{
    print(prog + ': ' + msg, file=sys.stderr)
  #}
#}

def MsgWarn(msg): #{
  print(prog + ': WARNING - ' + msg, file=sys.stderr)
#}

def MsgError(msg, stat): #{
  print(prog + ': ERROR - ' + msg, file=sys.stderr)
  if(stat > 0): #{
    CleanExit(stat)
  #}
#}

def ReadJsonStrFromURL(u):
  o = urllib2.urlopen(u)
  s = o.read()
  # Remove any trailing commas!
  s = re.sub(",[ \t\r\n]+}", "}", s)
  s = re.sub(",[ \t\r\n]+\]", "]", s)
  return(s)

def ConvertFileWlzToNii(wlz): #{
  # Just can't see why this doesn't work! Fine if nii replaced by tif
  # if run in a shell or even if run in python from REPL.
  stat = 0
  nii = re.sub('\.wlz$', '.nii', wlz)
  cmd = wlzextffconvert + ' -f wlz -F nii -o ' + nii + ' ' + wlz
  MsgVerbose(cmd)
  stat = subprocess.call(cmd, shell=True)
  if(stat > 0): #{
    MsgWarn('Command Failed (' + cmd + ')(' + str(stat) + ')')
  #}
  return(stat)
#}

def ParseArgs(): #{
  parser = argparse.ArgumentParser(description= \
  'Creates NIfTI format files containing EMAP anatomy as indexed image\n' +
  'files along with ITK-SnAP label description files.')
  parser.add_argument('-a', '--atlasviewer', \
      action='store_true', default=False, \
      help='Add information from the eAtlasViewer:\n'
           'Overwrite the EMAP service colors with those of the \n' + \
           'appropriate eAtlasViewer colors if appropriate and the \n' + \
           'eAtlasViewer data exists.\n' +
           'Set voxel sizes using information from the eAtlasViewer.')
  parser.add_argument('-d', '--destdir', \
      type=str, default='.', \
      help='Destination directory for the output files.')
  parser.add_argument('-i', '--noindex', \
      action='store_true', default=False, \
      help='Don\'t make index files.')
  parser.add_argument('-n', '--nofiles', \
      action='store_true', default=False, \
      help='Don\'t write any files (label files written stdout, mainly\n' + \
           'useful for debugging).')
  parser.add_argument('-s', '--stage', \
      type=str, default=None, \
      help='Single stage (eg TS14).')
  parser.add_argument('-t', '--directory-tree',
      action='store_true', default=False, \
      help='Place files in a stage/model/file directory tree.')
  parser.add_argument('-v', '--verbose', \
      action='store_true', default=False, \
      help='Verbose output (mainly useful for debugging).')
  args = parser.parse_args()
  return(args)
#}

def MakeFlatLabelTable(table, entry, model): #{
  name = entry['name']
  if('domain' in entry): #{
    try: #{
      row = Label(name, int(entry['domain']['accession'].split(':')[1]),
                  entry['domain']['color'], entry['domain']['index'])
      table.append(row)
    except: #}{
      MsgWarn('parse error for model ' + model + ' at ' + name)
    #}
  #}
  if('children' in entry): #{
    for ent in entry['children']: #{
      MakeFlatLabelTable(table, ent, model)
    #}
  #}
  return(table)
#}

def MergeTreeColors(table, anat_tree): #{
  for row in table: #{
    a_node = None
    t_EMAPA = row.accession
    t_name = row.name
    # Search for the matching EMAPA number in the anatomy tree
    for node in anat_tree: #{
      n_EMAPA = -1
      for eid in node['extId']: #{
        e = re.split('[^0-9a-zA-Z]', eid)
        try: #{
          if(e[0] == 'EMAPA'): #{
            n_EMAPA = int(e[1])
          #}
          if(n_EMAPA == t_EMAPA): #{
            a_node = node
            break
          #}
        except: #}{
          pass
        #}
      #}
      if(bool(a_node)): #{
        break
      #}
    #}
    # If EMAPA number not found search for exact name match
    if(not bool(a_node)): #{
      for node in anat_tree: #{
        n_name = node['name']
        if(t_name == n_name): #{
          a_node = node
          break
        #}
      #}
    #}
    if(bool(a_node)): #{
      try: #{
        a_color = a_node['domainData']['domainColour']
        MsgVerbose('Merged |' + str(t_EMAPA) + '|' +
                   row.name + '|' + a_node['name'] + '|' +
                   str(row.color) + '|' + str(a_color))
        row.color = [int(a_color[0]), int(a_color[1]), int(a_color[2])]
      except: #}{
        a_node = False
      #}
    #}
    if(not bool(a_node)): #{
      MsgWarn('failed to merge color for EMAPA:' + str(t_EMAPA) + ' ' +
              row.name)
    #}
  #}
#}

if __name__ == '__main__': #{
  #
  # Process the command line
  args = ParseArgs()
  prog = sys.argv[0];
  errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
  if(args.stage): #{
    match = re.match(r'^TS\d\d?$', args.stage.upper())
    if(match): #{
      single_stage = match.group(0)
    else: #}{
      MsgError('Bad stage, specify using TS[0-9][0-9]? (eg TS07\n' +
               'or TS17).', 1)
    #}
  #}
  MsgVerbose('args = ' + str(args))
  #
  # Fetch the EMAP service description file
  MsgVerbose('Fetching service description file (using url =\n' +
             emap_srv_desc_url + ').')
  try: #{
    s = ReadJsonStrFromURL(emap_srv_desc_url)
    srvDesc = json.loads(s)
  except: #}{
    MsgError('Failed to open service description file (url =\n' +
             emap_srv_desc_url + ').')
  #}
  #
  # For each entry of the EMAP service description file
  for entry in srvDesc: #{ service description entry
    idx = 0
    stage = entry['stage']
    if(stage and ((not single_stage) or (single_stage == stage))): #{ stage
      MsgVerbose('Processing stage ' + stage + '.')
      if 'voxelModel' in entry: #{
        vox_model = entry['voxelModel']
        # There may be many models ie vox_model[i] with i > 0. Here we
        # select the grey value object for each model in turn
        for mod in vox_model: #{ model
          vox_file = None
          model = None
          model_has_anat = True
          out_file_idx = None
          out_file_ref = None
          out_file_txt = None
          anat_table = None
          idx_obj = None
          idx_max = 0
          ref_obj = None
          voxel_sz = [1.0, 1.0, 1.0]
          out_base_dir = args.destdir
          for mfm in mod: #{
            MsgVerbose('Checking for model, name = ' + mfm['name'] + '.')
            if(mfm['type'] == 'grey'): #{
              match = re.match(r'(EMA\d+)', mfm['name'])
              if(match): #{
                model = match.group(1)
                vox_file = str(mfm['fileName'])
                MsgVerbose('At stage ' + stage + ' with model ' +
                           model + ', vox_file = ' + vox_file + '.')
              else: #}{
                MsgWarn('parse error for model in ' + mfm['name'])
              #}
            #}
          #}
          if(args.directory_tree): #{
            out_base_dir = out_base_dir + '/' + str(stage) + '/' + str(model)
            try: #{
              os.makedirs(out_base_dir)
            except OSError as ex: #}{
              if(ex.errno == errno.EEXIST): #{
                pass
              else: #}{
                MsgError('Failed to create output directory ' + out_base_dir,
                         1)
              #}
            #}
          #}
          file_base = out_base_dir + '/' + stage + '_' + model
          # Fetch the anatomy description file and then make a ITK SnAP
          # label file
          adu = emap_anat_desc_url + '&stage=' + stage + \
                '&model=' + model
          MsgVerbose('Fetching anatomy description file (using ' + adu + ').')
          try: #{
            s = ReadJsonStrFromURL(adu)
            anatomy_desc = json.loads(s)
          except: #}{
            model_has_anat = False
            MsgWarn('Assuming no anatomy since failed to load anatomy ' + 
                    'description file\n (url ' + adu + ').')
          #}
          if(model_has_anat): #{ model_has_anat
            #print(str(anatomy_desc))
            # Make a flat anatomy table from the anatomy description
            anat_table = MakeFlatLabelTable([], anatomy_desc, model)
            # Optionally overwrite the flat anatomy table colours
            # using atlasViewer tree colours.
            anat_tree = None
            if(args.atlasviewer): #{
              treeURL = emap_atlasviewer_url + model + emap_atlasviewer_tree
              MsgVerbose('Overwriting anatomy colors using file (url ' +
                         treeURL + ').')
              try: #{
                s = ReadJsonStrFromURL(treeURL)
                anat_tree = json.loads(s)
              except: #}{
                anat_tree = None
              #}
              if(bool(anat_tree)): #{
                #print(json.dumps(anat_tree))
                MergeTreeColors(anat_table, anat_tree)
              #}
            #}
            # Sort the label table by index
            anat_table = sorted(anat_table,
                                key=lambda label: label.index)
            # Build label descriptions
            label_desc = \
                '################################################\n' + \
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
                '################################################\n'
            label_desc = label_desc + \
                '      0    0    0    0 0.0      0 0 "Clear Label"\n'
            for row in anat_table: #{
              if(row.index > idx_max): #{
                idx_max = row.index
              #}
              label_desc = label_desc + \
                  ' % 6d' % row.index + \
                  ' % 4d' % row.color[0] + \
                  ' % 4d' % row.color[1] + \
                  ' % 4d' % row.color[2] + \
                  ' 1.0      1 1' + \
                  ' "EMAPA:' + str(row.accession) + ' ' + row.name + '"\n'
            #}
            # Print anatomy label table to file
            try: #{
              if(args.nofiles or args.noindex): #{
                f = sys.stdout
              else: #}{
                out_file_txt = file_base + '_anatomy.txt'
                f = open(out_file_txt, 'w')
              #}
              print(label_desc, end='', file=f)
              if(not (args.nofiles or args.noindex)): #{
                f.close()
              #}
            except: #}{
              MsgError('Failed to write label description to file ' +
                       out_file_txt + '.', 1)
            #}
          #} model_has_anat
          #
          # Get voxel size from eAtlasViewer
          if(args.atlasviewer): #{
            tim_desc = ''
            timURL = emap_atlasviewer_url + model + '/' + \
                     emap_atlasviewer_tiledimagemodel
            MsgVerbose('Getting true voxel size using file (url ' +
                       timURL + ').')
            try: #{
              s = ReadJsonStrFromURL(timURL)
              tim_desc = json.loads(s)
              v_sz = tim_desc['voxelSize']
              voxel_sz[0] = float(v_sz['x'])
              voxel_sz[1] = float(v_sz['y'])
              voxel_sz[2] = float(v_sz['z'])
            except: #}{
              voxel_sz = [1.0, 1.0, 1.0]
              MsgWarn('Assuming voxel size (1.0, 1.0, 1.0) since failed ' +
                      'to load tiled image model ' +
                      'from file\n (url ' + timURL + ').')
            #}
            MsgVerbose('Voxel size is ' + str(voxel_sz) + '.')
          #}
          #
          # Copy the grey scale Woolz reference object setting name and
          # text properties
          if(not args.nofiles): #{
            fple = c.CFUNCTYPE(c.c_void_p, c.c_void_p) \
                   (w.WlzFreePropertyListEntry)
            # Read the Woolz grey scale reference object and write it
            try: #{
              MsgVerbose('Reading voxel image from file ' + vox_file)
              errNum = w.enum__WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
              fp = libc.fopen(vox_file.encode('utf-8'), 'rb')
              ref_obj = w.WlzAssignObject( \
                        w.WlzReadObj(fp, c.byref(errNum)), None)
              libc.fclose(fp)
              ChkWoolzError(errNum, 0)
              MsgVerbose('Setting voxel image voxel size to ' + str(voxel_sz))
              w.WlzSetVoxelSize(ref_obj, voxel_sz[0], voxel_sz[1], voxel_sz[2])
              MsgVerbose('Creating reference object property list.')
              p_lst = w.WlzMakePropertyList(None)
              p_str = stage + '_' + model  + '_reference'
              p_nam = w.WlzMakeNameProperty(p_str.encode('utf-8'),
                                            c.byref(errNum))
              p_txt = w.WlzMakeTextProperty('Origin'.encode('utf-8'),
                                            origin.encode('utf-8'),
                                            c.byref(errNum))
              alcerr = w.AlcDLPListEntryAppend(p_lst.contents.list, None,
                                               p_nam, fple)
              alcerr = w.AlcDLPListEntryAppend(p_lst.contents.list, None,
                                               p_txt, fple)
              if(bool(alcerr)): #{
                errNum = w.enum__WlzErrorNum(w.WLZ_ERR_MEM_ALLOC)
              #}
              ChkWoolzError(errNum, 0)
              ref_obj.contents.plist = w.WlzAssignPropertyList(p_lst, None)
              out_file_ref = file_base + '_reference.wlz'
              MsgVerbose('Writing voxel image to file ' + out_file_ref)
              fp = libc.fopen(out_file_ref.encode('utf-8'), 'wb')
              errNum = w.WlzWriteObj(fp, ref_obj)
              libc.fclose(fp)
              ChkWoolzError(errNum, 0)
            except: #}{
              MsgError('Failed to copy grey scale Woolz voxel object (' +
                       vox_file + ')', 1)
            #}
          #}
          if(model_has_anat): #{
            #
            # Create compound object from the anatomy domains.
            if(not (args.nofiles or args.noindex)): #{
              try: #{
                nullValues = w.WlzValues(None)
                errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
                cpd_obj = w.WlzMakeCompoundArray( \
                    w.enum__WlzObjectType(w.WLZ_COMPOUND_ARR_2), \
                    1, idx_max + 1, None,
                    ref_obj.contents.type, c.byref(errNum)) 
                ChkWoolzError(errNum, 0)
              except: #}{
                MsgError('Failed to create compound object for ' +
                         stage + ' ' + model +
                         ' (' + w.WlzStringFromErrorNum(errNum, None) + ')', 1)
              #}
              for idx in range(0, len(anat_table)): #{
                row = anat_table[idx]
                domfile = glob.glob(db_image_area + stage.lower() +'/' + 
                                    model +'/' + 
                                    '/anatomy/EMAPA_' + str(row.accession) +
                                    '/*.wlz')
                if(len(domfile) == 0): #{
                  domfile = None
                  MsgWarn('Anatomy domain file not found for ' + stage +
                          ' ' + model + ' EMAPA_' + str(row.accession))
                elif(len(domfile) > 1): #}{
                  MsgError('Multiple domain files found for ' + stage +
                           ' ' + model + ' EMAPA_' + str(row.accession), 1)
                else: #}{
                  domfile = domfile[0]
                #}
                try: #{
                  MsgVerbose('Reading anatomy domain from ' + domfile)
                  errNum = w.enum__WlzErrorNum(w.WLZ_ERR_FILE_OPEN)
                  fp = libc.fopen(domfile.encode('utf-8'), 'rb')
                  cpd_obj.contents.o[row.index] = w.WlzAssignObject( \
                      w.WlzReadObj(fp, c.byref(errNum)), None)
                  libc.fclose(fp)
                  ChkWoolzError(errNum, 0)
                except: #}{
                  MsgError('Failed to read anatomy domain from' + domfile, 1)
                #}
              #}
              try: #{
                alcerr = [None, None, None]
                MsgVerbose('Creating index object from compound of domains.')
                errNum = w.enum__WlzErrorNum(w.WLZ_ERR_NONE)
                tmp_obj = w.WlzAssignObject(
                          w.WlzIndexObjFromCompound(cpd_obj, 
                                                    c.byref(errNum)), None)
                w.WlzFreeObj(cpd_obj)
                ChkWoolzError(errNum, 0)
                MsgVerbose('Cutting index image to a cuboid.')
                box = w.WlzBoundingBox3I(ref_obj, c.byref(errNum))
                ChkWoolzError(errNum, 0)
                idx_obj_grey_type = w.WLZ_GREY_INT
                if (idx_max < 255): #{
                  idx_obj_grey_type = w.WLZ_GREY_UBYTE
                elif (idx_max < 32768): #}{
                  idx_obj_grey_type = w.WLZ_GREY_SHORT
                #}
                idx_obj = w.WlzAssignObject(
                          w.WlzCutObjToValBox3D(tmp_obj, box,
                             w.enum__WlzGreyType(idx_obj_grey_type), None,
                             0, 0.0, 0.0, c.byref(errNum)), None)
                ChkWoolzError(errNum, 0)
                w.WlzFreeObj(tmp_obj)
                MsgVerbose('Setting index image voxel size to ' + str(voxel_sz))
                w.WlzSetVoxelSize(idx_obj, 
                                  voxel_sz[0], voxel_sz[1], voxel_sz[2])
                MsgVerbose('Creating index object property list.')
                p_lst = w.WlzMakePropertyList(None)
                p_str = stage + '_' + model  + '_anatomy'
                p_nam = w.WlzMakeNameProperty(p_str.encode('utf-8'),
                                              c.byref(errNum))
                ChkWoolzError(errNum, 0)
                p_org = w.WlzMakeTextProperty('Origin'.encode('utf-8'),
                                              origin.encode('utf-8'),
                                              c.byref(errNum))
                ChkWoolzError(errNum, 0)
                p_lab = w.WlzMakeTextProperty(
                            'Label Decriptions'.encode('utf-8'),
                            label_desc.encode('utf-8'),
                            c.byref(errNum))
                ChkWoolzError(errNum, 0)
                alcerr[0] = w.AlcDLPListEntryAppend(p_lst.contents.list, None,
                                                    p_nam, fple)
                alcerr[1] = w.AlcDLPListEntryAppend(p_lst.contents.list, None,
                                                    p_org, fple)
                alcerr[2] = w.AlcDLPListEntryAppend(p_lst.contents.list, None,
                                                    p_lab, fple)
                if(bool(alcerr[0]) or bool(alcerr[1]) or bool(alcerr[2])): #{
                  errNum = w.enum__WlzErrorNum(w.WLZ_ERR_MEM_ALLOC)
                #}
                ChkWoolzError(errNum, 0)
                idx_obj.contents.plist = w.WlzAssignPropertyList(p_lst, None)
              except Exception as e: #}{
                MsgError('Failed to create index object from compound object ' +
                         stage + model + ' (' + str(e) + ') ' + str(alcerr) +
                         ' (' + w.WlzStringFromErrorNum(errNum, None) + ')', 1)
              #}
            #}
          #} model_has_anat
          # Write index object to file
          if(model_has_anat and not (args.nofiles or args.noindex)): #{
            out_file_idx = file_base + '_anatomy.wlz'
            try: #{
              MsgVerbose('Writing index image to file ' + out_file_idx)
              fp = libc.fopen(out_file_idx.encode('utf-8'), 'wb')
              errNum = w.WlzWriteObj(fp, idx_obj)
              libc.fclose(fp)
              ChkWoolzError(errNum, 0)
            except: #}{
              MsgError('Failed to write anatomy index object to file' +
                       out_file_idx +
                       ' (' + w.WlzStringFromErrorNum(errNum, None) + ')', 1)
            #}
          #}
          if(not args.nofiles): #{
            if(model_has_anat and not (args.noindex)): #{
              ConvertFileWlzToNii(out_file_idx)
            #}
            ConvertFileWlzToNii(out_file_ref)
          #}
          w.WlzFreeObj(idx_obj)
          w.WlzFreeObj(ref_obj)
        #} model
      #}
    #} stage
  #} service description entry
  CleanExit(0)
#}
