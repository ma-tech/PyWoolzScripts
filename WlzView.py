#!/usr/bin/python

from __future__ import print_function
import os
import sys
import argparse
import logging
import ctypes as c
import numpy as np
import Wlz as w
import math as m
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui

libc = c.CDLL("libc.so.6")

libc.fopen.restype = c.POINTER(w.FILE)

class WlzError(Exception): #{
  pass
#}

class WlzView(QtGui.QMainWindow): #{

  args = None
  current_path = ''
  version = '0.0.1'
  errnum = c.c_int(w.WLZ_ERR_NONE)
  # numpy image representation
  img = None
  img_itm = None
  img_view_box = None
  # primary Woolz object
  obj = None
  # objects and properties derived from the primary object
  obj2d = None
  obj_gtype = c.c_int(w.WLZ_GREY_ERROR)
  obj2d_org = [0, 0]
  obj2d_sz = [0, 0]
  # 3D sectioning parameters
  pitch = 0.0
  yaw = 0.0
  roll = 0.0
  dist = 0.0
  view = None
  # histogram/profile plot
  plt = None
  plt_itm = None
  prog = 'WlzView'
  # image ROI
  roi = None
  roi_type = 'N'
  # tracking of image values at cursor position
  track_value = False
  track_proxy = None
  # file formats understood
  file_formats = ['wlz']

  def __init__(self, prog, args): #{
    super(WlzView, self).__init__()
    self.prog = prog
    self.args = args
    self.initUI()
  #}

  def initUI(self): #{
    logging.debug('initUI()')
    ## Menubar
    menubar = self.menuBar()
    fMenu = menubar.addMenu('&File')
    aMenu = menubar.addMenu('&Actions')
    mMenu = menubar.addMenu('&Measurement')
    vMenu = menubar.addMenu('&View')
    menubar.addSeparator()
    hMenu = menubar.addMenu('&Help')
    # File menu
    fOpen = QtGui.QAction(QtGui.QIcon('open.wlz'), 'Open', self)
    fOpen.setShortcut('Ctrl+O')
    fOpen.setStatusTip('Open new Woolz file.')
    fOpen.triggered.connect(self.openFile)
    fMenu.addAction(fOpen)
    # 
    fExit  = QtGui.QAction(QtGui.QIcon('exit.png'), 'Exit', self)
    fExit.setShortcut('Ctrl+Q')
    fExit.setStatusTip('Exit Application')
    fExit.triggered.connect(self.quit)
    fMenu.addAction(fExit)
    # View menu
    vViewAll = QtGui.QAction('View All', self)
    vViewAll.setShortcut('Ctrl+A')
    vViewAll.setStatusTip('View All')
    vViewAll.triggered.connect(self.viewAll)
    vMenu.addAction(vViewAll)
    # Measurement menu
    curVal = mMenu.addAction('Value at cursor')
    curVal.setCheckable(True)
    curVal.triggered.connect(self.setTrackCurVal)
    roiMenu = mMenu.addMenu('&ROI')
    roiNone = roiMenu.addAction('None')
    roiLine = roiMenu.addAction('Line')
    roiRect = roiMenu.addAction('Rectangle')
    roiGrp = QtGui.QActionGroup(self)
    roiGrp.addAction(roiNone)
    roiGrp.addAction(roiLine)
    roiGrp.addAction(roiRect)
    roiNone.setCheckable(True)
    roiLine.setCheckable(True)
    roiRect.setCheckable(True)
    roiNone.setActionGroup(roiGrp)
    roiLine.setActionGroup(roiGrp)
    roiRect.setActionGroup(roiGrp)
    roiNone.triggered.connect(self.setROITypeNone)
    roiLine.triggered.connect(self.setROITypeLine)
    roiRect.triggered.connect(self.setROITypeRect)
    # Help menu
    self.about = AboutDialog(self.version)
    hAbout = QtGui.QAction('About', self)
    hAbout.setShortcut('Ctrl+b')
    hAbout.setStatusTip('About WlzView.')
    hAbout.triggered.connect(self.showAbout)
    hMenu.addAction(hAbout)
    #
    self.statusBar()
    self.setWindowTitle(self.prog)
    #
    w0 = QtGui.QWidget(self)
    self.setCentralWidget(w0)
    #
    grd0 = QtGui.QGridLayout(w0)
    grd0.setSpacing(8)
    gl0 = pg.GraphicsLayoutWidget()
    pw = pg.PlotWidget()
    ctl = QtGui.QFrame()
    grd0.addWidget(gl0, 0, 0, 8, 4)
    grd0.addWidget(pw, 0, 6, 2, 2)
    grd0.addWidget(ctl, 4, 6, 2, 2)
    #
    grd1 = QtGui.QGridLayout(ctl)
    grd1.setSpacing(8)
    gl1 = pg.GraphicsLayoutWidget()
    dst_lab = QtGui.QLabel('Distance')
    dst_sld = QtGui.QSlider(QtCore.Qt.Horizontal)
    dst_sld.setRange(0.0, 100.0) # TODO
    dst_sld.setValue(0.0) # TODO
    dst_val = QtGui.QLineEdit('0.0')
    pit_lab = QtGui.QLabel('Pitch')
    yaw_lab = QtGui.QLabel('Yaw')
    rol_lab = QtGui.QLabel('Roll')
    grd1.addWidget(dst_lab, 0, 0, 1, 1)
    dst_lab.setAlignment(QtCore.Qt.AlignLeft)
    grd1.addWidget(dst_sld, 0, 1, 1, 6)
    grd1.addWidget(dst_val, 0, 7, 1, 1)
    dst_val.setAlignment(QtCore.Qt.AlignRight)
    grd1.addWidget(pit_lab, 1, 0, 1, 1)
    grd1.addWidget(yaw_lab, 2, 0, 1, 1)
    grd1.addWidget(rol_lab, 3, 0, 1, 1)
    #
    self.img_view_box = gl0.addViewBox(0, 0)
    self.img_itm = pg.ImageItem()
    self.img_view_box.setAspectLocked()
    self.img_view_box.addItem(self.img_itm)
    self.img_view_box.invertY()
    #
    self.plt_itm = pw.getPlotItem()
    #
    if bool(self.args.infile): #{
      logging.debug('file given on cmdline')
      self.addObjFromFile(self.args.infile)
      logging.debug('added file from cmdline')
    #}
    #
    self.show()
  #}

  def openFile(self): #{
    logging.debug('openFile()')
    c = ''
    s = 'Woolz files ('
    for f in self.file_formats: #{
      s = s + c + '*.' + f
      c = ' '
    #}
    s = s + ');; All files (*)'
    path = QtGui.QFileDialog.getOpenFileName(self, 'Open Woolz file',
               self.current_path, s)
    if(path): #{
      p = str(path)
      self.addObjFromFile(p)
    #}
  #}

  def addObjFromFile(self, f): #{
    logging.debug('addObjFromFile()')
    obj = self.readWlzObj(f)
    if bool(obj): #{
      self.addObj(obj)
    #}
  #}
  
  def addObj(self, o): #{
    logging.debug('addObj()')
    t = int(o.contents.type)
    self.errnum = c.c_int(w.WLZ_ERR_OBJECT_TYPE)
    if (t == int(w.WLZ_2D_DOMAINOBJ)) or (t == int(w.WLZ_3D_DOMAINOBJ)): #{
      if t == int(w.WLZ_2D_DOMAINOBJ): #{
        self.errnum = c.c_int(w.WLZ_ERR_NONE)
      else: #}{
        logging.debug('making view struct')
        v = w.WlzMake3DViewStruct(c.c_int(w.WLZ_3D_VIEW_STRUCT), \
                                  c.byref(self.errnum))
        if not bool(self.errnum): #{
          f = w.WlzDVertex3()
          v.vtX = c.c_double(0.0)
          v.vtY = c.c_double(0.0)
          v.vtZ = c.c_double(0.0)
          v.contents.theta = self.yaw   * m.pi / 180.0
          v.contents.phi   = self.pitch * m.pi / 180.0
          v.contents.zeta  = self.roll  * m.pi / 180.0
          v.contents.dist  = self.dist
          v.contents.fixed = f
          v.contents.view_mode = c.c_int(w.WLZ_UP_IS_UP_MODE)
          v.contents.scale = c.c_double(1.0)
          v.contents.voxelRescaleFlg = c.c_int(0)
          w.WlzFree3DViewStruct(self.view)
          self.view = w.WlzAssign3DViewStruct(v, None)
        #}
      #}
    #}
    if not bool(self.errnum): #{
      w.WlzFreeObj(self.obj)
      self.obj = w.WlzAssignObject(o, None)
      self.errnum = self.setObj2D()
    #}
    if not bool(self.errnum): #{
      ary = self.wlz2DToNP()
      logging.debug('setting image')
      self.img = ary.astype(np.float64).T
      self.img_itm.setImage(self.img)
      self.setROIType('N')
    #}
    if bool(self.errnum): #{
      self.warnWlzError('Failed to add object.')
    #}
  #}

  def setObj2D(self): #{
    logging.debug('setObj2D')
    errnum = c.c_int(w.WLZ_ERR_NONE)
    t = self.obj.contents.type
    if t == w.WLZ_2D_DOMAINOBJ: #{
      self.obj2d = w.WlzAssignObject(self.obj, None)
    elif t == w.WLZ_3D_DOMAINOBJ: #}{
      logging.debug('cutting section from 3D object')
      w.WlzInit3DViewStruct(self.view, self.obj)
      o2d = w.WlzGetSubSectionFromObject(self.obj, None,
                self.view, c.c_int(w.WLZ_INTERPOLATION_NEAREST),
                None, c.byref(errnum))
      print(o2d)
      print(errnum)
      if not bool(errnum): #{
        w.WlzFreeObj(self.obj2d)
        self.obj2d = w.WlzAssignObject(o2d, None)
      #}
    #}
    return(errnum)
  #}

  def wlz2DToNP(self): #{
    ary = None
    logging.debug('wlz2DToNP()')
    try: #{
      box = w.WlzBoundingBox3I(self.obj2d, c.byref(self.errnum))
      if(bool(self.errnum)): #{
        raise WlzError()
      #}
      sz = w.WlzIVertex2()
      sz.vtX = box.xMax - box.xMin + 1
      sz.vtY = box.yMax - box.yMin + 1
      org = w.WlzIVertex2()
      org.vtX = box.xMin
      org.vtY = box.yMin
      gtype = w.WlzGreyTypeFromObj(self.obj2d, c.byref(self.errnum))
      if(bool(self.errnum)): #{
        raise WlzError()
      #}
      if gtype == w.WLZ_GREY_INT: #{
        vtype = c.c_int
      elif gtype == w.WLZ_GREY_SHORT: #}{
        vtype = c.c_short
      elif gtype == w.WLZ_GREY_UBYTE: #}{
        vtype = c.c_ubyte
      elif gtype == w.WLZ_GREY_FLOAT: #}{
        vtype = c.c_float
      elif gtype == w.WLZ_GREY_DOUBLE: #}{
        vtype = c.c_double
      else: #}{
        self.errnum = c.c_int(w.WLZ_ERR_GREY_TYPE)
        raise WlzError()
      #}
      UPP = c.POINTER(c.POINTER(vtype))
      UPV = c.POINTER(c.c_void_p)
      aryc = c.cast(0,UPV)
      self.errnum = w.WlzToArray2D(c.byref(aryc), self.obj2d, sz, org, 0, \
                                   c.c_int(gtype))
      if(bool(self.errnum)): #{
        raise WlzError()
      #}
      aryc = c.cast(aryc, UPP)
      ary = np.ctypeslib.as_array(aryc.contents, (sz.vtY, sz.vtX))
      self.obj_gtype = gtype
      self.obj2d_sz = [sz.vtX, sz.vtY]
      self.obj2d_org = [org.vtX, org.vtY]
      w.AlcFree(aryc)
    except WlzError: #}{
      self.warnWlzError('Failed to extract numeric data from object.')
    #}
    return(ary)
  #}

  def setTrackCurVal(self, q): #{
    logging.debug('setTrackCurVal(' + str(q) + ')')
    if bool(q): #{
      self.track_value = True
      self.track_proxy = pg.SignalProxy(self.img_itm.scene().sigMouseMoved, \
          rateLimit=60, slot=self.trackCurVal)
      logging.debug('setTrackCurVal added event handler')
    else: #}{
      self.track_value = False
      self.track_proxy = None
      logging.debug('setTrackCurVal removed event handler')
    #}
  #}

  def trackCurVal(self, evt): #{
    logging.debug('trackCurVal(' + str(evt) + ')')
    pos = evt[0]
    if self.img_itm.sceneBoundingRect().contains(pos): #{
      q = self.img_view_box.mapSceneToView(pos)
      p = [int(q.x()), int(q.y())]
      v = self.img[p[0], p[1]]
      g = w.WlzStringFromGreyType(self.obj_gtype, None)
      msg = str(g) + ' ' + str(p) + ' ' + str(v)
      self.statusBar().showMessage(msg)
    #}
  #}

  def setROITypeNone(self): #{
    logging.debug('setROITypeNone()')
    self.setROIType('N')
  #}

  def setROITypeLine(self): #{
    logging.debug('setROITypeLine()')
    self.setROIType('L')
  #}

  def setROITypeRect(self): #{
    logging.debug('setROITypeRect()')
    self.setROIType('R')
  #}

  def updateROI(self): #{
    logging.debug('updateROI()')
    data = self.roi.getArrayRegion(self.img, self.img_itm)
    if self.roi_type == 'L': #{
      self.plt_itm.clear()
      self.plt_itm.plot(data)
    elif self.roi_type == 'R': #}{
      self.plt_itm.clear()
      y, x = np.histogram(data, 256)
      self.plt_itm.plot(x, y, stepMode=True)
    #}
  #}

  def setROIType(self, t): #{
    logging.debug('setROIType()')
    if bool(self.roi): #{
      logging.debug('removing current roi')
      self.img_view_box.removeItem(self.roi)
    #}
    logging.debug('roi is now of type ' + t)
    if (t == 'N'): #{
      self.roi = None
      self.plt_itm.clear()
      if bool(self.plt): #{
        self.plt_itm.clear()
        self.plt = self.plt_itm.plot()
      #}
      self.plt_itm.setTitle('Image Histogram')
      y, x = np.histogram(self.img, 256)
      self.plt_itm.plot(x, y, stepMode=True)
    elif (t == 'L'): #}{
      self.roi = pg.LineSegmentROI([[self.obj2d_org[0], self.obj2d_org[1]],
                                    [20, 20]], pen=(0,9))
      self.roi.sigRegionChanged.connect(self.updateROI)
      self.img_view_box.addItem(self.roi)
      if bool(self.plt): #{
        self.plt_itm.clear()
        self.plt = self.plt_itm.plot()
      #}
      self.plt_itm.setTitle('Line Profile')
    elif (t == 'R'): #}{
      self.roi = pg.RectROI([self.obj2d_org[0], self.obj2d_org[1]],
                            [20, 20], pen=(0,9))
      self.roi.sigRegionChanged.connect(self.updateROI)
      self.img_view_box.addItem(self.roi)
      if bool(self.plt): #{
        self.plt_itm.clear()
        self.plt = self.plt_itm.plot()
      #}
      self.plt_itm.setTitle('ROI Histogram')
    #}
    self.roi_type = t
  #}

  def readWlzObj(self, fname): #{
    logging.debug('readWlzObj()')
    obj = None
    logging.debug('attempting to read object from ' + fname)
    self.errnum = c.c_int(w.WLZ_ERR_FILE_OPEN)
    fp = c.cast(libc.fopen(fname, 'rb'), c.POINTER(w.FILE))
    if bool(fp): #{
      obj = w.WlzAssignObject( \
            w.WlzReadObj(fp, c.byref(self.errnum)), None)
      libc.fclose(fp)
      if bool(self.errnum): #{
        self.warnWlzError('Failed to read object from ' + fname);
      #}
    #}
    return(obj)
  #}

  def quit(self): #{
    logging.debug('quit()')
    sys.exit(0)
  #}

  def warnWlzError(self, msg): #{
    logging.debug('warnWlzError')
    logging.debug('error message is: ' + msg)
    err = str(w.WlzStringFromErrorNum(self.errnum, None))
    logging.debug('error code is: ' + err)
    QtGui.QMessageBox.critical(self,
        'Warning',
        msg + ' (' + err + ')',
        QtGui.QMessageBox.Close)
    self.errnum = c.c_int(w.WLZ_ERR_NONE);
  #}

  def viewAll(self): #{
    logging.debug('viewAll')
    if bool(self.obj2d): #{
      self.img_view_box.autoRange()
    #}
  #}

  def showAbout(self): #{
    logging.debug('showAbout')
    self.about.show()
  #}
#}

class AboutDialog(QtGui.QDialog): #{

  def __init__(self, version, parent = None): #{
    super(AboutDialog, self).__init__(parent)
    self.initUI(version)
  #}

  def initUI(self, version): #{
    msg = '<b>WlzView</b>' + \
          '<p>Version ' + version + '.</p>' + \
          '<p>A simple pyqtgraph based Woolz image viewer.</p>' + \
          '<p>ma-tech@igmm.ed.ac.uk</p>'
    lo = QtGui.QVBoxLayout(self)
    txt = QtGui.QLabel(msg)
    lo.addWidget(txt)
    but = QtGui.QDialogButtonBox(
              QtGui.QDialogButtonBox.Cancel,
              QtCore.Qt.Horizontal, self)
    but.rejected.connect(self.cancel)
    lo.addWidget(but)
  #}

  def cancel(self): #{
    self.hide()
  #}
#}


def parseArgs(): #{
  parser = argparse.ArgumentParser(description =
      'A simple interactive Woolz object viewer written using PyWoolz.')
  parser.add_argument(
      '-l', '--log',
      dest='logLevel',
      choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
      help='Set the logging level.')
  parser.add_argument('infile',
      help='Input Woolz file.')
  args = parser.parse_args()
  if args.logLevel: #{
    logging.basicConfig(level=logging.getLevelName(args.logLevel))
  #}
  return(args)
#}

def main(): #{
  prog = sys.argv[0]
  app = QtGui.QApplication(sys.argv)
  args = parseArgs()
  wv = WlzView(prog, args)
  sys.exit(app.exec_())
  
#}

if __name__ == '__main__': #{
  main()
#}
  


