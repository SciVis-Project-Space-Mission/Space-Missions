from typing import Any, List, Dict
import spiceypy as spice
import vtk
import argparse
import os
import numpy as np
from scipy import interpolate
import math
import time
import json
import faulthandler
import sys
from datetime import datetime
import matplotlib.pyplot as plt


def make_reference_frame(origin: np.ndarray, X: np.ndarray, Y: np.ndarray, Z: np.ndarray, size: float) -> vtk.vtkActor:
    vectors = vtk.vtkFloatArray()
    vectors.SetNumberOfComponents(3)
    vectors.SetNumberOfTuples(3)
    X *= size
    Y *= size
    Z *= size
    vectors.SetTuple3(0, X[0], X[1], X[2])
    vectors.SetTuple3(1, Y[0], Y[1], Y[2])
    vectors.SetTuple3(0, Z[0], Z[1], Z[2])

    scalars = vtk.vtkFloatArray()
    scalars.SetNumberOfComponents(1)
    scalars.SetNumberOfTuples(3)
    scalars.SetTuple1(0, 0)
    scalars.SetTuple1(1, 1)
    scalars.SetTuple1(2, 2)

    pos = vtk.vtkFloatArray()
    pos.SetNumberOfComponents(3)
    pos.SetNumberOfTuples(3)
    pos.SetTuple3(0, origin[0], origin[1], origin[2])
    pos.SetTuple3(1, origin[0], origin[1], origin[2])
    pos.SetTuple3(2, origin[0], origin[1], origin[2])
    pts = vtk.vtkPoints()
    pts.SetData(pos)

    data = vtk.vtkPolyData()
    data.SetPoints(pts)
    data.GetPointData().SetVectors(vectors)
    data.GetPointData().SetScalars(scalars)

    arrow = vtk.vtkArrowSource()
    glyph = vtk.vtkGlyph3D()
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetInputData(data)
    glyph.Update()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(glyph.GetOutput())

    ctf = vtk.vtkColorTransferFunction()
    ctf.AddRGBPoint(0, 0, 0, 1)
    ctf.AddRGBPoint(1, 0, 1, 0)
    ctf.AddRGBPoint(2, 1, 0, 0)
    mapper.SetLookupTable(ctf)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    return actor

'''
Observed issues:
- camera focal point correction not working if time is paused
- Unable to see Clipper when it is the focus object


Note to self:
    - time_scale is in minutes per second (m/s)
    - timeline position is in days (d)
    - refresh rate (in Hz) is fixed (Hz=1/s)
    - time step (in seconds) is time_scale * Units.minute / refresh rate 

Actors:
    - Planet 
    - Planet position
    - Planet orbit
    - Moon
    - Moon position
    - Clipper position
    - Clipper forward orbit
    - Clipper backward orbit
    - Clipper orbit
    - Ecliptic plane
    - Ecliptic frame

SimulationState.params:
    - time_step
    - clock 
    - cam_position_update_policy
    - came_focus_update policy
    - saved_relative_position
    - video_recording
    - video file basename 
    - snaphot file basename 
    

Needed cameras:
    - Cape Canaveral pointing at the sky
    - one anchored at Clipper
    - one pointing at Clipper 
    - tethered to one of the main bodies: Sun, 8 planets, Moon, 4 moons of Jupiter: 14 bodies
    - Ecliptic

    Each body can carry its own camera with its own default setting that the user can change through
    interaction.
'''

from PyQt5.QtWidgets import QApplication, QWidget, QMainWindow, QCheckBox, QGroupBox, QHBoxLayout, QVBoxLayout, QComboBox, QRadioButton, QGridLayout, QCalendarWidget, QMenuBar, QMenu, QStatusBar, QAction, QGridLayout, QLabel, QPushButton, QSpinBox, QDial, QSpacerItem, QSizePolicy, QDoubleSpinBox, QDateTimeEdit
import PyQt5.QtCore as QtCore
from PyQt5.QtCore import Qt, QTimer, QDate, QTime, QDateTime
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

verbose = False       

def toint(alist):
    return [ int(x) for x in alist ]

'''
Ui_MainWindow:
A Qt GUI that includes a vtkWidget with vtkRenderWindow and some control widgets for the 
visualization of the Clipper mission
'''
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        self.targets = { 0: 'None', 1: 'Sun', 2: 'Mercury', 3: 'Venus', 4: 'Earth', 5: 'Mars', 6: 'Jupiter', 7: 'Saturn', 8: 'Uranus', 9: 'Neptune', 10: 'Moon', 11: 'Europa', 12: 'Io', 13: 'Ganymede', 14: 'Callisto', 15: 'Clipper', 16: 'Sky'}
        self.anchors = { 0: 'None', 1: 'Cape Canaveral', 2: 'Clipper' }
        self.timesteps = [ '1 minute', '15 minutes', '30 minutes', '1 hour', '6 hours', '1 day', '1 week', '1 month' ]
        self.target_to_id = {}
        self.anchor_to_id = {}
        for i in self.targets.keys():
            self.target_to_id[self.targets[i]] = i
        for i in self.anchors.keys():
            self.anchor_to_id[self.anchors[i]] = i

        self.events = ['No Event', 'launch', 'Mars assist', 'Earth assist',
                       'Jupiter capture', 'Europa flyby']

        MainWindow.setObjectName("MainWindow")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.mainGridLayout = QGridLayout(self.centralwidget)
        self.mainGridLayout.setObjectName("mainGridLayout")

        self.vtkWidget = QVTKRenderWindowInteractor(self.centralwidget)
        self.vtkWidget.setObjectName("vtkWindow")
        self.mainGridLayout.addWidget(self.vtkWidget, 0, 0, 1, 5)

        self.spatialControlBox = QGroupBox(self.centralwidget)
        #
        self.spatialControlBox.setObjectName("spatialControlBox")
        self.spaceControlLayout = QVBoxLayout(self.spatialControlBox)
        self.spaceControlLayout.setObjectName("spaceControlLayout")
        self.cameraTargetGroupBox = QGroupBox(self.spatialControlBox)
        self.cameraTargetGroupBox.setObjectName("cameraTargetGroupBox")
        self.cameraTargetLayout = QHBoxLayout(self.cameraTargetGroupBox)
        self.cameraTargetLayout.setObjectName("cameraTargetLayout")
        self.cameraTargetComboBox = QComboBox(self.cameraTargetGroupBox)
        self.cameraTargetComboBox.setObjectName("cameraTargetComboBox")
        for i in range(len(self.targets)):
            self.cameraTargetComboBox.addItem('')
        self.cameraTargetLayout.addWidget(self.cameraTargetComboBox)
        self.resetCameraTarget = QCheckBox(self.cameraTargetGroupBox)
        self.cameraTargetLayout.addWidget(self.resetCameraTarget)
        self.spaceControlLayout.addWidget(self.cameraTargetGroupBox)
        #
        self.cameraPositionGroupBox = QGroupBox(self.spatialControlBox)
        self.cameraPositionGroupBox.setObjectName("cameraPositionGroupBox")
        self.cameraPositionLayout = QHBoxLayout(self.cameraPositionGroupBox)
        self.cameraPositionLayout.setObjectName("cameraPositionLayout")
        self.cameraPositionComboBox = QComboBox(self.cameraPositionGroupBox)
        self.cameraPositionComboBox.setObjectName("cameraPositionComboBox")
        for i in range(len(self.anchors)):
            self.cameraPositionComboBox.addItem('')
        self.cameraPositionLayout.addWidget(self.cameraPositionComboBox)
        self.resetCameraPosition = QCheckBox(self.cameraPositionGroupBox)
        self.cameraPositionLayout.addWidget(self.resetCameraPosition)
        #
        self.spaceControlLayout.addWidget(self.cameraPositionGroupBox)
        #
        self.scaleGroupBox = QGroupBox(self.spatialControlBox)
        self.scaleGroupBox.setObjectName("scaleGroupBox")
        self.scaleGroupLayout = QHBoxLayout(self.scaleGroupBox)
        self.scaleGroupLayout.setObjectName("scaleGroupLayout")
        self.scaleLabel = QLabel(self.scaleGroupBox)
        self.scaleLabel.setObjectName("scaleLabel")
        self.scaleGroupLayout.addWidget(self.scaleLabel)
        self.scaleDial = QDial(self.scaleGroupBox)
        self.scaleDial.setObjectName("scaleDial")
        self.scaleGroupLayout.addWidget(self.scaleDial)
        self.scaleSpinBox = QSpinBox(self.scaleGroupBox)
        self.scaleSpinBox.setObjectName("scaleSpinBox")
        self.scaleGroupLayout.addWidget(self.scaleSpinBox)
        self.spaceControlLayout.addWidget(self.scaleGroupBox)
        self.mainGridLayout.addWidget(self.spatialControlBox, 1, 0, 1, 1)

        self.timeControlBox = QGroupBox(self.centralwidget)
        self.timeControlBox.setObjectName("timeControlBox")
        self.timeControlLayout = QVBoxLayout(self.timeControlBox)
        self.timeControlLayout.setObjectName("timeControlLayout")
        self.timeStepLabel = QLabel(self.timeControlBox)
        self.timeStepLabel.setObjectName("timeStepLabel")
        self.timeControlLayout.addWidget(self.timeStepLabel)
        self.pausePushButton = QPushButton(self.timeControlBox)
        self.pausePushButton.setObjectName("pausePushButton")
        self.timeControlLayout.addWidget(self.pausePushButton)
        self.timestepRadioButton = []
        for i, t in enumerate(self.timesteps):
            self.timestepRadioButton.append(QRadioButton(self.timeControlBox))
            self.timestepRadioButton[i].setObjectName(f"timestepRadioButton[{i}]")
            self.timeControlLayout.addWidget(self.timestepRadioButton[i])
        self.mainGridLayout.addWidget(self.timeControlBox, 1, 1, 1, 1)

        self.timePointBox = QGroupBox(self.centralwidget)
        self.timePointBox.setObjectName("timePointBox")
        self.timePointLayout = QVBoxLayout(self.timePointBox)
        self.timePointLayout.setObjectName("timePointLayout")
        # self.calendar = QCalendarWidget(self.centralwidget)
        self.calendar = QDateTimeEdit(self.timePointBox)
        self.calendar.setObjectName("calendar")
        self.timePointLayout.addWidget(self.calendar)
        self.eventComboBox = QComboBox(self.timePointBox)
        self.eventComboBox.setObjectName("eventComboBox")
        for i in range(len(self.events)):
            self.eventComboBox.addItem('')
        self.timePointLayout.addWidget(self.eventComboBox)
        self.mainGridLayout.addWidget(self.timePointBox, 1, 2, 1, 1)

        self.cameraControl = QGroupBox(self.centralwidget)
        self.cameraControl.setObjectName('cameraControl')
        self.cameraControlLayout = QVBoxLayout(self.cameraControl)
        self.cameraControlLayout.setObjectName('cameraControlLayout')
        #
        self.viewAngleControl = QGroupBox(self.cameraControl)
        self.viewAngleLayout = QHBoxLayout(self.viewAngleControl)
        self.viewAngleLabel = QLabel(self.viewAngleControl)
        self.viewAngleDial = QDial(self.viewAngleControl)
        self.viewAngleSpinBox = QSpinBox(self.viewAngleControl)
        self.viewAngleLayout.addWidget(self.viewAngleLabel)
        self.viewAngleLayout.addWidget(self.viewAngleDial)
        self.viewAngleLayout.addWidget(self.viewAngleSpinBox)
        self.cameraControlLayout.addWidget(self.viewAngleControl)
        #
        self.clippingControl = QGroupBox(self.cameraControl)
        self.clippingLayout = QHBoxLayout(self.clippingControl)
        self.clippingNearSpinBox = QDoubleSpinBox(self.clippingControl)
        self.clippingNearLabel = QLabel(self.clippingControl)
        self.clippingFarSpinBox = QDoubleSpinBox(self.clippingControl)
        self.clippingFarLabel = QLabel(self.clippingControl)
        self.clippingLayout.addWidget(self.clippingNearLabel)
        self.clippingLayout.addWidget(self.clippingNearSpinBox)
        self.clippingLayout.addWidget(self.clippingFarLabel)
        self.clippingLayout.addWidget(self.clippingFarSpinBox)
        self.cameraControlLayout.addWidget(self.clippingControl)
        #
        self.mainGridLayout.addWidget(self.cameraControl, 1, 3, 1, 1)

        self.recordingControl = QGroupBox(self.centralwidget)
        self.recordingControl.setObjectName("recordingControl")
        self.recordingLayout = QVBoxLayout(self.recordingControl)
        self.recordingLayout.setObjectName("recordingLayout")
        self.snapshotPushButton = QPushButton(self.recordingControl)
        self.snapshotPushButton.setObjectName("snapshotPushButton")
        self.recordingLayout.addWidget(self.snapshotPushButton)
        self.startRecordingPushButton = QPushButton(self.recordingControl)
        self.startRecordingPushButton.setObjectName("startRecordingPushButton")
        self.recordingLayout.addWidget(self.startRecordingPushButton)
        self.pauseRecordingPushButton = QPushButton(self.recordingControl)
        self.pauseRecordingPushButton.setObjectName("pauseRecordingPushButton")
        self.recordingLayout.addWidget(self.pauseRecordingPushButton)
        self.stopRecordingPushButton = QPushButton(self.recordingControl)
        self.stopRecordingPushButton.setObjectName("stopRecordingPushButton")
        self.recordingLayout.addWidget(self.stopRecordingPushButton)
        self.quitPushButton = QPushButton(self.recordingControl)
        self.quitPushButton.setObjectName("quitPushButton")
        self.recordingLayout.addWidget(self.quitPushButton)
        self.mainGridLayout.addWidget(self.recordingControl, 1, 4, 1, 1)
        # self.spacerLayout = QVBoxLayout(self.centralwidget)
        self.spacer = QSpacerItem(10, 800, QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.mainGridLayout.addItem(self.spacer, 0, 5, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        # self.menubar.setGeometry(QtCore.QRect(0, 0, 1023, 24))
        self.menubar.setObjectName("menubar")
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuView = QMenu(self.menubar)
        self.menuView.setObjectName("menuView")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionSelect_Output_Path = QAction(MainWindow)
        self.actionSelect_Output_Path.setObjectName("actionSelect_Output_Path")
        self.actionQuit = QAction(MainWindow)
        self.actionQuit.setObjectName("actionQuit")
        self.actionFull_Screen = QAction(MainWindow)
        self.actionFull_Screen.setObjectName("actionFull_Screen")
        self.menuFile.addAction(self.actionSelect_Output_Path)
        self.menuFile.addAction(self.actionQuit)
        self.menuView.addAction(self.actionFull_Screen)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuView.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def index_to_target(self, index: int) -> str:
        return self.targets[index]
    
    def target_to_index(self, target: str) -> int:
        return self.target_to_id[target]

    def index_to_anchor(self, index: int) -> str:
        return self.anchors[index]

    def anchor(self, anchor: str) -> int:
        return self.anchor_to_id[anchor]

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "LayoutTest"))
        self.cameraTargetGroupBox.setTitle(
            _translate("MainWindow", "Camera Focus Point"))
        for i, t in self.targets.items():
            self.cameraTargetComboBox.setItemText(i, _translate('MainWindow', t))
        self.resetCameraTarget.setText(_translate('MainWindow', 'Reset'))
        self.cameraPositionGroupBox.setTitle(
            _translate("MainWindow", "Camera Position"))
        for i, a in self.anchors.items():
            self.cameraPositionComboBox.setItemText(
                i, _translate("MainWindow", a))
        self.resetCameraPosition.setText(_translate('MainWindow', 'Reset'))
        self.scaleLabel.setText(_translate("MainWindow", "Scale"))

        for i, e in enumerate(self.events):
            self.eventComboBox.setItemText(i, _translate('MainWindow', e))

        self.timeStepLabel.setText(_translate("MainWindow", "Time Step"))
        self.timeStepLabel.setProperty(
            "", _translate("MainWindow", "Time step"))
        self.pausePushButton.setText(_translate("MainWindow", "Pause"))
        self.pausePushButton.setProperty("", _translate("MainWindow", "Pause"))
        for i, t in enumerate(self.timesteps):
            self.timestepRadioButton[i].setText(_translate("MainWindow", t))
            self.timestepRadioButton[i].setProperty(
                "", _translate("MainWindow", t))
            
        self.viewAngleLabel.setText(_translate('MainWindow', 'View angle'))
        self.viewAngleLabel.setProperty("", _translate("MainWindow", 'View angle'))
        self.clippingNearLabel.setText(_translate('MainWindow', 'Near clip (in km)'))
        self.clippingNearLabel.setProperty("", _translate("MainWindow", 'Near clip (in km)'))
        self.clippingFarLabel.setText(_translate('MainWindow', 'Far clip (in AU)'))
        self.clippingFarLabel.setProperty("", _translate("MainWindow", 'Far clip (in AU)'))
            
        self.snapshotPushButton.setText(_translate("MainWindow", "Snapshot"))
        self.snapshotPushButton.setProperty(
            "", _translate("MainWindow", "Snapshot"))
        self.startRecordingPushButton.setText(
            _translate("MainWindow", "Start Recording"))
        self.startRecordingPushButton.setProperty(
            "", _translate("MainWindow", "Start Recording"))
        self.pauseRecordingPushButton.setText(
            _translate("MainWindow", "Pause Recording"))
        self.pauseRecordingPushButton.setProperty(
            "", _translate("MainWindow", "Pause Recording"))
        self.stopRecordingPushButton.setText(
            _translate("MainWindow", "Stop Recording"))
        self.stopRecordingPushButton.setProperty(
            "", _translate("MainWindow", "Stop Recording"))
        self.quitPushButton.setText(_translate("MainWindow", "Quit"))
        self.quitPushButton.setProperty("", _translate("MainWindow", "Quit"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuView.setTitle(_translate("MainWindow", "View"))
        self.actionSelect_Output_Path.setText(
            _translate("MainWindow", "Select Output Path"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))
        self.actionFull_Screen.setText(_translate("MainWindow", "Full Screen"))

'''
Units:
Handy conversions
'''
class Units:

    class date:
        def __init__(self, seconds=0, minutes=0, hours=0, days=0):
            self.t = seconds + minutes*60 + hours*3600 + days*86400
        
        def time(self):
            return self.t 
        
    def interval(d1: date, d2: date) -> int:
        return d2.time()-d1.time()

    def duration(seconds=0, minutes=0, hours=0, days=0, weeks=0, months=0, years=0, output='seconds') -> float:
        delta = seconds + minutes * 60 + hours * 3600 + \
            days * 86400 + weeks * 604800 + months * 2592000 + years * 31557600
        if output == 'seconds':
            return delta 
        elif output == 'minutes':
            return delta/60.
        elif output == 'hours':
            return delta/3600.
        elif output == 'days':
            return delta/86400.
            
    # handy unit conversions
    def au_to_km(au: float) -> float:
        return au * 149597871.0

    def km_to_au(km: float) -> float:
        return km / 149597871.0

    def deg_to_rad(deg: float) -> float:
        return deg/180.*math.pi
    
    def rad_to_deg(rad: float) -> float:
        return rad/math.pi*180.

    def to_unit_rgb(rgb: np.ndarray) -> List:
        r = float(rgb[0]/255.)
        g = float(rgb[1]/255.)
        b = float(rgb[2]/255.)
        return [r, g, b]

    def monthStr2Int(month: str) -> int:
        months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
                'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
        return months.index(month) + 1

    def int2monthStr(monthid: int) -> str:
        months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
                'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
        return months[monthid-1]

    def time2QDate(et: float) -> QDate:
        timestr = spice.timout(int(et), 'MON/DD/YYYY')
        return QDate(int(timestr[7:]), Units.monthStr2Int(timestr[0:3]), int(timestr[4:6]))
    
    def time2QDateTime(et: float) -> QDateTime:
        timestr = spice.timout(int(et), 'MON/DD/YYYY HR:MN:SC')
        date = QDate(int(timestr[7:11]), Units.monthStr2Int(timestr[0:3]), int(timestr[4:6]))
        time = QTime(int(timestr[12:14]), int(timestr[15:17]), int(timestr[18:20]))
        return QDateTime(date, time)

    def QDate2time(date: QDate) -> float:
        timestr = f'{date.year()} {Units.int2monthStr(date.month())}, {date.day()} 00:00:00'
        print(f'convert {timestr} to epoch')
        return spice.utc2et(timestr)
    
    def QDateTime2time(datetime: QDateTime) -> float:
        timestr = f'{datetime.date().year()} {Units.int2monthStr(datetime.date().month())}, {datetime.date().day()} {datetime.time().hour()}:{datetime.time().minute()}:{datetime.time().second()}'
        print(f'convert {timestr} to epoch')
        return spice.utc2et(timestr)

# set frame and corrections
default_frame = 'ECLIPJ2000'
abcorr = 'NONE'

# last modification, interrupted, attempted to add the ability to control
# the time range of the visualization, its time step, and the automatic
# recording of each individual frame, as well as the specification of the
# camera setting, including in terms of the bodies and Clipper's position

def myprint(verbose: bool, text: str):
    if verbose:
        print(text)

'''
Body:
A class to represent a massive body (sun, planet, or moon), including shape, size, position, 
orientation, coordinates transformation, and graphical representation
'''
class Body(object):
    def __init__(self, name: str, texture_name: str, meridian: float, radii: np.ndarray, color: np.ndarray):
        self.fullname = name
        self.name = name.split()[0]
        self.img_filename = texture_name #texture_files[self.name]
        sphere = vtk.vtkTexturedSphereSource()
        sphere.SetRadius(1)
        sphere.SetThetaResolution(40)
        sphere.SetPhiResolution(40)

        shape = vtk.vtkTransform()
        shape.RotateZ(meridian) #-meridian_rotation[self.name])
        self.radii = radii # get_radii(self.name)
        # print(f'self.radii={self.radii}')
        shape.Scale(self.radii[0], self.radii[1], self.radii[2])
        transform = vtk.vtkTransformPolyDataFilter()
        transform.SetInputConnection(sphere.GetOutputPort())
        transform.SetTransform(shape)
        transform.Update()
        self.no_point = True
        self.ellipsoid = transform.GetOutput()
        # no (scale-independent) point representation for moons
        if name != 'Moon' and name != 'Callisto' and name != 'Io' \
            and name != 'Europa' and name != 'Ganymede':
            self.no_point = False
        if not self.no_point:
            self.point_actor = VTKUtils.make_point_actor([0,0,0], [1,1,1], 4)
            self.point_actor.SetObjectName(f'{self.name} point')
            self.point_actor.GetProperty().RenderPointsAsSpheresOn()
            self.point_actor.GetProperty().SetColor(color[0], color[1], color[2])
        else:
            self.point_actor = vtk.vtkActor()
            self.point_actor.SetObjectName(f'{self.name} point')
        if name == 'Sun':
            self.point_actor.GetProperty().SetAmbient(1000)
            self.point_actor.SetObjectName(f'Sun point')

        self.local_frame = 'IAU_' + self.name.upper()
        self.states = None
        self.times = None
        self.actor = None
        self.orbit = None
        self.orbit_frame = None

    def get_max_radius(self) -> float:
        return np.amax(self.radii)

    def set_states(self, states: List[np.ndarray], ets: np.ndarray, frame: str):
        self.states = states
        self.times = ets
        self.orbit_frame = frame
        self.orbit = interpolate.make_interp_spline(self.times, self.states)
        self.deriv = self.orbit.derivative()

    def get_actor(self) -> vtk.vtkActor:
        # Body in its own reference frame
        if self.actor is None:
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(self.ellipsoid)
            mapper.ScalarVisibilityOff()
            self.actor = vtk.vtkActor()
            self.actor.SetObjectName(self.name)
            self.actor.SetMapper(mapper)
            self.actor.GetProperty().SetDiffuse(0.1)
            self.actor.GetProperty().SetSpecular(0.05)
            reader = vtk.vtkJPEGReader()
            reader.SetFileName(self.img_filename)
            self.texture = vtk.vtkTexture()
            self.texture.SetInputConnection(reader.GetOutputPort())
            self.actor.SetTexture(self.texture)
        return self.actor

    def scale_actor(self, scaling: float):
        self.actor.SetScale(scaling, scaling, scaling)

    def get_position(self, et: float, frame=default_frame) -> np.ndarray:
        if self.orbit is not None:
            if self.orbit_frame == frame:
                return self.orbit(et)
            else:
                myprint(verbose=verbose, text='doing coordinate transformation')
                return DataLoader.coordinate_transform(self.orbit(et), self.orbit_frame, frame, et)
        else:
            myprint(verbose=verbose, text='Invalid request: get_position: no states available')
        return np.array([0.,0.,0.])
    
    def get_tangent(self, et: float, frame=default_frame) -> np.ndarray:
        if self.deriv is not None:
            if self.orbit_frame == frame:
                return self.deriv(et)
            else:
                myprint(verbose=verbose, text='doing coordinate transformation')
                return DataLoader.coordinate_transform(self.deriv(et), self.orbit_frame, frame, et)
        else:
            myprint(verbose=verbose,
                    text='Invalid request: get_position: no states available')
        return np.array([0., 0., 0.])

    def get_rotation(self, et: float, frame=default_frame):
        mat = DataLoader.body_orientation(self.name, et, frame)
        return mat

    def get_matrix(self, et: float, frame=default_frame) -> vtk.vtkMatrix4x4:
        mat = self.get_rotation(et, frame)
        pos = self.get_position(et, frame)
        matrix = vtk.vtkMatrix4x4()
        matrix.Identity()
        for row in range(3):
            matrix.SetElement(row, 0, mat[row, 0])
            matrix.SetElement(row, 1, mat[row, 1])
            matrix.SetElement(row, 2, mat[row, 2])
            matrix.SetElement(row, 3, pos[row])
        return matrix;

    def get_local_position(self, pos: np.ndarray, et: float) -> np.ndarray:
        mat = self.get_rotation(et)
        center =self.get_position(et)
        dp = pos-center
        return np.matmul(dp, mat)

    def local_to_global_position(self, pos: np.ndarray, et: float) -> np.ndarray:
        mat = self.get_rotation(et)
        center = self.get_position(et)
        return center + np.matmul(mat, pos)

    def get_equatorial_plane(self, et: float, frame=default_frame) -> np.ndarray:
        mat = self.get_rotation(et, frame)
        p = np.array(self.get_position(et, frame))
        e1 = mat[:,0]
        e2 = mat[:,1]
        return [p, e1, e2]

    def move_actor(self, et: float, frame=default_frame) -> vtk.vtkActor:
        actor = self.get_actor()
        actor.SetUserMatrix(self.get_matrix(et, frame))
        return actor

    def move_point_actor(self, et: float, frame=default_frame) -> vtk.vtkActor:
        actor = self.point_actor
        actor.SetUserMatrix(self.get_matrix(et, frame))
        return actor
    
    def is_intersected(self, et: float, obs: int, pos: np.ndarray, target: int):
        target_pos, et0 = spice.spkezp(targ=target, et=et, ref=default_frame, abcorr='XLT', obs=obs)
        ray = target_pos - pos
        intcp, et1, dir, found = spice.sincpt(method='ELLIPSOID', target=self.fullname, et=et, fixref='IAU_'+self.name, abcorr='XLT', obsrv=obs, dref=default_frame, dvec=ray)
        if found:
            print(f'{self.name} was intersected by ray from {obs} to {target} at time {et}')
            print(f'Time to {target}: {et0}, time to intercept: {et1}')
        return found

def duplicate_camera(camera: vtk.vtkCamera) -> vtk.vtkCamera:
    cam = vtk.vtkCamera()
    cam.DeepCopy(camera)
    return cam

def copy_camera(from_camera: vtk.vtkCamera, to_camera: vtk.vtkCamera) -> vtk.vtkCamera:
    to_camera.DeepCopy(from_camera)
    return to_camera

'''
Schedule:
A class to encapsulate the date of the various events of the mission and the time intervals 
in betweeen
''' 
class Schedule:
    def __init__(self, init: float, arrival: float, end: float):
        self.initial_time = init
        self.launch_time = init + Units.duration(days=0.5)
        self.launch_duration = Units.duration(hours=2)
        print(f'launch time is {spice.timout(int(self.launch_time), "AP: MN: SC AMPM Month DD, YYYY")}')
        self.arrival_time = arrival
        self.final_time = end
        # set up time array
        # total duration of the covered portion of Clipper in days
        self.full_timeline = np.arange(init, end, Units.duration(days=1))
        # 1.5 day from launch in minutes intervals
        self.prelaunch_timeline = np.arange(init, self.launch_time, Units.duration(minutes=1))
        self.launch_phase = np.arange(self.launch_time, self.launch_time+self.launch_duration)
        # interplanetary travel in days intervals (includes time to launch)
        self.pretour_timeline = np.arange(init, arrival, Units.duration(days=1))
        # moon tour in hours
        self.tour_timeline = np.arange(arrival, end, Units.duration(hours=1))
        self.transfer_timeline = np.arange(self.launch_time+self.launch_duration, arrival, Units.duration(days=1))
        self.mission_timeline= np.concatenate((self.prelaunch_timeline, self.launch_phase, self.transfer_timeline, self.tour_timeline))
        self.active_mission_timeline = np.concatenate((self.launch_phase, self.transfer_timeline, self.tour_timeline))

'''
Clipper:
A class to represent the Clipper spacecraft, including position, orbit, and graphical representation
'''
class Clipper:
    def make_orbit(self, clock):
        forward_points = []
        backward_points = []
        if clock < self.launch:
            return
        # first 2 hours after launch in minutes, then hours
        if clock < self.launch + Units.duration(hours=2):
            timeline = np.arange(self.launch, clock,
                                    Units.duration(minutes=1))
        else:
            timeline_launch = np.arange(self.launch, Units.duration(
                hours=2), Units.duration(minutes=1))
            timeline_after = np.arange(
                self.launch + Units.duration(hours=2), clock, Units.duration(hours=1))
            timeline = np.concatenate((timeline_launch, timeline_after))
        for t in timeline:
            p_clipper = self.intp(t)
            forward_points.append(p_clipper)
            # p_launch_site = self.earth.local_to_global_position(
            #     self.launch_site, t)
            # backward_points.append(p_launch_site)
        # print('there are {} points in foward points and {} points in backward points'.format(len(forward_points), len(backward_points)))
        self.forward_orbit_actor = VTKUtils.make_curve_actor(
            forward_points, [0.5, 0.5, 0], actor=self.forward_orbit_actor)
        # self.backward_orbit_actor = VTKUtils.make_curve_actor(
        #     backward_points, [0.5, 0, 0.5], actor=self.backward_orbit_actor)

    def __init__(self, data, clock=None):
        self.earth = data.bodies['Earth']
        self.schedule = data.schedule
        self.has_launched = False
        # print(f'mission timeline = {self.schedule.mission_timeline}')
        # print(f'clipper steps = {data.full_paths["clipper"]}')
        self.intp = interpolate.make_interp_spline(x=self.schedule.mission_timeline, y=data.full_paths['clipper'])
        self.deriv = self.intp.derivative()
        self.deriv2 = self.deriv.derivative()
        self.launch = self.schedule.launch_time
        print('clipper launch set to {} / '.format(self.launch, 
            spice.timout(int(self.launch), 'AP:MN:SC AMPM Month DD, YYYY')))
        p_clipper = self.intp(self.schedule.initial_time)
        self.launch_site = self.earth.get_local_position(p_clipper, self.schedule.initial_time)
        p_clipper[2] += 10 # raise launch site above Earth surface
        # compute local position of launching site in Earth rotating coordinate system

        self.position_actor = VTKUtils.make_point_actor(p_clipper, [1,1,0], 10)
        self.position_actor.SetObjectName('Clipper position')
        self.position_actor.GetProperty().RenderPointsAsSpheresOn()
        self.forward_orbit_actor = vtk.vtkActor()
        self.forward_orbit_actor.SetObjectName('Clipper forward orbit')
        # self.backward_orbit_actor = vtk.vtkActor()
        # self.backward_orbit_actor.SetObjectName('Clipper backward orbit')
        if clock is not None and clock >= self.launch:
            print('Clipper launched at initialization')
            self.has_launched = True
            self.make_orbit(clock)
        self.forward_orbit_actor.GetProperty().RenderLinesAsTubesOn()
        self.forward_orbit_actor.GetProperty().BackfaceCullingOff()
        self.forward_orbit_actor.GetProperty().FrontfaceCullingOff()
        # self.backward_orbit_actor.GetProperty().RenderLinesAsTubesOn()
        # self.backward_orbit_actor.GetProperty().FrontfaceCullingOff()
        # self.backward_orbit_actor.GetProperty().FrontfaceCullingOff()

    def update(self, et: float):
        p_clipper = self.intp(et)
        self.position_actor.SetPosition(p_clipper)
        if not self.has_launched and et >= self.launch:
            print(f'Clipper just launched at time {et}!!!')
            self.has_launched = True
            self.make_orbit(et)
        elif self.has_launched:
            print(f'adding one point to curve at t={et}')
            myprint(verbose=verbose, text='adding one point to curve')
            self.forward_orbit_actor = VTKUtils.add_point_to_curve(self.forward_orbit_actor, p_clipper)
            # p_earth = self.earth.local_to_global_position(self.launch_site, et)
            # self.backward_orbit_actor = VTKUtils.add_point_to_curve(self.backward_orbit_actor, p_earth)
            self.forward_orbit_actor.Modified()
            # self.backward_orbit_actor.Modified()

    def get_position(self, et: float) -> np.ndarray:
        return self.intp(et)

    def get_velocity(self, et: float) -> np.ndarray:
        return self.deriv(et)

    def get_acceleration(self, et: float) -> np.ndarray:
        return self.deriv2(et)
    
    def get_tangent(self, et: float) -> np.ndarray:
        return self.deriv(et)
    
    def get_actors(self):
        return { 'Clipper position': self.position_actor, #'Clipper backward orbit': self.backward_orbit_actor, 
                 'Clipper forward orbit': self.forward_orbit_actor }

'''
Frame:
A simple static frame with conversion to and from reference coordinates
'''
class Frame:
    def __init__(self, origin: np.ndarray, X: np.ndarray, Y: np.ndarray):
        self.origin = np.array(origin)
        self.X = X/np.linalg.norm(X)
        self.Y = Y - np.dot(Y, self.X)*self.X
        self.Y = self.Y/np.linalg.norm(self.Y)
        self.Z = np.cross(self.X, self.Y)
        self.actor = vtk.vtkActor()
        self.size = 1

    def __repr__(self) -> str:
        return 'Frame'
    
    def __str__(self) -> str:
        return f'origin:\n\t{self.origin}\nX:\n\t{self.X},\nY:\n\t{self.Y}\nZ:\n\t{self.Z}'

    def global_to_local_vector(self, world_vec: np.ndarray, et=0, frame=None) -> np.ndarray:
        world_vec = np.array(world_vec)
        x = np.dot(world_vec, self.X)
        y = np.dot(world_vec, self.Y)
        z = np.dot(world_vec, self.Z)
        return np.array([x, y, z])
    
    def get_frame(self, et=0):
        return self

    def local_to_global_vector(self, local_vec: np.ndarray, et=0, frame=None) -> np.ndarray:
        return local_vec[0]*self.X + local_vec[1]*self.Y + local_vec[2]*self.Z

    def global_to_local_coordinates(self, world_pos: np.ndarray, et=0, frame=None) -> np.ndarray:
        dp = np.array(world_pos) - self.origin
        return self.global_to_local_vector(dp, et)

    def local_to_global_coordinates(self, local_pos: np.ndarray, et=0, frame=None) -> np.ndarray:
        p = self.origin + self.local_to_global_vector(local_pos, et)
        return p
    
    def set_actor_size(self, size: float):
        self.size = size

    def get_actor(self) -> vtk.vtkActor:
        self.actor = make_reference_frame(origin=self.origin, X=self.X, Y=self.Y, Z=self.Z, size=self.size)
        return self.actor

'''
DynamicFrame:
A time-dependent, quasi-inertial, non-rotating, barycenter based extension of Frame
'''
class DynamicFrame(Frame):
    def __init__(self, object, center: Body):
        self.object = object
        self.center = center
        self.size = 1
        self.actor = vtk.vtkActor()

    def __repr__(self) -> str:
        return 'DynamicFrame'

    def __str__(self) -> str:
        return f'object:\n{self.object}\ncenter:\n{self.center}'

    def get_frame(self, et: float) -> Frame:
        origin = self.object.get_position(et)
        X = self.center.get_position(et) - origin
        xsave = X
        X[2] = 0
        # cheap tangent approximation
        # Y = origin - self.object.get_position(et+Units.duration(hours=1))
        Z = np.array([0, 0, 1])
        Y = np.cross(Z, X) 
        print(f'in get frame: X={X}, Y={Y}')
        return Frame(origin, X, Y)

    def global_to_local_vector(self, world_vec: np.ndarray, et: float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
        return f.global_to_local_vector(world_vec)

    def local_to_global_vector(self, local_vec: np.ndarray, et: float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
        return f.local_to_global_vector(local_vec)

    def global_to_local_coordinates(self, position: np.ndarray, et: float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
            print(f'local frame is:\n{f}')
            print(f'position is {position}')
        return f.global_to_local_coordinates(position)

    def local_to_global_coordinates(self, position: np.ndarray, et: float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
            print(f'local frame is:\n{f}')
            print(f'position is {position}')
        return f.local_to_global_coordinates(position)
    
    def get_actor(self, et:float, frame=None) -> vtk.vtkActor:
        if frame is None:
            frame = self.get_frame(et)
        frame.set_actor_size(self.size)
        self.actor = frame.get_actor()
        return self.actor

'''
ClipperBodyDynamicFrame:
An extension of DynamicFrame to handle the special case of a frame that keeps the relative position
of Clipper to a reference body constant. 
'''
class ClipperBodyDynamicFrame(DynamicFrame):
    def __init__(self, clipper: Clipper, target: Body, reference: Body):
        self.clipper = clipper
        self.object = self.clipper
        self.target = target
        self.center = self.target
        self.reference = reference
        self.size = 1
        self.actor = vtk.vtkActor()

    def get_frame(self, et: float) -> Dict[str, np.ndarray]:
        origin = self.clipper.get_position(et)
        targetp = self.target.get_position(et)
        Y = targetp - origin
        Ynorm = np.linalg.norm(Y)
        Z = np.array([0, 0, 1], dtype=np.float64)
        Y /= Ynorm
        Z -= np.dot(Z, Y)*Y
        Z /= np.linalg.norm(Z)
        X = np.cross(Y, Z)
        return { 'origin': origin, 'X': X, 'Y': Y, 'Z': Z, 'Ynorm': Ynorm }

    def global_to_local_vector(self, world_vec: np.ndarray, et: float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
        vx = np.dot(world_vec, f['X'])
        vy = np.dot(world_vec, f['Y'])/f['Ynorm']
        vz = np.dot(world_vec, f['Z'])
        return np.array([vx, vy, vz])
    
    def local_to_global_vector(self, local_vec: np.ndarray, et: float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
        return local_vec[0]*f['X'] + local_vec[1]*f['Ynorm']*f['Y'] + local_vec[2]*f['Z']

    def global_to_local_coordinates(self, position: np.ndarray, et: float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
        dp = position - f['origin']
        return self.global_to_local_vector(dp, et, f)
    
    def local_to_global_coordinates(self, position: np.ndarray, et: float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
        return f['origin'] + self.local_to_global_vector(position, et, f)
    
class ClipperFrame(DynamicFrame):
    def __init__(self, clipper: Clipper):
        self.clipper = clipper
        self.size = 1
        self.actor = vtk.vtkActor()
        self.size = 1

    def get_frame(self, et: float) -> Dict[str, np.ndarray]:
        origin = self.clipper.get_position(et)
        forward = self.clipper.get_position(et + 100)
        X = forward - origin
        X /= np.linalg.norm(X)
        Z = np.array([0, 0, 1], dtype=np.float64)
        Y = np.cross(Z, X)
        return { 'origin': origin, 'X': X, 'Y': Y, 'Z': Z }

'''
Conversions between local and global reference frames for camera settings
'''
def local2global(frame: Frame, localcam: vtk.vtkCamera, et: float) -> vtk.vtkCamera:
    globalcam = vtk.vtkCamera()
    f = frame.get_frame(et)
    print(f'local2global: frame:\n{f}')
    globalcam.SetPosition(frame.local_to_global_coordinates(localcam.GetPosition(), et, f))
    globalcam.SetFocalPoint(frame.local_to_global_coordinates(localcam.GetFocalPoint(), et, f))
    globalcam.SetViewUp(frame.local_to_global_vector(localcam.GetViewUp(), et, f))
    globalcam.SetViewAngle(localcam.GetViewAngle())
    globalcam.SetClippingRange(localcam.GetClippingRange())
    return globalcam

def global2local(frame: Frame, globalcam: vtk.vtkCamera, et: float) -> vtk.vtkCamera:
    localcam = vtk.vtkCamera()
    f = frame.get_frame(et)
    print(f'global2local: frame:\n{f}')
    localcam.SetPosition(frame.global_to_local_coordinates(globalcam.GetPosition(), et, f))
    localcam.SetFocalPoint(frame.global_to_local_coordinates(globalcam.GetFocalPoint(), et, f))
    localcam.SetViewUp(frame.global_to_local_vector(globalcam.GetViewUp(), et, f))
    localcam.SetViewAngle(globalcam.GetViewAngle())
    localcam.SetClippingRange(globalcam.GetClippingRange())
    return localcam


'''
CapeCanaveralFrame:
A dynamic frame anchored at Cape Canaveral and rotating with the Earth
'''
class CapeCanaveralFrame(DynamicFrame):
    def __init__(self, earth: Body, cc: np.ndarray):
        self.earth = earth
        self.cc = cc
        self.size=1
        self.actor = vtk.vtkActor()

    def get_frame(self, et: float) -> Frame:
        origin = self.earth.local_to_global_position(self.cc, et)
        center = self.earth.get_position(et)
        up = origin-center 
        up /= np.linalg.norm(up)
        ploc = self.cc + np.array([1.0, 0, 0])
        X = self.earth.local_to_global_position(ploc, et)
        X -= np.dot(X, up)*up
        X /= np.linalg.norm(X)
        Y = np.cross(up, X)
        return Frame(origin=origin, X=X, Y=Y)

    def local_to_global_vector(self, vector: np.ndarray, et: float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
        return f.local_to_global_vector(vector)
    
    def local_to_global_coordinates(self, position: np.ndarray, et: float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
        return f.local_to_global_coordinates(position)
    
    def global_to_local_vector(self, vector: np.ndarray, et:float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
        return f.global_to_local_vector(vector)
    
    def global_to_local_coordinates(self, position: np.ndarray, et:float, f=None) -> np.ndarray:
        if f is None:
            f = self.get_frame(et)
        return f.global_to_local_coordinates(position)

'''
DataLoader:
A class to import and store all the NAIF information about the mission as well as media assets
used for texture mapping. This class also stores the precomputed orbits of all massive bodies
and of Clipper along with their associated frames
'''
class DataLoader:
    # Static NAIF/SPICE helper functions
    def get_radii(body_name: str) -> List[float]:
        radii = spice.bodvrd( body_name.split()[0], 'RADII', 3 )
        myprint(verbose=verbose, text='radii={}'.format(radii))
        return radii[1]
 
    def coordinate_transform(coord: List[float], from_frame:str, to_frame: str, et: float) -> List[float]:
        mat = spice.pxform(from_frame, to_frame, et)
        return spice.mxv(mat, coord)
  
    def body_orientation(body_name: str, et: float, to_frame: str, from_frame=None) -> np.ndarray:
        if from_frame is None:
            from_frame = 'IAU_' + body_name
        return spice.pxform(from_frame, to_frame, et)

    # Import all the necessary NAIF data and sample trajectories
    def __init__(self, media_path='./media', naif_path='./naif', scale=1):
        # load the data
        myprint(verbose=verbose, text='loading the data from NAIF files...')
        # spice.furnsh(os.path.join(naif_path, '19F23_VEEGA_L230511_A290930_LP01_V2_scpse.bsp')) # clipper
        # spice.furnsh(os.path.join(naif_path, '21F31_MEGA_L241010_A300411_LP01_V4_postLaunch_scpse.bsp'))
        spice.furnsh(os.path.join(naif_path, '21F31_MEGA_L241010_A300411_LP01_V5_pad_scpse.bsp'))
        spice.furnsh(os.path.join(naif_path, 'naif0012.tls')) # bodies' dynamics
        spice.furnsh(os.path.join(naif_path, 'pck00010.tpc')) # bodies' constant values and orientation
        myprint(verbose=verbose, text='...done')

        # set clipper SPICE ID
        self.spice_id    = -159
        self.spice_idstr = '-159'

        self.planet_names_spice = \
            [ 'Mercury', 'Venus', 'Earth', 'Mars Barycenter', \
              'Jupiter', 'Saturn barycenter', 'Uranus barycenter', \
              'Neptune barycenter' ]
        self.planet_names = [ name.split()[0] for name in self.planet_names_spice ]
        self.jupiter_moon_names = [ 'Io', 'Europa', 'Ganymede', 'Callisto' ]
        self.earth_moon_name = [ 'Moon' ]
        self.moon_names = self.earth_moon_name + self.jupiter_moon_names
        self.all_body_names = ['Sun'] + self.planet_names + self.moon_names
        self.all_body_names_spice = ['Sun'] + self.planet_names_spice + self.moon_names

        self.saturn_rings = { # inner and outer radii or Saturn's rings
                'D': (66900, 74510), 'C': (74658, 92000),
                'B': (92000, 117580), 'Cassini': (117580, 122170),
                'A': (122170, 136775), 'Roche': (136775, 139380),
                'F': (140183, 140680)
                }

        self.texture_files = {
                'Sun':     os.path.join(media_path, '2k_sun_bright.jpg'),
                'Mercury': os.path.join(media_path, '2k_mercury.jpg'),
                'Venus':   os.path.join(media_path, '2k_venus.jpg'),
                # 'Earth':   os.path.join(media_path, '2k_earth_daymap.jpg'),
                'Earth':   os.path.join(media_path, 'earth/world.topo.200412.3x5400x2700.jpg'),
                'Mars':    os.path.join(media_path, '2k_mars.jpg'),
                'Jupiter': os.path.join(media_path, '2k_jupiter.jpg'),
                'Saturn':  os.path.join(media_path, '2k_saturn.jpg'),
                'Uranus':  os.path.join(media_path, '2k_uranus.jpg'),
                'Neptune': os.path.join(media_path, '2k_neptune.jpg'),
                'Stars':   os.path.join(media_path, 'starry_background.jpg'),
                'Rings':   os.path.join(media_path, 'ring_colors2d.jpg'),
                'Io':      os.path.join(media_path, '2k_io.jpg'),
                'Europa':  os.path.join(media_path, '2k_europa.jpg'),
                'Ganymede':os.path.join(media_path, '2k_ganymede.jpg'),
                'Callisto':os.path.join(media_path, '2k_callisto.jpg'),
                'Moon':    os.path.join(media_path, '2k_moon.jpg')
                }

        self.meridian_rotation = {
                'Sun': 0,
                'Mercury': 201,
                'Venus': 180,
                'Earth': 180,
                'Mars': 180,
                'Jupiter': 232.5,
                'Saturn': 0,
                'Uranus': 0,
                'Neptune': 0,
                'Io': 0,
                'Europa': 0,
                'Ganymede': 0,
                'Callisto': 0,
                'Moon': 180
                }

        self.orbital_period = {
                'Sun': 0,
                'Mercury': Units.duration(days=88),
                'Venus': Units.duration(days=225),
                'Earth': Units.duration(days=365.256),
                'Mars': Units.duration(days=686.971),
                'Jupiter': Units.duration(days=4332.59),
                'Saturn': Units.duration(days=10759.22),
                'Uranus': Units.duration(days=30688.5),
                'Neptune': Units.duration(days=60182),
                'Io': Units.duration(days=4332.59),
                'Europa': Units.duration(days=4332.59),
                'Ganymede': Units.duration(days=4332.59),
                'Callisto': Units.duration(days=4332.59),
                'Moon': Units.duration(days=365.256),
                 }

        # colors to be used for single-point representation of distant bodies
        self.planet_colors = {
                'Sun':      Units.to_unit_rgb([255,167,0]),
                'Mercury':  Units.to_unit_rgb([112,109,113]),
                'Venus':    Units.to_unit_rgb([240,218,166]),
                'Earth':    Units.to_unit_rgb([44,49,115]),
                'Mars':     Units.to_unit_rgb([244,108,76]),
                'Jupiter':  Units.to_unit_rgb([240,210,145]),
                'Saturn':   Units.to_unit_rgb([200,170,100]),
                'Uranus':   Units.to_unit_rgb([193,238,238]),
                'Neptune':  Units.to_unit_rgb([81,106,161]),
                'Io':       Units.to_unit_rgb([252,240,113]),
                'Europa':   Units.to_unit_rgb([148,115,68]),
                'Ganymede': Units.to_unit_rgb([81,70,59]),
                'Callisto': Units.to_unit_rgb([86,77,53]),
                'Moon':     Units.to_unit_rgb([128, 128, 128])
                }

        self.radii = {}
        for name in self.all_body_names:
            self.radii[name] = DataLoader.get_radii(name)

        # set frame and corrections
        self.default_frame = 'ECLIPJ2000'
        self.abcorr = 'NONE'
 
        myprint(verbose=verbose, text='Setting temporal domain for Clipper mission...')
        # coverage dates for Clipper
        etb = spice.cell_double(10000)
        # entire duration of the mission
        # spice.spkcov(os.path.join(naif_path, '19F23_VEEGA_L230511_A290930_LP01_V2_scpse.bsp'), self.spice_id, etb)
        # spice.spkcov(os.path.join(naif_path, '21F31_MEGA_L241010_A300411_LP01_V4_postLaunch_scpse.bsp'), self.spice_id, etb)
        spice.spkcov(os.path.join(naif_path, '21F31_MEGA_L241010_A300411_LP01_V5_pad_scpse.bsp'), self.spice_id, etb)
        # arrival time
        arrival_time = spice.str2et('2030 APR 11 00:00:00 TDB')
        init_time = etb[0]
        final_time = etb[1]
        # selected_interval = [init_time, final_time]
        print(f'init time={init_time}, final_time={final_time}')
        print(f'init_time={spice.etcal(init_time)}, final_time={spice.etcal(final_time)}')
        self.schedule = Schedule(init=init_time, arrival=arrival_time, end=final_time)

        # query states during entire mission for clipper and relevant planets
        # states are given in (km, km/s)
        '''
        Bodies orbits during mission time frame
        '''
        self.full_paths = {}
        self.interpl_paths = {}
        self.tour_paths = {}
        myprint(verbose=verbose, text='Acquiring bodies\' paths across time domain...')
        for id, name in enumerate(self.all_body_names):
            myprint(verbose=verbose, text=' * {}'.format(name))
            spice_name = self.all_body_names_spice[id]
            self.full_paths[name], not_used = spice.spkezr(spice_name, self.schedule.full_timeline, default_frame, abcorr, 'Sun')
            self.interpl_paths[name], not_used = spice.spkezr(spice_name, self.schedule.transfer_timeline, default_frame, abcorr, 'Sun')
            self.full_paths[name] = scale*np.array(self.full_paths[name])[:,:3]
            self.interpl_paths[name] = scale*np.array(self.interpl_paths[name])[:,:3]

        myprint(verbose=verbose, text='Acquiring Jupiter\'s moons\' path across time domain')
        for name in self.jupiter_moon_names:
            myprint(verbose=verbose, text=' * {}'.format(name))
            self.tour_paths[name], not_used = spice.spkezr(name, self.schedule.mission_timeline, default_frame, abcorr, 'Jupiter')
            self.tour_paths[name] = scale*np.array(self.tour_paths[name])[:,:3]

        # clipper trajectory
        myprint(verbose=verbose, text=' * Clipper')
        self.clipper_full, not_used = spice.spkezr(self.spice_idstr, self.schedule.mission_timeline, default_frame, abcorr, 'Sun')
        self.full_paths['clipper'] = scale*np.array(self.clipper_full)[:,:3]
        myprint(verbose=verbose, text='...done')

        objreader = vtk.vtkOBJReader()
        objreader.SetFileName(os.path.join(media_path, 'clipper_spacecraft.obj'))
        objreader.Update()
        clipper_model = objreader.GetOutput()
        clipper_mapper = vtk.vtkPolyDataMapper()
        clipper_mapper.SetInputData(clipper_model)
        clipper_mapper.ScalarVisibilityOff()
        self.clipper_model_actor = vtk.vtkActor()
        self.clipper_model_actor.SetMapper(clipper_mapper)
        self.clipper_model_actor.GetProperty().BackfaceCullingOff()
        self.clipper_model_actor.GetProperty().FrontfaceCullingOff()
        self.clipper_model_actor.GetProperty().SetAmbient(0.4)
        # self.clipper_model_actor.GetProperty().SetInterpolationToPBR()
        # self.clipper_model_actor.GetProperty().SetMetallic(1.0)
        # self.clipper_model_actor.GetProperty().SetColor(0.5, 0.5, 0.5)
        # spacecraft dimensions are in dm (?!), convert to km
        # self.clipper_model_actor.SetScale(0.0001, 0.0001, 0.0001)
                              
        # def __init__(self, name, texture_name, meridian, radii, color)
        myprint(verbose=verbose, text='Initializing body objects...')
        self.bodies = {}
        self.tags = {}
        for name in self.all_body_names:
            self.bodies[name] = Body(name, self.texture_files[name], self.meridian_rotation[name], self.radii[name], self.planet_colors[name])
            self.bodies[name].set_states(self.full_paths[name], self.schedule.full_timeline, default_frame)

        self.clipper = Clipper(self)
        self.dynbodyframe = {}
        for name in self.planet_names:
            self.dynbodyframe[('None', name)] = DynamicFrame(self.bodies[name], self.bodies['Sun'])
            self.dynbodyframe[('None', name)].set_actor_size(1.5*self.bodies[name].get_max_radius())
            self.dynbodyframe[('Clipper', name)] = ClipperBodyDynamicFrame(clipper=self.clipper, target=self.bodies[name], reference=self.bodies['Sun'])
            self.dynbodyframe[('Clipper', name)].set_actor_size(400)

        self.dynbodyframe[('None', 'Moon')] = DynamicFrame(self.bodies['Moon'], self.bodies['Earth'])
        self.dynbodyframe[('None', 'Moon')].set_actor_size(1.5*self.bodies['Moon'].get_max_radius())
        self.dynbodyframe[('Clipper', 'Moon')] = ClipperBodyDynamicFrame(clipper=self.clipper, target=self.bodies['Moon'], reference=self.bodies['Earth'])
        self.dynbodyframe[('Clipper', 'Moon')].set_actor_size(1.5*self.bodies['Moon'].get_max_radius())
        for name in self.jupiter_moon_names:
            self.dynbodyframe[('None', name)] = DynamicFrame(self.bodies[name], self.bodies['Jupiter'])
            self.dynbodyframe[('None', name)].set_actor_size(1.5*self.bodies[name].get_max_radius()) 
            self.dynbodyframe[('Clipper', name)] = ClipperBodyDynamicFrame(clipper=self.clipper, target=self.bodies[name], reference=self.bodies['Jupiter'])
            self.dynbodyframe[('Clipper', name)].set_actor_size(1.5*self.bodies[name].get_max_radius())
        self.dynbodyframe[('None', 'Sun')] = Frame(origin=[0,0,0], X=[0,0,1], Y=[0,1,0])
        self.dynbodyframe[('None', 'Sun')].set_actor_size(1.5*self.bodies['Sun'].get_max_radius())
        self.dynbodyframe[('Clipper', 'None')] = DynamicFrame(self.clipper, self.bodies['Sun'])
        self.dynbodyframe[('Clipper', 'None')].set_actor_size(400)
        self.dynbodyframe[('None', 'Clipper')] = DynamicFrame(self.clipper, self.bodies['Sun'])
        self.dynbodyframe[('None', 'Clipper')].set_actor_size(400)
        self.dynbodyframe[('Clipper', 'Clipper')] = DynamicFrame(self.clipper, self.bodies['Sun'])
        self.dynbodyframe[('Clipper', 'Clipper')].set_actor_size(400)
        self.dynbodyframe[('Cape Canaveral', 'Sky')] = CapeCanaveralFrame(earth=self.bodies['Earth'], cc=self.clipper.launch_site)
        self.dynbodyframe[('Cape Canaveral', 'Sky')].set_actor_size(1000)

    def distance(self, body1: str, body2: str, et: float, frame=default_frame) -> float:
        p1 = self.bodies[body1].get_position(et, frame)
        p2 = self.bodies[body2].get_position(et, frame)
        return np.linalg.norm(p1-p2)

    def reference_body(self, name: str) -> str:
        if name == 'Sun' or name in self.planet_names:
            return 'Sun'
        elif name in self.jupiter_moon_names:
            return 'Jupiter'
        elif name == 'Moon':
            return 'Earth'
        elif name == 'Clipper':
            return 'Earth'
        else:
            myprint(verbose=verbose, text='unrecognized body name: {}'.format(name))

def what_day(data, et):
    return (et-data.schedule.initial_time)/Units.duration(days=1)

'''
VTK helper functions
'''
class VTKUtils:
    def print_actors(actors):
        print('\nactors\n')
        for k in actors.keys():
            v = actors[k]
            if isinstance(v, dict):
                names = [sk for sk in v.keys()]
                print(f'{k}: {names}')
            else:
                print(k)

    def get_renderer(window: vtk.vtkRenderWindow) -> vtk.vtkRenderer:
        return window.GetRenderers().GetFirstRenderer()

    def make_point_actor(pos=[0,0,0], color=[1,0,0], size=5):
        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        verts = vtk.vtkCellArray()
        verts.InsertNextCell(1)
        verts.InsertCellPoint(0)
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetVerts(verts)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)
        mapper.ScalarVisibilityOff()
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.SetPosition(pos[0], pos[1], pos[2])
        actor.GetProperty().SetColor(color[0], color[1], color[2])
        actor.GetProperty().SetPointSize(size)
        return actor

    def make_curve_actor(pts: List[np.ndarray], color=[0.5, 0.5, 0], size=3, actor=None):
        points = vtk.vtkPoints();
        for p in pts:
            points.InsertNextPoint(p[0], p[1], p[2])
        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(len(pts))
        for i in range(len(pts)):
            polyline.GetPointIds().SetId(i, i)
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(polyline)
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetLines(lines)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)
        mapper.ScalarVisibilityOff()
        if actor is None:
            actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetLineWidth(size)
        actor.GetProperty().SetColor(color[0], color[1], color[2])
        return actor

    def add_point_to_curve(curve_actor: vtk.vtkActor, point: np.ndarray) -> vtk.vtkActor:
        # print('adding point to actor')
        polydata = curve_actor.GetMapper().GetInput()
        points = polydata.GetPoints()
        n = points.GetNumberOfPoints()
        # print('there are currently {} points in polydata'.format(n))
        points.InsertNextPoint(point[0], point[1], point[2])
        points.GetData().Modified()
        lines = polydata.GetLines()
        # print('adding lines')
        lines.InsertNextCell(2)
        lines.InsertCellPoint(n-1)
        lines.InsertCellPoint(n)
        lines.GetData().Modified()
        # print('resetting points')
        # polydata.SetPoints(points)
        # polydata.SetLines(lines)
        # print('about to invoke set input data')
        # curve_actor.GetMapper().SetInputData(polydata)
        # print('reset actor input about to activate update')
        curve_actor.GetMapper().Update()
        # print('Update done')
        return curve_actor

    def copy_transformation(actor_from: vtk.vtkActor, actor_to: vtk.vtkActor):
        actor_to.SetUserMatrix(actor_from.GetUserMatrix())

    '''
    Default viewpoints for various camera location / target combinations
    '''
    def vantage_point(camera_name: str, data: DataLoader, camera_pos=None, camera_view_up=None) -> vtk.vtkCamera:
        
        camera = vtk.vtkCamera()
        far = 10000000

        if camera_name == ('Clipper', 'Clipper') or camera_name == ('None', 'Clipper'):
            # tethered to spacecraft, about 30m (0.03km) wide
            camera.SetPosition(0, 10, 0)
            camera.SetFocalPoint(0, 0, 0)
            camera.SetViewAngle(90)
            camera.SetViewUp(0, 0, 1)
            camera.SetClippingRange(0.01, far)
        elif camera_name == ('Clipper', 'None'):
            # anchored at spacecraft, looking ahead
            camera.SetPosition(0, 10, 0)
            camera.SetFocalPoint(0, 0, 0)
            camera.SetViewUp(0, 0, 1)
            camera.SetViewAngle(45)
            camera.SetClippingRange(0.05, far)
        elif camera_name == ('Cape Canaveral', 'Sky'):
            camera.SetPosition(0, 0, 0)
            camera.SetFocalPoint(0, 0, 100)
            camera.SetViewAngle(45)
            camera.SetViewUp(0, 1, 0)
            camera.SetClippingRange(0.1, far)
        elif camera_name[0] == 'None':
            # tethered to a planet, following trajectory
            bodyname = camera_name[1]
            r = data.bodies[bodyname].get_max_radius()
            if camera_pos is not None:
                camera.SetPosition(*camera_pos)
            if camera_view_up is not None:
                camera.SetViewUp(*camera_view_up)
            else:
                camera.SetPosition(5*r, 0, 0)
                camera.SetViewUp(0, 0, 1)
            camera.SetViewAngle(45)
            camera.SetFocalPoint(0, 0, 0)
            camera.SetClippingRange(r, 50*r)
        elif camera_name[0] == 'Clipper':
            # anchored at spacecraft 
            camera.SetPosition(0, -0.05, 0)
            camera.SetFocalPoint(0, 0.5*Units.au_to_km(1), 0)
            camera.SetViewUp(0, 0, 1)
            camera.SetViewAngle(45)
            camera.SetClippingRange(-0.04, far)
        else:
            print(f'ERROR: unrecognized camera name: {camera_name}')
            raise ValueError(f'ERROR: unrecognized camera name: {camera_name}')
        
        return camera

    def clock_as_str(cl: float) -> str:
        return spice.timout(int(cl), 'AP:MN:SC AMPM Month DD, YYYY')

    def sun_mercury_distance(data: DataLoader) -> float:
        return data.distance('Mercury', 'Sun', data.schedule.initial_time)

    def sun_radius(data: DataLoader) -> float:
        return data.bodies['Sun'].get_max_radius()

    def get_camera(window: vtk.vtkRenderWindow) -> vtk.vtkRenderer:
        return VTKUtils.get_renderer(window).GetActiveCamera()

    def camera_description(camera: vtk.vtkCamera) -> str:
        pos = camera.GetPosition()
        pos_str = '({:06.2f}, {:06.2f}, {:06.2f})'.format(pos[0], pos[1], pos[2])
        foc = camera.GetFocalPoint()
        foc_str = '({:06.2f}, {:06.2f}, {:06.2f})'.format(foc[0], foc[1], foc[2])
        up = camera.GetViewUp()
        up_str = '({:06.2f}, {:06.2f}, {:06.2f})'.format(up[0], up[1], up[2])
        angle = camera.GetViewAngle()
        angle_str = '{:06.2f}'.format(angle)
        arange = camera.GetClippingRange()
        return 'p={}, f={}, u={}, a={}, range={}'.format(pos_str, foc_str, up_str, angle_str, arange)

    # Return: ring_actors 
    # actors corresponding to major ring groups
    def make_rings(data: DataLoader) -> List[vtk.vtkActor]:
        ring_reader = vtk.vtkJPEGReader()
        ring_reader.SetFileName(data.texture_files['Rings'])
        ring_texture = vtk.vtkTexture()
        ring_texture.SetInputConnection(ring_reader.GetOutputPort())

        rings_info = data.saturn_rings
        minR = rings_info['C'][0]
        maxR = rings_info['A'][1]
        ring_actors = []
        for ring_name in ['C', 'B', 'A']:
            inner, outer = rings_info[ring_name]
            tc = vtk.vtkFloatArray()
            tc.SetNumberOfComponents(2)
            source = vtk.vtkDiskSource()
            source.SetInnerRadius(inner)
            source.SetOuterRadius(outer)
            source.SetCircumferentialResolution(90)
            source.SetRadialResolution(2)
            source.Update()
            aring = source.GetOutput()
            # compute texture coordinates
            points = aring.GetPoints()
            for i in range(points.GetNumberOfPoints()):
                p = points.GetPoint(i)
                r = math.sqrt(p[0]*p[0] + p[1]*p[1])
                a = math.atan2(p[1], p[0])
                y = (r-minR)/(maxR-minR)
                x = 0.5*a/math.pi
                tc.InsertNextTuple((y,x))
            aring.GetPointData().SetTCoords(tc)
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(aring)
            actor = vtk.vtkActor()
            actor.SetObjectName(f'Saturn\'s ring {ring_name}')
            actor.SetMapper(mapper)
            actor.SetTexture(ring_texture)
            ring_actors.append(actor)
        return ring_actors

    # Return: actor_plane, actor_perimeter
    # Ecliptic disk and surrounding perimeter
    # TODO: command line option: AttributeError: 'MainWindow' object has no attribute 'graphics'
    def make_ecliptic_plane(data: DataLoader) -> List[vtk.vtkActor]:
        ecliptic_actors = []
        ecliptic = vtk.vtkRegularPolygonSource()
        ecliptic.SetCenter(0,0,0)
        ecliptic.SetNumberOfSides(500)
        ecliptic.SetRadius(data.distance('Neptune', 'Sun', data.schedule.initial_time))
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(ecliptic.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetObjectName('ecliptic')
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0,0.2,1)
        actor.GetProperty().SetOpacity(0.05)
        actor.GetProperty().SetAmbient(1)
        ecliptic_actors.append(actor)
        frame = vtk.vtkRegularPolygonSource()
        frame.SetCenter(0,0,0)
        frame.SetNumberOfSides(500)
        frame.SetRadius(data.distance('Neptune', 'Sun', data.schedule.initial_time))
        frame.GeneratePolygonOff()
        f_mapper = vtk.vtkPolyDataMapper()
        f_mapper.SetInputConnection(frame.GetOutputPort())
        f_actor = vtk.vtkActor()
        f_actor.SetObjectName('ecliptic frame')
        f_actor.SetMapper(f_mapper)
        f_actor.GetProperty().SetColor(0,0.2,1)
        f_actor.GetProperty().SetOpacity(0.1)
        f_actor.GetProperty().SetAmbient(1)
        f_actor.GetProperty().SetLineWidth(2)
        f_actor.GetProperty().RenderLinesAsTubesOn()
        f_actor.GetProperty().BackfaceCullingOff()
        f_actor.GetProperty().FrontfaceCullingOff()
        return actor, f_actor

    # Return: orbit_actors
    # Dictionary associated body name to trajectory depiction actor
    def make_orbits(data: DataLoader) -> Dict[str, vtk.vtkActor]:
        orbit_actors = {}
        sched = data.schedule
        body_groups = [ [ data.planet_names, Units.duration(days=1) ], [ data.jupiter_moon_names, Units.duration(hours=1) ] ]
        for names, period in body_groups:
            for name in names:
                myprint(verbose=verbose, text=' * {}'.format(name))
                body = data.bodies[name]
                points = []
                for t in np.arange(sched.initial_time, sched.initial_time + min(1.01*data.orbital_period[name], sched.final_time-sched.initial_time), period):
                    points.append(body.get_position(t))
                actor = VTKUtils.make_curve_actor(points, data.planet_colors[name], 2)
                actor.SetObjectName(f'orbit of {name}')
                orbit_actors[name] = actor
        myprint(verbose=verbose, text=' * Clipper')
        intp = interpolate.make_interp_spline(x=data.schedule.mission_timeline, y=data.full_paths['clipper'])
        points = []
        for t in np.arange(sched.initial_time, sched.final_time, Units.duration(hours=1)):
            points.append(intp(t))
        actor = VTKUtils.make_curve_actor(points, [0.3, 0, 0], 0.25)
        actor.SetObjectName('orbit of Clipper')
        orbit_actors['Clipper'] = actor
        return orbit_actors

    # Return vector from view point to focal point
    def focal_vector(camera: vtk.vtkCamera) -> np.ndarray:
        p = np.array(camera.GetPosition())
        f = np.array(camera.GetFocalPoint())
        return f-p

    # Print settings of provided camera
    def print_camera(camera: vtk.vtkCamera):
        settings = {
            'from' : camera.GetPosition(),
            'to'   : camera.GetFocalPoint(),
            'up'   : camera.GetViewUp(),
            'clip' : camera.GetClippingRange(),
            'view' : camera.GetViewAngle()
        }
        txt = json.dump(settings)
        print(txt)

    # Return camera imported from file
    def read_camera(name: str) -> vtk.vtkActor:
        settings = json.load(str)
        camera = vtk.vtkCamera()
        if 'from' in settings:
            camera.SetPosition(settings['from'])
        if 'to' in settings:
            camera.SetFocalPoint(settings['to'])
        if 'up' in settings:
            camera.SetViewUp(settings['up'])
        if 'clip' in settings:
            camera.SetClippingRange(settings['clip'])
        if 'view' in settings:
            camera.SetViewAngle(settings['view'])
        return camera

    # Return length of focal vector
    def focal_distance(camera: vtk.vtkCamera) -> float:
        return np.linalg.norm(VTKUtils.focal_vector(camera))

    def make_lights(actors: Dict[str, vtk.vtkActor], sun_radius: float) -> vtk.vtkLightCollection:
        print('sun radius={}'.format(sun_radius))
        print('actors = {}'.format(actors))
        collection = vtk.vtkLightCollection()
        for name in actors.keys():
            if name == 'Sun':
                continue
            a = actors[name]
            m = a.GetUserMatrix()
            p = np.array([m.GetElement(0,3), m.GetElement(1,3), m.GetElement(2,3)])
            sp = sun_radius * p/np.linalg.norm(p)
            print(f'p={p}')
            light = vtk.vtkLight()
            light.SetPosition(0.005*sp+0.995*p)
            light.SetFocalPoint(p)
            # light.SetConeAngle(180)
            light.SetPositional(0)
            light.SetIntensity(0.15)
            collection.AddItem(light)
        return collection
    
'''
SimulationState:
Class to encapsulate various aspects of the simulation state
'''
class SimulationState:
    def __init__(self, params: Dict[str, Any], time_step: float, clock: float):
        self.params = params
        if self.params.plt_metrics is not None:
            self.timesteps_plt = []
            self.clipper_europa_dist_plt = []
            self.clipper_vel_plt = []
            self.clipper_acc_plt = []
        self.time_step = time_step
        self.clock = clock
        self.tether_changed = False
        # self.focus_to_target = None
        self.anchor_changed = False 
        # self.position_to_anchor = None
        # self.focus_vector = None
        self.what_tether = 'N/A'
        self.what_anchor = 'N/A'
        self.camera_modified = False
        self.active_frame = None
        self.default_camera = {}
        self.current_camera = {}
        self.user_paused = False
        self.skip_next_update = False

'''
GraphicsObjects:
Class to encapsulate the various graphical elements of the visualization
'''
class GraphicsObjects:
    def __init__(self, window: vtk.vtkRenderWindow, all_actors: Dict[str, List[vtk.vtkActor]], resize: vtk.vtkImageResize):
        self.window = window
        self.all_actors = all_actors
        self.resize = resize
        renderers = None
        self.frame_counter = 0

'''
Overarching class that contains and controls all aspects of the simulation
and visualization
'''
class Simulation:
    def __init__(self, state: SimulationState, data: DataLoader, graphics: GraphicsObjects):
        self.state = state
        self.data = data
        self.graphics = graphics
        self.render_window = self.graphics.window
        self.state.params.record_video = False
        for name in data.all_body_names:
            name1 = ('None', name)
            name2 = ('Clipper', name)
            self.state.default_camera[name1] = VTKUtils.vantage_point(camera_name=name1, data=self.data)
            self.state.default_camera[name2] = VTKUtils.vantage_point(camera_name=name2, data=self.data)
        self.state.default_camera[('Clipper', 'Clipper')] = VTKUtils.vantage_point(('Clipper', 'Clipper'), data=self.data)
        self.state.default_camera[('Cape Canaveral', 'Sky')] = VTKUtils.vantage_point(('Cape Canaveral', 'Sky'), data=self.data)
        self.state.default_camera[('None', 'Clipper')] = VTKUtils.vantage_point(('None', 'Clipper'), data=self.data)

        for name in self.state.default_camera.keys():
            self.state.current_camera[name] = duplicate_camera(camera=self.state.default_camera[name])

        self.state.active_frame = ('None', 'Sun')

    # function to calculate the disctance between Clipper and any Celestial body defined by the User/ Programmer
    def getClipper2BodyDist(self, body=str):
        p1 = self.data.bodies[body].get_position(self.state.clock)
        p2 = self.data.clipper.get_position(self.state.clock)
        return np.linalg.norm(p1-p2)

    # function called at every step of the internal clock of the simulation
    def timer_cb(self, obj, value):
        # no state update when paused
        if self.state.params.paused:
            return
        # access current time information
        clock = self.state.clock
        timestr = VTKUtils.clock_as_str(clock)
        # no update if past end of mission time
        if clock >= self.data.schedule.final_time:
            return
        
        print('in timer callback: anchor={}, tether={}'.format(self.state.active_frame[0], self.state.active_frame[1]))
        print('active frame: {}'.format(self.state.active_frame))
        f_p = self.state.current_camera[self.state.active_frame].GetPosition()
        i_p = [ int(f_p[0]), int(f_p[1]), int(f_p[2])]
        print('active camera at: {}'.format(i_p))

        # move camera to same coordinates in updated frame 
        if not self.state.camera_modified:
            frame = self.data.dynbodyframe[self.state.active_frame]
            curcam = self.state.current_camera[self.state.active_frame]
            self.graphics.renderers['bodies'].SetActiveCamera(local2global(frame, curcam, clock))
        else:
            print('Camera modified: doing nothing')

        for name in self.data.all_body_names:
            # print(f'updating position of {name}')
            body = self.data.bodies[name]
            # Update position and rotation
            self.graphics.all_actors['bodies'][name] = body.move_actor(clock)
            self.graphics.all_actors['points'][name] = body.move_point_actor(clock)
            if name == 'Saturn':
                for ra in self.graphics.all_actors['rings']:
                    # ra.SetPosition(self.graphics.all_actors['bodies']['Saturn'].GetPosition())
                    VTKUtils.copy_transformation(body.get_actor(), ra)

        # calling the function to calculate the distance between Clipper and Europa at every time step
        dClipper2Body = int(self.getClipper2BodyDist('Europa'))
        # print('Distance from Clipper spacecraft to Europa is:', dClipper2Body, 'km')
        clipper_velocity_vec = self.data.clipper.get_velocity(self.state.clock)
        ClipperVelocity = np.linalg.norm(clipper_velocity_vec)
        # print('The Clipper velocity is', ClipperVelocity)
        # exit()
        clipper_acceleration_vec = self.data.clipper.get_acceleration(self.state.clock)
        vel_acc_proj = np.dot(clipper_velocity_vec, clipper_acceleration_vec)
        ClipperAcceleration = math.copysign(np.linalg.norm(clipper_acceleration_vec), vel_acc_proj)

        if self.state.params.plt_metrics is not None:
            self.state.timesteps_plt.append(clock)
            self.state.clipper_europa_dist_plt.append(dClipper2Body)
            self.state.clipper_vel_plt.append(ClipperVelocity)
            self.state.clipper_acc_plt.append(ClipperAcceleration)

        # extend Clipper's orbit
        self.data.clipper.update(clock)
        self.data.clipper_model_actor.SetPosition(self.data.clipper.position_actor.GetPosition())
        # rotate Clipper model to align with trajectory
        tangent = self.data.clipper.get_tangent(clock)
        print(f'tangent={tangent}')
        tangent[2] = 0
        tangent /= np.linalg.norm(tangent)
        theta = Units.rad_to_deg(math.atan2(-tangent[0], tangent[1]))
        print(f'theta={theta}')
        self.data.clipper_model_actor.SetOrientation(0, 0, theta)
        self.graphics.all_actors['bodies']['Clipper model'] = self.data.clipper_model_actor
        self.graphics.all_actors['text']['time'].SetInput(timestr)
        if self.data.clipper.has_launched == True:
            self.graphics.all_actors['text']['distance'].SetInput('Clipper-Europa Distance (km): {:,}'.format(dClipper2Body))
            self.graphics.all_actors['text']['velocity'].SetInput('Clipper Velocity (km/s): {:.3f}'.format(ClipperVelocity))
            self.graphics.all_actors['text']['launch'].SetInput('Liftoff! Clipper spacecraft has launched successfully')
            self.graphics.all_actors['text']['launch'].GetTextProperty().SetColor(0, 1.0, 0)
            self.graphics.all_actors['text']['acceleration'].SetInput('Clipper Acceleration (km/s^2): {:.7f}'.format(ClipperAcceleration))
        else:
            self.graphics.all_actors['text']['distance'].SetInput('Clipper-Europa distance Calculating...')
            self.graphics.all_actors['text']['velocity'].SetInput('Stand by for Clipper Velocity...')
            self.graphics.all_actors['text']['launch'].SetInput('Launch Sequence in progress...')
            self.graphics.all_actors['text']['launch'].GetTextProperty().SetColor(1.0, 0, 0)
            self.graphics.all_actors['text']['acceleration'].SetInput('Stand by for Clipper Acceleration...')
        self.render_window.Render()
        self.state.clock += self.state.time_step

        

'''
The central class of the program that creates all the other ones and manages user interaction
'''
class MainWindow(QMainWindow):
    def __init__(self, parent = None):
        QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.state = None
        self.data = None
        self.clipper_mission = None
        print('Main window initialized')
        self.view_angle_locked = False
        self.near_clipping_locked = False
        self.far_clipping_locked = False
    
    def change_planet_focus(self, anchor, target, ref_camera=None):

        print('anchor or target changed...')
        print('anchor={}, target={}'.format(anchor, target))
        selected_frame = None
        if anchor == 'Cape Canaveral' or target == 'Sky':
            skyid = self.ui.target_to_index('Sky')
            capeid = self.ui.anchor_to_index('Cape Canaveral')
            # force consistency between anchor and target
            self.ui.cameraTargetComboBox.blockSignals(True)
            self.ui.cameraTargetComboBox.setCurrentIndex(skyid)
            self.ui.cameraTargetComboBox.blockSignals(False)
            self.ui.cameraPositionComboBox.blockSignals(True)
            self.ui.cameraPositionComboBox.setCurrentIndex(capeid) 
            self.ui.cameraPositionComboBox.blockSignals(False)
            selected_frame = ('Cape Canaveral', 'Sky')
            self.state.what_anchor = 'Cape Canaveral'
            self.state.what_tether = 'Sky'
        else:
            selected_frame = (anchor, target)
            self.state.what_anchor = anchor
            self.state.what_tether = target

        print('selected frame = {}'.format(selected_frame))
        self.graphics.all_actors['text']['bodies'].SetInput('Current Focus: {}'.format(selected_frame[1]))
        
        if self.state.what_anchor == 'None' and self.state.what_tether == 'None':
            return
        
        self.state.active_frame=(self.state.what_anchor, self.state.what_tether)

        self.ui.viewAngleDial.blockSignals(True)
        self.ui.viewAngleDial.setValue(int(self.state.current_camera[self.state.active_frame].GetViewAngle()))
        self.ui.viewAngleDial.blockSignals(False)
        self.ui.viewAngleSpinBox.blockSignals(True)
        self.ui.viewAngleSpinBox.setValue(int(self.state.current_camera[self.state.active_frame].GetViewAngle()))
        self.ui.viewAngleSpinBox.blockSignals(False)
        self.ui.clippingNearSpinBox.blockSignals(True)
        self.ui.clippingNearSpinBox.setValue(self.state.current_camera[self.state.active_frame].GetClippingRange()[0])
        self.ui.clippingNearSpinBox.blockSignals(False)
        self.ui.clippingFarSpinBox.blockSignals(True)
        self.ui.clippingFarSpinBox.setValue(self.state.current_camera[self.state.active_frame].GetClippingRange()[1])
        self.ui.clippingFarSpinBox.blockSignals(False)

        if ref_camera is None:
            ref_camera = self.state.default_camera[selected_frame]

        if self.ui.resetCameraPosition.checkState()==2 and self.ui.resetCameraTarget.checkState()==2:
            self.state.current_camera[selected_frame] = duplicate_camera(ref_camera)
        elif self.ui.resetCameraPosition.checkState()==2:
            self.state.current_camera[selected_frame].SetPosition(ref_camera.GetPosition())
        elif self.ui.resetCameraTarget.checkState()==2:
            self.state.current_camera[selected_frame].SetFocalPoint(ref_camera.GetFocalPoint())
        # otherwise nothing to be done


        print('selected frame={}'.format(self.data.dynbodyframe[selected_frame]))
        print('selected camera={}'.format(self.state.current_camera[selected_frame]))
        f = self.data.dynbodyframe[selected_frame]
        c = self.state.current_camera[selected_frame]
        # self.graphics.renderers['bodies'].SetActiveCamera(local2global(self.data.dynbodyframe[selected_frame], self.state.current_camera[selected_frame], self.state.clock))
        self.graphics.renderers['bodies'].SetActiveCamera(local2global(f, c, self.state.clock))
        self.graphics.window.Render()

    ''' 
    Callback functions from Simulation
    '''
    def anchor_or_target_changed_cb(self, unused_id):
        anchor = self.ui.anchors[self.ui.cameraPositionComboBox.currentIndex()]
        target = self.ui.targets[self.ui.cameraTargetComboBox.currentIndex()]

        self.change_planet_focus(anchor, target)
 
    def quit_cb(self):
        if self.state.params.plt_metrics is not None:
            self.save_metrics_plts(self.state.params.plt_metrics)
        sys.exit(0)

    def anyevent_cb(self, obj, event):
        print('{}'.format(event))

    def key_pressed_cb(self, obj, event):
        new_key = obj.GetKeySym()
        print(f'key pressed: {new_key}')
        if new_key.isdigit():
            key = int(new_key)
            if key <= 8 and self.data.all_body_names[key] != self.state.params.planet_focus:
                myprint(verbose=True, text='changing focus from {} to {}'.format(
                    self.state.params.planet_focus, self.data.all_body_names[key]))
                self.state.params.planet_focus = self.data.all_body_names[key]
                self.state.planet_focus_changed = True
                self.state.params.do_tether = True
                cam_setting = VTKUtils.vantage_point(
                    target_name=self.state.params.planet_focus, data=self.data, et=self.state.clock)
                VTKUtils.get_camera(self.graphics.window).SetPosition(
                    cam_setting['pos'])
                VTKUtils.get_camera(self.graphics.window).SetFocalPoint(
                    cam_setting['focal'])
                VTKUtils.get_camera(self.graphics.window).SetViewUp(
                    cam_setting['up'])
                self.graphics.all_actors['text']['bodies'].SetInput(
                    'Current Focus: {}'.format(self.state.params.planet_focus))
        # TODO: bugged, replace with change_planet_focus()

            # Overrides the default vtk behavior for keypress '3'
            if int(new_key) == 3:
                self.graphics.window.StereoRenderOn()
        elif new_key == 'minus':
            self.state.params.planet_focus = 'N/A'
            self.state.planet_focus_changed = False
            self.graphics.all_actors['text']['bodies'].SetInput(
                'Current Focus: None')
            self.state.params.do_tether = False
        elif new_key == 'o' or new_key == 'O': # TODO: orbits bugged, circuit eventually not complete
            self.state.params.show_orbits = not self.state.params.show_orbits
            if self.state.params.show_orbits:
                myprint(verbose=True, text='now showing orbits')
                for name in self.graphics.all_actors['orbits']:
                    if self.state.clock < self.data.schedule.arrival_time and name in self.data.jupiter_moon_names:
                        continue
                    myprint(verbose=True, text='added orbit of {}'.format(name))
                    myprint(verbose=True, text='there are {} vertices in this orbit'.format(
                        self.graphics.all_actors['orbits'][name].GetMapper().GetInput().GetPoints().GetNumberOfPoints()))
                    VTKUtils.get_renderer(self.graphics.window).AddActor(
                        self.graphics.all_actors['orbits'][name])
            else:
                myprint(verbose=True, text='now hiding orbits')
                for name in self.graphics.all_actors['orbits']:
                    VTKUtils.get_renderer(self.graphics.window).Remove_actor(
                        self.graphics.all_actors['orbits'][name])
            self.graphics.window.Render()
        elif new_key == 'f':
            if self.graphics.window.GetFullScreen():
                self.graphics.window.FullScreenOff()
            self.graphics.window.Render()
        elif new_key == 'exclam' or new_key == 'at' or new_key == 'numbersign' \
                or new_key == 'dollar' or new_key == 'asciicircum' \
                or new_key == 'plus' or new_key == 'equal' or new_key == 'ampersand' \
                or new_key == 'v' or new_key == 'a' or new_key == 'b' or new_key == 'c' \
                or new_key == 'd' or new_key == 'g' or new_key == 'semicolon' \
                or new_key == 'quoteright' or new_key == 'k':
            self.state.new_focus_planet = self.state.params.planet_focus
            if new_key == 'exclam':
                self.state.new_focus_planet = 'Io'
            elif new_key == 'at':
                self.state.new_focus_planet = 'Europa'
            elif new_key == 'numbersign':
                self.state.new_focus_planet = 'Ganymede'
            elif new_key == 'dollar':
                self.state.new_focus_planet = 'Callisto'
            elif new_key == 'asciicircum':
                self.state.new_focus_planet = 'Moon'
            elif new_key == 'plus':
                self.state.new_focus_planet = 'Clipper Forward'
            elif new_key == 'equal':
                self.state.new_focus_planet = 'Clipper Backward'
            elif new_key == 'ampersand':
                self.state.new_focus_planet = 'Clipper'
            elif new_key == 'v':
                self.state.new_focus_planet = 'ecliptic'
            elif new_key == 'a':
                self.state.new_focus_planet = 'Clipper-Earth'
            elif new_key == 'b':
                self.state.new_focus_planet = 'Clipper-Mars'
            elif new_key == 'c':
                self.state.new_focus_planet = 'Clipper-Jupiter'
            elif new_key == 'd':
                self.state.new_focus_planet = 'Clipper-Europa'
            elif new_key == 'g':
                self.state.new_focus_planet = 'Clipper-Venus'
            elif new_key == 'semicolon':
                self.state.params.zoom_factor /= 2.0
            elif new_key == 'quoteright':
                self.state.params.zoom_factor *= 2.0
            elif new_key == 'k':
                VTKUtils.print_camera(self.graphics.window)
                return
            if self.state.new_focus_planet != self.state.params.planet_focus:
                self.state.params.planet_focus = self.state.new_focus_planet.split()[
                    0]
                self.state.planet_focus_changed = True
                self.state.saved_relative_position = False
                self.state.params.do_tether = True
                camera = VTKUtils.get_camera(self.graphics.window)
                cam_setting = VTKUtils.vantage_point(
                    target_name=self.state.params.planet_focus, data=self.data, et=self.state.clock)
                camera.SetPosition(cam_setting['pos'])
                camera.SetFocalPoint(cam_setting['focal'])
                camera.SetViewUp(cam_setting['up'])
                self.graphics.all_actors['text']['title'].SetInput(
                    'Current Focus: {}'.format(self.state.params.planet_focus))
            camera = VTKUtils.get_camera(self.graphics.window)
            if self.state.params.zoom_factor != 1:
                camera.Zoom(self.state.params.zoom_factor)
                self.state.params.zoom_factor = 1
                myprint(verbose=True, text='camera view angle: {}'.format(
                    camera.GetViewAngle()))
        elif new_key == 't':
            # tether mode
            self.state.params.do_tether = not self.state.params.do_tether
            if self.state.params.do_tether:
                myprint(verbose=True, text='tethering activated')
            else:
                myprint(verbose=True, text='tethering deactivated')
        elif new_key == 'Z' and self.state.params.do_tether:
            print('Zoom factor 2')
            camera = VTKUtils.get_camera(self.graphics.window)
            camera.Zoom(2)
            print(
                f'current view angle is {VTKUtils.get_camera(self.graphics.window).GetViewAngle()}')
        elif new_key == 'z' and self.state.params.do_tether:
            print('Zoom factor 1/2')
            camera = VTKUtils.get_camera(self.graphics.window)
            camera.Zoom(0.5)
            print(
                f'current view angle is {VTKUtils.get_camera(self.graphics.window).GetViewAngle()}')
        elif new_key == 'Up' or new_key == 'Down':
            cam = self.state.current_camera[self.state.active_frame]
            f = np.array(cam.GetFocalPoint())
            p = np.array(cam.GetPosition())
            direction = f - p
            if new_key == 'Up':
                p += 0.25*direction
                print(f'moving camera forward by {0.25*direction}, local position is {p}')
            else:
                p -= 0.25*direction
                print(f'moving camera backward by {0.25*direction}, local position is {p}')
            self.state.current_camera[self.state.active_frame].SetPosition(p[0], p[1], p[2])
        elif new_key == 'h':
            myprint(verbose=True, text='List of keyboard commands:')
            myprint(verbose=True, text='\'0\', \'1\', \'2\', ..., \'8\': Point the camera to the Sun, Mercury, Venus, ..., Neptune and track')
            myprint(
                verbose=True, text='\'!\', \'@\', \'#\', \'$\':  Point the camera to Jupiter\'s Galilean moons')
            myprint(
                verbose=True, text='\'^\':                     Point the camera to the Moon')
            myprint(
                verbose=True, text='\'&\':                     Point the camera to Clipper\'s position')
            myprint(verbose=True,
                    text='\'-\':                     Turn off tracking')
            myprint(
                verbose=True, text='\'f\':                     Toggle full screen on and off')
            myprint(
                verbose=True, text='\'a\', ..., \'d\':           Point camera from Clipper to Earth, Mars, Jupiter, Europa')
            myprint(
                verbose=True, text='\'Ctrl-9\':                Point camera from Clipper to the Moon')
            myprint(
                verbose=True, text='\'Ctrl-!\', \'Ctrl-@\':      Point camera from Clipper to Jupiter\'s moon')
            myprint(
                verbose=True, text='\'t\':                     Toggle tethering on and off')
            myprint(
                verbose=True, text='\'v\':                     Point the camera to the Sun\'s North pole')
            myprint(
                verbose=True, text='\'o\' (\'O\'):               Turn on/off depiction of planets\' orbits')
            myprint(
                verbose=True, text='\'Z\':                     Zoom in on current focus object')
            myprint(
                verbose=True, text='\'z\':                     Zoom out on current focus')
            myprint(verbose=True,
                    text='\'h\':                     Print this information')
        else:
            myprint(verbose=verbose,
                    text='unrecognized entered key is {}'.format(new_key))

    def time_step_changed_1minute_cb(self):
        self.time_step_change_cb(value=1)

    def time_step_changed_15minutes_cb(self):
        self.time_step_change_cb(value=15)

    def time_step_changed_30minutes_cb(self):
        self.time_step_change_cb(value=30)

    def time_step_changed_1hour_cb(self):
        self.time_step_change_cb(
            value=Units.duration(hours=1, output='minutes'))

    def time_step_changed_6hours_cb(self):
        self.time_step_change_cb(
            value=Units.duration(hours=6, output='minutes'))

    def time_step_changed_1day_cb(self):
        self.time_step_change_cb(
            value=Units.duration(days=1, output='minutes'))

    def time_step_changed_1week_cb(self):
        self.time_step_change_cb(
            value=Units.duration(weeks=1, output='minutes'))

    def time_step_changed_1month_cb(self):
        self.time_step_change_cb(
            value=Units.duration(months=1, output='minutes'))

    def date_activated_cb(self, date: QDate):
        print(f'date activated, date={date}')

    def date_page_changed_cb(self, year: int, month: int):
        print(f'date page changed: year={year}, month={month}')
        self.date_change_cb(QDate(year, month, 1))

    def date_selection_changed_cb(self):
        print('date selection changed')

    # returns time step in minutes per second
    def time_step_change_cb(self, value):
        print('time step change')
        value = Units.duration(minutes=value)
        self.state.time_step = float(value)/self.state.params.refresh_rate
        self.state.params.time_scale = float(value)
        print(f'time step is now {self.state.time_step}')

    def planet_scale_change_cb(self, value):
        print('planet scale change')
        self.state.params.planet_scale = value
        for name in self.data.bodies.keys():
            body = self.data.bodies[name]
            if name == 'Sun' and self.state.params.planet_scale > self.state.params.max_sun_scale:
                body.scale_actor(self.state.params.max_sun_scale)
            else:
                body.scale_actor(self.state.params.planet_scale)
        # adjust Saturn rings accordingly
        for actor in self.graphics.all_actors['rings']:
            actor.SetScale(self.state.params.planet_scale,
                           self.state.params.planet_scale, self.state.params.planet_scale)
        self.graphics.window.Render()

    def window_size_change_cb(self, obj, event):
        print('window size change')
        width, height = self.graphics.window.GetSize()
        if self.state.params.resolution[0] != width or self.state.params.resolution[1] != height:
            self.state.params.resolution[0] = width
            self.state.params.resolution[1] = height
            if not self.state.params.show_shadows:
                self.graphics.resize.SetOutputDimensions(
                    self.state.params.resolution[0], self.state.params.resolution[1], 0)
                self.graphics.resize.Update()
                myprint(verbose=verbose, text='texture has been resized')
            self.graphics.window.Render()

    def anyevent_cb(self, caller, event):
        print('anyevent: event={}'.format(event))
    
    def short_camera_description(self, camera, frame_name, camera2=None):
        if camera2 is None:
            p = toint(camera.GetPosition())
            f = toint(camera.GetFocalPoint())
            camg = local2global(frame=self.data.dynbodyframe[frame_name], localcam=camera, et=self.state.clock)
            pg = toint(camg.GetPosition())
            fg = toint(camg.GetFocalPoint())
            return f'local:\n\tp={p},\n\tf={f},\nglobal:\n\tp={pg},\n\tf={fg}'
        else:
            dp = np.array(camera2.GetPosition()) - np.array(camera.GetPosition())
            df = np.array(camera2.GetFocalPoint()) - np.array(camera.GetFocalPoint())
            camg1 = local2global(frame=self.data.dynbodyframe[frame_name], localcam=camera, et=self.state.clock)
            camg2 = local2global(frame=self.data.dynbodyframe[frame_name], localcam=camera2, et=self.state.clock)
            dpg = np.array(camg2.GetPosition()) - np.array(camg1.GetPosition())
            dfg = np.array(camg2.GetFocalPoint()) - np.array(camg1.GetFocalPoint())
            return f'local correction:\n\tdp=\n\t{dp}\n\tdf=\n\t{df}\n\tglobal:\n\tdp=\n\t{dpg}\n\tdf=\n\t{dfg}'
        
    def camera_change_cb(self, caller, event):
        print('in camera change callback function, event={}'.format(event))
        self.ui.clippingNearSpinBox.blockSignals(True)
        self.ui.clippingFarSpinBox.blockSignals(True)
        self.ui.viewAngleDial.blockSignals(True)
        self.ui.viewAngleSpinBox.blockSignals(True)
        camera = self.graphics.renderers['bodies'].GetActiveCamera()
        self.ui.clippingNearSpinBox.setValue(camera.GetClippingRange()[0])
        self.ui.clippingFarSpinBox.setValue(Units.km_to_au(camera.GetClippingRange()[1]))
        self.ui.viewAngleDial.setValue(int(camera.GetViewAngle()))
        self.ui.viewAngleSpinBox.setValue(int(camera.GetViewAngle()))
        self.ui.clippingNearSpinBox.blockSignals(False)
        self.ui.clippingFarSpinBox.blockSignals(False)
        self.ui.viewAngleDial.blockSignals(False)
        self.ui.viewAngleSpinBox.blockSignals(False)
        if self.state.active_frame[1] == 'None':
            print('no tethering, exiting')
            return
        # print(VTKUtils.camera_description(self.graphics.renderers['bodies'].GetActiveCamera()))
        # print('\n!!!!camera change: event is {}!!!!\n'.format(event))
        if self.state.params.paused or event=='MouseWheelBackwardEvent' or event=='MouseWheelForwardEvent':
            print('updating local camera position')
            print(f'global camera position is {self.graphics.renderers["bodies"].GetActiveCamera().GetPosition()}')
            # force expected camera behavior under user manipulation:
            # rotation: does not change the focal point or focal distance, just the location
            # forward and backward roll: does not change the focal point, but focal distance and location along focal vector
            # zooming: should only change view angle but accessible by mouse?
            frame_name = self.state.active_frame
            ref_cam = self.state.current_camera[frame_name]
            focal_vector = np.array(ref_cam.GetFocalPoint()) - np.array(ref_cam.GetPosition())
            # refcam = self.state.current_camera[frame_name]
            # old_cam_str = self.short_camera_description(self.state.current_camera[frame_name], frame_name)
            tmp_cam = global2local(frame=self.data.dynbodyframe[frame_name], globalcam=self.graphics.renderers['bodies'].GetActiveCamera(), et=self.state.clock)
            l_ref = np.linalg.norm(focal_vector)
            # print(f'before focal point correction, p= {self.state.current_camera[frame_name].GetPosition()}, f={self.state.current_camera[frame_name].GetFocalPoint()}')
            p = np.array(tmp_cam.GetPosition())
            f = np.array(tmp_cam.GetFocalPoint())
            up = np.array(tmp_cam.GetViewUp())
            far = ref_cam.GetClippingRange()[1]
            near = ref_cam.GetClippingRange()[0]
            if self.state.params.paused:
                # camera rotation assumed: preserve distance to focal point
                # first place focal point back in the center of the coordinate system
                p -= f
                f = np.array([0, 0, 0])
                l_new = np.linalg.norm(p)
                p *= l_ref/l_new
                self.state.current_camera[frame_name].SetPosition(p)
                self.state.current_camera[frame_name].SetFocalPoint(f)
                self.state.current_camera[frame_name].SetViewUp(up)
                self.state.current_camera[frame_name].SetClippingRange(near, far)
            else:
                # camera forward/backward motion assumes: preserve focal vector direction
                f = np.array([0, 0, 0])
                dp = -1.*p
                p = np.dot(-p, focal_vector)/l_ref*focal_vector
                self.state.current_camera[frame_name].SetPosition(p)
                self.state.current_camera[frame_name].SetFocalPoint(f)
                self.state.current_camera[frame_name].SetViewUp(up)
                self.state.current_camera[frame_name].SetClippingRange(near, far)
            # newcam = self.state.current_camera[frame_name]
            # save_cam = vtk.vtkCamera()
            # save_cam.DeepCopy(newcam)
            # new_cam_str = self.short_camera_description(newcam, frame_name)

            # restore focal point
            # f = np.array(refcam.GetFocalPoint())
            # print(f'focal point corrected from {toint(newcam.GetFocalPoint())} to {toint(f)}')
            # newcam.SetFocalPoint(f)
            # dnew = f - np.array(newcam.GetPosition())
            # dref = f - np.array(refcam.GetPosition())

            # if self.state.params.paused:
            #     # camera rotation assumed: restore focal distance
            #     l = np.linalg.norm(dref)
            #     pnew = f - l/np.linalg.norm(dnew)*dnew
            # else:
            #     # mouse wheel forward or backward: restore focal axis
            #     l = np.linalg.norm(dnew)
            #     pnew = f - l/np.linalg.norm(dref)*dref

            # print(f'position corrected from {toint(newcam.GetPosition())} to {toint(pnew)}')
            # newcam.SetPosition(pnew)
            # new_cam_str = self.short_camera_description(newcam, frame_name)
            # corr_str = self.short_camera_description(save_cam, frame_name, newcam)
            # self.state.current_camera[frame_name] = newcam
            # gnewcam = local2global(frame=self.data.dynbodyframe[frame_name], localcam=newcam, et=self.state.clock)
            # self.graphics.renderers['bodies'].SetActiveCamera(gnewcam)
            # # print('event={}, camera moved:\nfrom {}\nto\n{}'.format(event, old_cam_str, new_cam_str))
            # print(f'event={event}\n{corr_str}')
            # self.state.camera_changed = True
            # if event.find('MouseWheel') != -1:
            #     self.state.skip_next_update = True
        # else:
        #     print('No change to current camera setting because not paused')

        # if self.state.params.paused and self.state.what_tether != 'None':
        #     print('interaction detected while tethering, modifying camera')
        #     self.state.camera_changed = True

    def interaction_start_cb(self, caller, event):
        # print('interaction starts')
        if self.state.what_tether != 'None':
            # print('paused is set to True because tethering is active')
            self.state.params.paused = True
        # else:
        #     print('tether is False. nothing to do')

    def interaction_end_cb(self, caller, event):
        # print('interaction ends')
        if self.state.what_tether != 'None':
            self.state.params.paused = False
            self.state.skip_next_update = True
    
    def date_change(self, et):

        current = self.state.clock
        self.state.clock = et

        renderer = VTKUtils.get_renderer(self.graphics.window)

        # actors = renderer.GetActors()
        # actors.InitTraversal()
        # actor = actors.GetNextActor()
        # print('in date change: list of actors:')
        # # VTKUtils.print_actors(self.graphics.all_actors)
        # while actor is not None:
        #     print(f'{actor.GetObjectName()}')
        #     actor = actors.GetNextActor()

        if self.state.params.show_clipper_orbit:
            renderer.RemoveActor(self.graphics.all_actors['Clipper']['Clipper forward orbit'])
        # renderer.RemoveActor(
        #     self.graphics.all_actors['Clipper']['Clipper backward orbit'])
        renderer.RemoveActor(self.graphics.all_actors['Clipper']['Clipper position'])

        # print('after removal, actors are')
        # # VTKUtils.print_actors(self.graphics.all_actors)
        # actors = renderer.GetActors()
        # actors.InitTraversal()
        # actor = actors.GetNextActor()
        # while actor is not None:
        #     print(f'{actor.GetObjectName()}')
        #     actor = actors.GetNextActor()

        self.data.clipper = Clipper(self.data, self.state.clock)
        self.graphics.all_actors['Clipper'] = self.data.clipper.get_actors()
        renderer.AddActor(
            self.graphics.all_actors['Clipper']['Clipper position'])
        
        if self.state.params.show_clipper_orbit:
            renderer.AddActor(self.graphics.all_actors['Clipper']['Clipper forward orbit'])
        # renderer.AddActor(
        #     self.graphics.all_actors['Clipper']['Clipper backward orbit'])

        # print('after re-addition, actors are')
        # # VTKUtils.print_actors(self.graphics.all_actors)
        # actors = renderer.GetActors()
        # actors.InitTraversal()
        # actor = actors.GetNextActor()
        # while actor is not None:
        #     print(f'{actor.GetObjectName()}')
        #     actor = actors.GetNextActor()

    # returns date in days since start
    def date_change_cb(self, datetime: QDateTime):
        et = Units.QDateTime2time(datetime)

        self.date_change(et)
    
    def change_time_step_1minute(self):
        self.ui.timestepRadioButton[0].setChecked(True)
        self.time_step_changed_1minute_cb()
    
    def change_time_step_15minutes(self):
        self.ui.timestepRadioButton[1].setChecked(True)
        self.time_step_changed_15minutes_cb()
    
    def change_time_step_1hour(self):
        self.ui.timestepRadioButton[3].setChecked(True)
        self.time_step_changed_1hour_cb()
    
    def change_time_step_1month(self):
        self.ui.timestepRadioButton[7].setChecked(True)
        self.time_step_changed_1month_cb()

    # TODO: need UI to update clipping plane locs
    def play_event_cb(self, event_id):

        event_name = self.ui.events[event_id]

        if event_name == 'No Event':
            return

        events_time_st = {
            'launch': '2024 10, 10 15:50:00',
            'Mars assist': '2025 02, 28 00:00:00',
            'Earth assist': '2026 12, 01 12:00:00',
            'Jupiter capture': '2028 08, 01 00:00:00',
            'Europa flyby': '2031 05, 27 06:00:00'
        }
        events_time_scale = {
            'launch': self.change_time_step_1minute,
            'Mars assist': self.change_time_step_1hour,
            'Earth assist': self.change_time_step_1hour,
            'Jupiter capture': self.change_time_step_1month,
            'Europa flyby': self.change_time_step_15minutes
        }
        events_frame = {
            'launch': ('None', 'Earth'),
            'Mars assist': ('None', 'Mars'),
            'Earth assist': ('None', 'Earth'),
            'Jupiter capture': ('None', 'Jupiter'),
            'Europa flyby': ('None', 'Europa')
        }
        events_camera_pos = {
            'launch': (0, 0, 5*self.data.bodies['Earth'].get_max_radius()),
            'Mars assist': (0, 0, 100*self.data.bodies['Mars'].get_max_radius()),
            'Earth assist': (0, 0, 100*self.data.bodies['Earth'].get_max_radius()),
            'Jupiter capture': (0, 0, 8000*self.data.bodies['Jupiter'].get_max_radius()),
            'Europa flyby': (0, 0, 100*self.data.bodies['Europa'].get_max_radius())
        }
        events_camera_view_up = {
            'launch': (0, 1, 0),
            'Mars assist': (0, 1, 0),
            'Earth assist': (0, 1, 0),
            'Jupiter capture': (0, 1, 0),
            'Europa flyby': (0, 1, 0)
        }

        time_st = events_time_st[event_name]
        et_st = spice.utc2et(time_st)
        self.date_change(et_st)

        events_time_scale[event_name]()

        frame = events_frame[event_name]
        anchor = frame[0]
        target = frame[1]
        pos = events_camera_pos[event_name]
        up = events_camera_view_up[event_name]
        event_camera = VTKUtils.vantage_point(camera_name=frame, data=self.data,
                                              camera_pos=pos, camera_view_up=up)
        self.change_planet_focus(anchor, target, ref_camera=event_camera)

        # self.time_end = events_time_ed[event_name]
        # need to poll for when self.state.clock >= self.time_end

    def planet_scale_dial_changed_cb(self, value):
        self.ui.scaleSpinBox.blockSignals(True)
        self.ui.scaleSpinBox.setValue(value)
        self.ui.scaleSpinBox.blockSignals(False)
        self.planet_scale_change_cb(value)

    def planet_scale_spinbox_changed_cb(self, value):
        self.ui.scaleDial.blockSignals(True)
        self.ui.scaleDial.setValue(value)
        self.ui.scaleDial.blockSignals(False)
        self.planet_scale_change_cb(value)

    def start(self):
        self.state.params.paused = False

    def stop(self):
        self.state.params.paused = True

    def quit(self):
        if self.state.params.plt_metrics is not None:
            self.save_metrics_plts(self.state.params.plt_metrics)
        sys.exit(0)

    def timer_callback(self):
        # print('in timer callback')
        et = self.state.clock
        cam = self.state.current_camera[self.state.active_frame]
        print(f'In timer, camera position is {cam.GetPosition()},\nfocal point is {cam.GetFocalPoint()},\nUp vector is {cam.GetViewUp()}\nat t={self.state.clock}')
        # today = Units.time2QDate(et)
        now = Units.time2QDateTime(et)
        self.ui.calendar.blockSignals(True)
        self.ui.calendar.setDate(now.date())
        self.ui.calendar.setTime(now.time())
        self.ui.calendar.blockSignals(False)
        if self.state.params.record_video:
            self.win2img.Modified()
            self.png_writer.SetFileName(f'{self.recording_basename}{int(et)}.png')
            self.png_writer.Write()
        for key in self.data.dynbodyframe.keys():
            if isinstance(self.graphics.renderers['bodies'], Frame):
                self.graphics.renderers['bodies'].RemoveActor(self.data.dynbodyframe[key].actor)
                self.graphics.renderers['bodies'].AddActor(self.data.dynbodyframe[key].get_actor(et))
        if not self.state.params.paused and not self.state.user_paused and not self.state.skip_next_update:
            self.clipper_mission.timer_cb(obj=None, value=self.state.clock)
        elif self.state.skip_next_update:
            self.state.skip_next_update = False

    def screenshot_cb(self):
        # ---------------------------------------------------------------
        # Save current contents of render window to PNG file
        # ---------------------------------------------------------------
        file_name = self.state.params.video_name + \
            str(self.state.clock) + ".png"
        image = vtk.vtkWindowToImageFilter()
        image.SetInput(self.graphics.window)
        png_writer = vtk.vtkPNGWriter()
        png_writer.SetInputConnection(image.GetOutputPort())
        png_writer.SetFileName(file_name)
        self.graphics.window.Render()
        png_writer.Write()

    def nowAsString(self):
        now = datetime.now()
        return now.strftime("%Y_%m_%d_%H.%M.%S")

    def start_recording_cb(self):
        self.state.params.record_video = True
        # open a new file output file for movie export
        self.win2img = vtk.vtkWindowToImageFilter()
        self.win2img.SetInput(self.graphics.window)
        self.win2img.SetInputBufferTypeToRGBA()
        self.win2img.ReadFrontBufferOff()
        self.win2img.Update()  # needed, why?
        name = self.state.params.video_name + '-recorded_at=' + self.nowAsString() + '-start_et=' + VTKUtils.clock_as_str(self.state.clock).replace(' ', '_')
        os.mkdir(name)
        self.recording_basename = os.path.join(name, 'frame_')
        self.png_writer = vtk.vtkPNGWriter()
        self.png_writer.SetInputConnection(self.win2img.GetOutputPort())
        print('recording basename =' +
              f'{self.recording_basename}{int(self.state.clock)}.png')
        print(f'directory name = {name}')
        self.png_writer.SetFileName(f'{self.recording_basename}{int(self.state.clock)}.png')
        self.png_writer.Write()

    def stop_recording_cb(self):
        self.state.params.record = False

    def pause_recording_cb(self):
        self.state.params.record = False

    def pause_clock_cb(self):
        if self.state.user_paused is False:
            self.state.user_paused = True
            self.ui.pausePushButton.setDown(True)
        else:
            self.state.user_paused = False
            self.ui.pausePushButton.setDown(False)

    def camera_view_angle_changed_cb(self, value):
        print('changing camera view angle')
        # enforce consistency between dial and spinbox
        self.ui.viewAngleSpinBox.blockSignals(True)
        self.ui.viewAngleDial.blockSignals(True)
        self.ui.viewAngleSpinBox.setValue(value)
        self.ui.viewAngleDial.setValue(value)
        self.ui.viewAngleSpinBox.blockSignals(False)
        self.ui.viewAngleDial.blockSignals(False)
        self.state.current_camera[self.state.active_frame].SetViewAngle(value)
        self.state.camera_modified = False

    def camera_near_clipping_changed_cb(self, value):
        amax = self.state.current_camera[self.state.active_frame].GetClippingRange()[1]
        self.state.current_camera[self.state.active_frame].SetClippingRange(value, amax)
        self.state.camera_modified = False
    
    def camera_far_clipping_changed_cb(self, value):
        amin = self.state.current_camera[self.state.active_frame].GetClippingRange()[0]
        amax = Units.au_to_km(value)
        self.state.current_camera[self.state.active_frame].SetClippingRange(amin, amax)
        self.state.camera_modified = False

    def save_metrics_plts(self, plt_metrics_path):
        tsteps = self.state.timesteps_plt

        plt.figure()
        plt.plot(tsteps, self.state.clipper_europa_dist_plt, color=(1.0, 1.0, 0))
        plt.xlabel("clock time (s)")
        plt.ylabel("Clipper-Europa distance (km)")
        plt.title("Clipper-Europa Distance as a Function of Time")
        plt.savefig(os.path.join(plt_metrics_path, "clipper_europa_dist.png"))

        plt.figure()
        plt.plot(tsteps, self.state.clipper_vel_plt, color=(0.016, 0.85, 1.0))
        plt.xlabel("clock time (s)")
        plt.ylabel("Clipper velocity (km/s)")
        plt.title("Clipper Velocity as a Function of Time")
        plt.savefig(os.path.join(plt_metrics_path, "clipper_vel.png"))

        plt.figure()
        plt.plot(tsteps, self.state.clipper_acc_plt, color=(1, 0.063, 0.96))
        plt.ticklabel_format(axis='y', style='sci')
        plt.xlabel("clock time (s)")
        plt.ylabel("Clipper acceleration (km/s$^2$)")
        plt.title("Clipper Acceleration as a Function of Time")
        plt.savefig(os.path.join(plt_metrics_path, "clipper_acc.png"))

        fig, ax1 = plt.subplots()
        l1, = ax1.plot(tsteps, self.state.clipper_acc_plt, color=(1, 0.063, 0.96))
        ax1.set_xlabel("clock time (s)")
        ax1.set_ylabel("unitless for comparison (scaled to acceleration)")
        ax1.set_title("Combined Metric Plots as a Function of Time")
        ax2 = ax1.twinx()
        l2, = ax2.plot(tsteps, self.state.clipper_europa_dist_plt, color=(1.0, 1.0, 0))
        ax2.set_yticks([])
        ax3 = ax1.twinx()
        l3, = ax3.plot(tsteps, self.state.clipper_vel_plt, color=(0.016, 0.85, 1.0))
        ax1.legend([l2, l3, l1], ["dist", "vel", "acc"])
        ax3.set_yticks([])

        fig.savefig(os.path.join(plt_metrics_path, "clipper_combined.png"))

    def main(self, args):
        verbose = args.verbose
        ''' Mission Data '''
        print('creating mission data')
        self.data = DataLoader(media_path=args.media_path, naif_path=args.naif_path, scale=args.scale)
        print('mission data created')

        # Timer Frequency
        timer_period = 1./args.refresh_rate
        args.time_step = timer_period * Units.duration(minutes=1)
        args.time_scale = args.time_step * args.refresh_rate
        args.zoom_factor = 1

        # global state variables
        args.planet_scale = 1
        args.ambient = 0.4
        args.paused = False

        # interaction state variables
        args.planet_focus_changed = False
        args.do_tether = False
        args.saved_relative_position = None
        args.in_progress = False
        args.cam_changed = True

        # Render size set to HD Display

        ''' Mission State '''
        print('creating mission state')
        self.state = SimulationState(params=args, time_step=args.time_step, clock=self.data.schedule.initial_time)
        print('mission state created')

        myprint(verbose=verbose, text='time step is {}'.format(args.time_step))

        # Create render window
        render_window = self.ui.vtkWidget.GetRenderWindow()
        if not args.show_shadows:
            render_window.SetNumberOfLayers(2)
        if args.resolution[0] > 0 and args.resolution[1] > 0:
            render_window.SetSize(args.resolution[0], args.resolution[1])  # Set the window size
        else:
            resolution = render_window.GetScreenSize()
            render_window.SetSize(resolution)
        render_window.SetWindowName("Solar System")
        render_window.StereoRenderOff()

        all_actors = {}
        myprint(verbose=verbose, text='Creating Jupiter\'s rings...')
        all_actors['rings'] = VTKUtils.make_rings(self.data)
        myprint(verbose=verbose, text='...done.')

        myprint(verbose=verbose, text='Creating ecliptic plane...')
        all_actors['ecliptic'] = VTKUtils.make_ecliptic_plane(self.data)
        myprint(verbose=verbose, text='...done.')

        myprint(verbose=verbose, text='Creating bodies\' orbit curves...')
        all_actors['orbits'] = VTKUtils.make_orbits(self.data)
        myprint(verbose=verbose, text='...done.')

        myprint(verbose=verbose, text='Creating Clipper object...')
        all_actors['Clipper'] = self.data.clipper.get_actors()
        all_actors['Clipper model'] = self.data.clipper_model_actor
        all_actors['Clipper model'].SetPosition(self.data.clipper.get_position(self.state.clock))

        idx = 0
        myprint(verbose=verbose, text='Initializing VTK actors...')
        body_actors = {}
        body_point_actors = {}
        for name in self.data.all_body_names:
            print(f'processing {name}')
            body = self.data.bodies[name]

            # Create Actors and add textures, rotations, and orbit position
            body_actors[name] = body.move_actor(self.state.clock)
            body_actors[name].GetProperty().SetAmbient(self.state.params.ambient)
            body_actors[name].GetProperty().SetSpecular(0) #.01)
            body_actors[name].GetProperty().SetDiffuse(0.65)

            body_point_actors[name] = body.move_point_actor(self.state.clock)

            if name == 'Saturn':
                myprint(verbose=verbose, text='moving rings')
                for actor in all_actors['rings']:
                    VTKUtils.copy_transformation(body.get_actor(), actor)
                    actor.GetProperty().SetAmbient(self.state.params.ambient)
        myprint(verbose=verbose, text='...done.')
        all_actors['bodies'] = body_actors
        all_actors['points'] = body_point_actors
        resize_image = vtk.vtkImageResize()

        print('after creating all actors, we have:\n{}'.format(all_actors))

        ''' Mission Graphics '''
        self.graphics = GraphicsObjects(window=render_window, all_actors=all_actors, resize=resize_image)

        ''' Mission '''
        self.clipper_mission = Simulation(state=self.state, data=self.data, graphics=self.graphics)

        # Create Lighting
        sunlight1 = vtk.vtkLight()
        sunlight1.SetPosition(0, 0, 0)
        sunlight1.SetFocalPoint(1.0e+8, 0, 0)  # Right hemisphere lighting
        sunlight1.SetConeAngle(180)
        sunlight1.SetPositional(1)

        sunlight2 = vtk.vtkLight()
        sunlight2.SetPosition(0, 0, 0)
        sunlight2.SetFocalPoint(-1.0e+8, 0, 0) # Left hemisphere lighting
        sunlight2.SetConeAngle(180)
        sunlight2.SetPositional(1)

        # Add Background Starry Image
        if not args.show_shadows:
            myprint(verbose=verbose, text='Loading background image...')
            jpeg_reader = vtk.vtkJPEGReader()
            jpeg_reader.SetFileName(self.data.texture_files['Stars'])
            jpeg_reader.Update()

            # Resize Image to Match Render Window Size
            resize_image.SetInputData(jpeg_reader.GetOutput())
            resize_image.SetResizeMethod(0)
            if args.resolution[0] > 0 and args.resolution[1] > 0:
                resize_image.SetOutputDimensions(args.resolution[0], args.resolution[1], 0)
            resize_image.Update()

            background_mapper = vtk.vtkImageMapper()
            background_mapper.SetInputData(resize_image.GetOutput())
            background_actor = vtk.vtkActor2D()
            background_actor.SetObjectName('background')
            background_actor.SetMapper(background_mapper)
            all_actors['background'] = background_actor
            myprint(verbose=verbose, text='...done')

        # Light up the sun no matter what
        sun_props = vtk.vtkProperty()
        sun_props.SetLighting(0)
        body_actors['Sun'].SetProperty(sun_props)

        # Create a renderer and add the actors to the scene
        myprint(verbose=verbose, text='Rendering Planets...')
        renderers = {}
        renderers['bodies'] = vtk.vtkRenderer()
        renderers['bodies'].SetBackground(0, 0, 0)  # black background
        renderers['bodies'].AddLight(sunlight1)
        renderers['bodies'].AddLight(sunlight2)
        myprint(verbose=verbose, text='Adding all actors to renderer...')
        for name in self.data.all_body_names:
            renderers['bodies'].AddActor(body_actors[name])
            renderers['bodies'].AddActor(body_point_actors[name])
        for actor in all_actors['rings']:
            renderers['bodies'].AddActor(actor)
        if args.show_ecliptic:
            for actor in self.graphics.ecliptic_actors:
                renderers['bodies'].AddActor(actor)
        if args.show_orbits:
            myprint(verbose=verbose, text='orbit_actors={}'.format(self.graphics.all_actors['orbits']))
            for body, actor in all_actors['bodies'].items():
                myprint(verbose=verbose, text='body={}, actor={}'.format(body, actor))
                renderers['bodies'].AddActor(actor)
        if args.show_clipper_orbit:
            renderers['bodies'].AddActor(all_actors['Clipper']['Clipper forward orbit'])
            # renderers['bodies'].AddActor(all_actors['Clipper']['Clipper backward orbit'])
        renderers['bodies'].AddActor(all_actors['Clipper']['Clipper position'])
        renderers['bodies'].AddActor(all_actors['Clipper model'])
        myprint(verbose=verbose, text='...done.')

        # Create the camera focused on the Sun
        camera = local2global(frame=self.data.dynbodyframe[('None','Sun')], localcam=self.state.default_camera[('None','Sun')], et=self.data.schedule.initial_time)
        self.state.active_frame = ('None', 'Sun')
        self.ui.cameraPositionComboBox.setCurrentIndex(0)
        self.ui.cameraTargetComboBox.setCurrentIndex(self.ui.target_to_index('Sun'))
        renderers['bodies'].SetActiveCamera(camera)
        if not args.show_shadows:
            renderers['bodies'].SetLayer(1)  # In front of background image

        if args.show_shadows:
            myprint(verbose=verbose, text='adding shadows')
            render_window.SetMultiSamples(0)
            shadows = vtk.vtkShadowMapPass()
            sequence = vtk.vtkSequencePass()
            passes = vtk.vtkRenderPassCollection()
            passes.AddItem(shadows.GetShadowMapBakerPass())
            passes.AddItem(shadows)
            sequence.SetPasses(passes)
            camera_pass = vtk.vtkCameraPass()
            camera_pass.SetDelegatePass(sequence)

            renderers['bodies'].SetPass(camera_pass)


        # Add useful text to the renderer
        # Render Title
        text_actors = {}
        myprint(verbose=verbose, text='Setting text actors...')
        print(f'resolution={args.resolution}')
        text_actors['title'] = vtk.vtkTextActor()
        text_actors['title'].SetInput("The Solar System")
        text_actors['title'].SetPosition(40, args.resolution[1] - 90)
        text_actors['title'].GetTextProperty().SetFontSize(15)
        text_actors['title'].GetTextProperty().SetColor(1.0, 1.0, 1.0)

        # Current planet of focus
        text_actors['bodies'] = vtk.vtkTextActor()
        text_actors['bodies'].SetInput('Current Focus: Sun')
        text_actors['bodies'].GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['bodies'].GetPositionCoordinate().SetValue(.80, .95)
        text_actors['bodies'].GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['bodies'].GetPosition2Coordinate().SetValue(.99, .95)
        text_actors['bodies'].GetTextProperty().SetFontSize(10)
        text_actors['bodies'].GetTextProperty().SetColor(1.0, 1.0, 1.0)

        text_actors['time'] = vtk.vtkTextActor()
        text_actors['time'].SetInput('{}'.format(spice.timout(self.data.schedule.initial_time, 'AP:MN:SC AMPM Month DD, YYYY')))
        text_actors['time'].GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['time'].GetPositionCoordinate().SetValue(0.01, 0.95)
        text_actors['time'].GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['time'].GetPosition2Coordinate().SetValue(0.25, 0.95)
        text_actors['time'].GetTextProperty().SetFontSize(10)
        text_actors['time'].GetTextProperty().SetColor(1.0, 1.0, 1.0)

        # Adding the distance visualization (Clipper-Europa distance calculations)
        text_actors['distance'] = vtk.vtkTextActor()
        text_actors['distance'].SetInput('Clipper-Europa distance Calculating...')
        text_actors['distance'].GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['distance'].GetPositionCoordinate().SetValue(0.01, 0.75)
        text_actors['distance'].GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['distance'].GetPosition2Coordinate().SetValue(0.25, 0.75)
        text_actors['distance'].GetTextProperty().SetFontSize(12)
        text_actors['distance'].GetTextProperty().SetColor(1.0, 1.0, 0) 

        # Adding the Clipper velocity visualization
        text_actors['velocity'] = vtk.vtkTextActor()
        text_actors['velocity'].SetInput('Stand by for Clipper Velocity...')
        text_actors['velocity'].GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['velocity'].GetPositionCoordinate().SetValue(0.01, 0.70)
        text_actors['velocity'].GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['velocity'].GetPosition2Coordinate().SetValue(0.25, 0.70)
        text_actors['velocity'].GetTextProperty().SetFontSize(12)
        text_actors['velocity'].GetTextProperty().SetColor(0.016, 0.85, 1.0)
        
         # Adding the Clipper acceleration visualization
        text_actors['acceleration'] = vtk.vtkTextActor()
        text_actors['acceleration'].SetInput('Stand by for Clipper Acceleration...')
        text_actors['acceleration'].GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['acceleration'].GetPositionCoordinate().SetValue(0.01, 0.65)
        text_actors['acceleration'].GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['acceleration'].GetPosition2Coordinate().SetValue(0.25, 0.65)
        text_actors['acceleration'].GetTextProperty().SetFontSize(12)
        text_actors['acceleration'].GetTextProperty().SetColor(1, 0.063, 0.96)

        # Adding a sign for users to know when Clipper has launched
        text_actors['launch'] = vtk.vtkTextActor()
        text_actors['launch'].SetInput('Launch Sequence in progress...')
        text_actors['launch'].GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['launch'].GetPositionCoordinate().SetValue(0.01, 0.85)
        text_actors['launch'].GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['launch'].GetPosition2Coordinate().SetValue(0.25, 0.85)
        text_actors['launch'].GetTextProperty().SetFontSize(12)
        text_actors['launch'].GetTextProperty().SetColor(1.0, 0, 0)
        



        all_actors['text'] = text_actors
        myprint(verbose=verbose, text='...done')

        renderers['bodies'].AddActor2D(text_actors['time'])
        renderers['bodies'].AddActor2D(text_actors['bodies'])
        renderers['bodies'].AddActor2D(text_actors['distance'])
        renderers['bodies'].AddActor2D(text_actors['velocity'])
        renderers['bodies'].AddActor2D(text_actors['launch'])
        renderers['bodies'].AddActor2D(text_actors['acceleration'])

        render_window.AddRenderer(renderers['bodies'])

        # Create a background renderer
        if not args.show_shadows:
            background_renderer = vtk.vtkRenderer()
            background_renderer.SetLayer(0)
            background_renderer.InteractiveOff()
            background_renderer.AddActor(background_actor)
            render_window.AddRenderer(background_renderer)

        # visualize ecliptic plane
        if args.show_ecliptic:
            for actor in all_actors['ecliptic']:
                renderers['bodies'].AddActor(actor)

        self.clipper_mission.graphics.renderers = renderers

        # Set-up interactor
        render_window_interactor = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        style = vtk.vtkInteractorStyleTrackballCamera()
        style.AutoAdjustCameraClippingRangeOff()
        render_window_interactor.SetInteractorStyle(style)
        render_window_interactor.Initialize()
        render_window.Render()
        render_window.GetRenderers().GetFirstRenderer().GetActiveCamera().AddObserver('AnyEvent', self.camera_change_cb)

        render_window.AddObserver('WindowResizeEvent', self.window_size_change_cb)

        myprint(verbose=verbose, text='Creating sliders...')

        # Setting up widgets
        # space control column:
        self.ui.cameraTargetComboBox.setCurrentIndex(0)
        self.ui.cameraPositionComboBox.setCurrentIndex(0)
        self.ui.scaleDial.setMinimum(1)
        self.ui.scaleDial.setMaximum(50)
        self.ui.scaleDial.setSliderPosition(1)
        self.ui.scaleSpinBox.setValue(1)
        self.ui.cameraPositionComboBox.currentIndexChanged.connect(self.anchor_or_target_changed_cb)
        self.ui.cameraTargetComboBox.currentIndexChanged.connect(self.anchor_or_target_changed_cb)
        self.ui.scaleDial.valueChanged.connect(self.planet_scale_dial_changed_cb)
        self.ui.scaleSpinBox.valueChanged.connect(self.planet_scale_spinbox_changed_cb)
        
        # time control column
        self.ui.pausePushButton.pressed.connect(self.pause_clock_cb)
        self.ui.timestepRadioButton[0].pressed.connect(self.time_step_changed_1minute_cb)
        self.ui.timestepRadioButton[1].pressed.connect(self.time_step_changed_15minutes_cb)
        self.ui.timestepRadioButton[2].pressed.connect(self.time_step_changed_30minutes_cb)
        self.ui.timestepRadioButton[3].pressed.connect(self.time_step_changed_1hour_cb)
        self.ui.timestepRadioButton[4].pressed.connect(self.time_step_changed_6hours_cb)
        self.ui.timestepRadioButton[5].pressed.connect(self.time_step_changed_1day_cb)
        self.ui.timestepRadioButton[6].pressed.connect(self.time_step_changed_1week_cb)
        self.ui.timestepRadioButton[7].pressed.connect(self.time_step_changed_1month_cb)

        self.ui.timestepRadioButton[0].setChecked(True)

        # calendar column
        self.ui.calendar.dateTimeChanged.connect(self.date_change_cb)
        # self.ui.calendar.activated.connect(self.date_activated_cb)
        # self.ui.calendar.currentPageChanged.connect(self.date_page_changed_cb)
        # self.ui.calendar.selectionChanged.connect(self.date_selection_changed_cb)

        self.ui.eventComboBox.currentIndexChanged.connect(self.play_event_cb)

        # camera control
        self.ui.viewAngleDial.setMinimum(5)
        self.ui.viewAngleDial.setMaximum(120)
        self.ui.viewAngleDial.setSliderPosition(45)
        self.ui.viewAngleDial.valueChanged.connect(self.camera_view_angle_changed_cb)
        self.ui.viewAngleSpinBox.setMinimum(5)
        self.ui.viewAngleSpinBox.setMaximum(120)
        self.ui.viewAngleSpinBox.setValue(45)
        self.ui.viewAngleSpinBox.valueChanged.connect(self.camera_view_angle_changed_cb)
        self.ui.clippingNearSpinBox.setMinimum(0.0001)
        self.ui.clippingNearSpinBox.setMaximum(10000.0)
        self.ui.clippingNearSpinBox.setValue(1)
        self.ui.clippingNearSpinBox.setDecimals(4)
        self.ui.clippingNearSpinBox.valueChanged.connect(self.camera_near_clipping_changed_cb)
        self.ui.clippingFarSpinBox.setMinimum(1.0)
        self.ui.clippingFarSpinBox.setMaximum(100.0)
        self.ui.clippingFarSpinBox.setValue(10)
        self.ui.clippingFarSpinBox.valueChanged.connect(self.camera_far_clipping_changed_cb)

        # recording control column
        self.ui.startRecordingPushButton.pressed.connect(self.start_recording_cb)
        self.ui.pauseRecordingPushButton.pressed.connect(self.pause_recording_cb)
        self.ui.stopRecordingPushButton.pressed.connect(self.stop_recording_cb)
        self.ui.quitPushButton.pressed.connect(self.quit_cb)

        final_day = int(what_day(data=self.data, et=self.data.schedule.final_time))
        start_day = int(what_day(data=self.data, et=self.data.schedule.initial_time))
        print(f'clock start: {start_day}, clock end: {final_day}, duration: {final_day-start_day} days')
    
        myprint(verbose=verbose, text='Initializing interactor and adding observers')
        interactor = render_window.GetInteractor()
        interactor.AddObserver('KeyPressEvent', self.key_pressed_cb)
        # catch all camera manipulation events
        interactor.AddObserver('StartInteractionEvent', self.interaction_start_cb)
        interactor.AddObserver('EndInteractionEvent', self.interaction_end_cb)
        interactor.AddObserver('InteractionEvent', self.camera_change_cb)
        interactor.AddObserver('MouseWheelForwardEvent', self.camera_change_cb)
        interactor.AddObserver('MouseWheelBackwardEvent', self.camera_change_cb)
        self.timer = QTimer(self)
        print(f'time step is {args.time_step} s.')

        myprint(verbose=verbose, text='...done')
        render_window.Render()
        myprint(verbose=verbose, text='Starting timer')
        self.timer.timeout.connect(self.timer_callback)
        self.timer.start(int(1000.*timer_period))
        myprint(verbose=verbose, text='...done')

        myprint(verbose=verbose, text='Starting app...')
        app_return_code = app.exec_()
        if self.state.params.plt_metrics is not None:
            self.save_metrics_plts(self.state.params.plt_metrics)
        sys.exit(app_return_code)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Visualize Clipper Spacecraft and Solar System using NAIF Data',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--show_ecliptic', action='store_true', help='Show ecliptic plane')
    parser.add_argument('--show_clipper_orbit', action='store_true', help='Display Clipper trajectory as it progresses')
    parser.add_argument('--show_orbits', action='store_true', help='Show planets\' orbits')
    parser.add_argument('-r', '--resolution', metavar='int', default=[2560, 1440], type=int, nargs=2, help='window resolution (-1: full screen)')
    parser.add_argument('--max_sun_scale', metavar='float', default='40', type=float, help='Max sun radius (in multiples of Mercury semi-major axis)')
    parser.add_argument('--max_planet_scale', metavar='float', default=1000.,
                        type=float, nargs=2, help='Range of scaling factors for planets')
    parser.add_argument('--planet_scale', metavar='float', default=1.,
                        type=float, nargs=2, help='Initial scaling factors for planets')
    parser.add_argument('--planet_focus', metavar='body', default='Sun', help='Initial body focus')
    parser.add_argument('-x', '--scale', metavar='float', default=1, type=float, help='Uniform space scaling factor')
    parser.add_argument('--time_range', metavar='float', default=[-1, -1], type=float, nargs=2, help='Temporal range to consider in days (using mission start as reference)')
    parser.add_argument('--video_name', metavar='filename', default='ephemeris', type=str, help='Filename to store each individual frame of the animation')
    parser.add_argument('-c', '--camera', metavar='json', type=str, help='Json file containing camera setting')
    parser.add_argument('--refresh_rate', type=float, default=5.0, help='Simulation refresh rate')
    parser.add_argument('--show_shadows', action='store_true', help='Use shadow map pass')
    parser.add_argument('--verbose', action='store_true', help='Toggle verbose mode')
    parser.add_argument('--media_path', type=str, default='media', help='Path to media files')
    parser.add_argument('--naif_path', type=str, default='naif', help='Path to naif files')
    parser.add_argument('--plt_metrics', type=str, nargs='?', const='', help='directory to plot metrics from simulation run')
    args = parser.parse_args()

    app = QApplication(sys.argv)
    window = MainWindow()
    window.setWindowState(Qt.WindowMaximized)
    window.show()
    window.main(args)
