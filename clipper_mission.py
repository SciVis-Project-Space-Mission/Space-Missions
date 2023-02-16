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

'''
Note to self:
    - time_scale is in hours per second (h/s)
    - timeline position is in days (d)
    - refresh rate (in Hz) is fixed (Hz=1/s)
    - time step (in hours) is time scale / refresh rate 
'''

from PyQt5.QtWidgets import QApplication, QWidget, QMainWindow, QSlider, QGridLayout, QLabel, QPushButton, QTextEdit, QDoubleSpinBox, QDial, QProgressBar
import PyQt5.QtCore as QtCore
from PyQt5.QtCore import Qt, QTimer
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

verbose = False

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName('The Main Window')
        MainWindow.setWindowTitle('Clipper Mission')
        # in Qt, windows are made of widgets.
        # centralWidget will contains all the other widgets
        self.centralWidget = QWidget(MainWindow)
        # we will organize the contents of our centralWidget
        # in a grid / table layout
        self.gridlayout = QGridLayout(self.centralWidget)
        # vtkWidget is a widget that encapsulates a vtkRenderWindow
        # and the associated vtkRenderWindowInteractor. We add
        # it to centralWidget.
        self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
        # Sliders
        self.timescale_spinbox = QDoubleSpinBox()
        self.date_slider = QSlider()
        self.planetscale_dial = QDial()
        # Push buttons
        self.push_screenshot = QPushButton()
        self.push_screenshot.setText('Save screenshot')
        self.push_quit = QPushButton()
        self.push_quit.setText('Quit')
        '''
        [----------][----------][----------][----------][----------][----------][----------][----------][-----------]

            0              1          2           3          4           5            6          7          8
        [Time scale][----------][----------][Planet scale][----------][--Timeline][----------][---------][-snapshot-]  
        [QSpingBox-][----slider--][---value--][--------slider--------][--value--][--quit----] 
        '''
        self.gridlayout.addWidget(self.vtkWidget, 0, 0, 5, 8)
        self.gridlayout.addWidget(QLabel("Time Scale (hours/second)"), 5, 0, 1, 1)
        self.gridlayout.addWidget(self.timescale_spinbox, 5, 1, 1, 1)
        self.gridlayout.addWidget(QLabel("Planet Scale"), 5, 3, 1, 1)
        self.gridlayout.addWidget(self.planetscale_dial, 5, 4, 1, 1)
        self.gridlayout.addWidget(QLabel("Date (days)"), 5, 5, 1, 1)
        self.gridlayout.addWidget(self.date_slider, 5, 6, 1, 2)
        self.gridlayout.addWidget(self.push_screenshot, 6, 0, 1, 1)
        self.gridlayout.addWidget(self.push_quit, 6, 1, 1, 1)
        MainWindow.setCentralWidget(self.centralWidget)

class Units:
    # set up constants
    day = 86400
    hour = 3600
    minute = 60
    second = 1

    # handy unit conversions
    def au_to_km(au: int) -> int:
        return au * 149597871

    def km_to_au(km: int) -> int:
        return km / 149597871

    def deg_to_rad(deg: float) -> float:
        return deg/180.*math.pi

    def to_unit_rgb(rgb: List[int]):
        r = float(rgb[0]/255.)
        g = float(rgb[1]/255.)
        b = float(rgb[2]/255.)
        return [r, g, b]

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

class Body(object):
    def __init__(self, name: str, texture_name: str, meridian: float, radii: List[float], color: List[int]):
        self.name = name.split()[0]
        self.img_filename = texture_name #texture_files[self.name]
        sphere = vtk.vtkTexturedSphereSource()
        sphere.SetRadius(1)
        sphere.SetThetaResolution(40)
        sphere.SetPhiResolution(40)

        shape = vtk.vtkTransform()
        # shape.RotateZ(meridian) #-meridian_rotation[self.name])
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
            self.point_actor.GetProperty().RenderPointsAsSpheresOn()
            self.point_actor.GetProperty().SetColor(color[0], color[1], color[2])
        else:
            self.point_actor = vtk.vtkActor()
        if name == 'Sun':
            self.point_actor.GetProperty().SetAmbient(1000)

        self.local_frame = 'IAU_' + self.name.upper()
        self.states = None
        self.times = None
        self.actor = None
        self.orbit = None
        self.orbit_frame = None

    def get_max_radius(self) -> float:
        return np.amax(self.radii)

    def set_states(self, states, ets, frame):
        self.states = states
        self.times = ets
        self.orbit_frame = frame
        self.orbit = interpolate.make_interp_spline(self.times, self.states)

    def get_actor(self) -> vtk.vtkActor:
        # Body in its own reference frame
        if self.actor is None:
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(self.ellipsoid)
            mapper.ScalarVisibilityOff()
            self.actor = vtk.vtkActor()
            self.actor.SetMapper(mapper)
            reader = vtk.vtkJPEGReader()
            reader.SetFileName(self.img_filename)
            self.texture = vtk.vtkTexture()
            self.texture.SetInputConnection(reader.GetOutputPort())
            self.actor.SetTexture(self.texture)
        return self.actor

    def scale_actor(self, scaling: float):
        self.actor.SetScale(scaling, scaling, scaling)

    def get_position(self, et: float, frame=default_frame) -> List[float]:
        if self.orbit is not None:
            if self.orbit_frame == frame:
                return self.orbit(et)
            else:
                myprint(verbose=verbose, text='doing coordinate transformation')
                return DataLoader.coordinate_transform(self.orbit(et), self.orbit_frame, frame, et)
        else:
            myprint(verbose=verbose, text='Invalid request: get_position: no states available')
        return [0.,0.,0.]

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

    def get_local_position(self, pos: List[float], et: float) -> List[float]:
        mat = self.get_rotation(et)
        center =self.get_position(et)
        dp = pos-center
        d = np.linalg.norm(dp)
        e1 = mat[:,0]
        e2 = mat[:,1]
        e3 = mat[:,2]
        x = np.dot(dp, e1)
        y = np.dot(dp, e2)
        z = np.dot(dp, e3)
        return [x, y, z]

    def get_equatorial_plane(self, et: float, frame=default_frame) -> List[float]:
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

class Schedule:
    def __init__(self, init: float, arrival: float, end: float):
        self.initial_time = init
        self.launch_time = init + 1.5*Units.day
        self.arrival_time = arrival
        self.final_time = end
        # set up time array
        # total duration of the covered portion of Clipper in days
        self.full_timeline = np.arange(
            self.initial_time, self.final_time, Units.day)
        # 1.5 day from launch in minutes
        self.prelaunch_timeline = self.initial_time + \
            np.arange(0, Units.day*1.5, Units.minute)
        # interplanetary travel in days (includes time to launch)
        self.pretour_timeline = np.arange(
            self.initial_time, self.arrival_time, Units.day)
        # moon tour in hours
        self.tour_timeline = np.arange(self.arrival_time, self.final_time, Units.hour)
        self.transfer_timeline = np.arange(
            self.launch_time, self.arrival_time, Units.day)
        self.mission_timeline= np.concatenate((self.prelaunch_timeline, self.transfer_timeline, self.tour_timeline))

class Clipper:
    def __init__(self, data, clock=None):
        self.earth = data.bodies['Earth']
        self.schedule = data.schedule
        self.has_launched = False
        # print(f'mission timeline = {self.schedule.mission_timeline}')
        # print(f'clipper steps = {data.full_paths["clipper"]}')
        self.intp = interpolate.make_interp_spline(x=self.schedule.mission_timeline, y=data.full_paths['clipper'])
        self.launch = self.schedule.launch_time
        p_clipper = self.intp(self.schedule.initial_time)
        p_clipper[2] += 10 # raise launch site above Earth surface
        self.position_actor = VTKUtils.make_point_actor(p_clipper, [1,1,0], 10)
        self.position_actor.GetProperty().RenderPointsAsSpheresOff()

        if clock is not None and clock >= self.launch:
            self.has_launched = True
            points = []
            for t in np.arange(self.launch, clock, Units.hour):
                p_clipper = self.intp(t)
                p_clipper[2] += 10
                points.append(p_clipper)
            self.orbit_actor = VTKUtils.make_curve_actor(points, [0.5, 0.5, 0])
        else:
            self.orbit_actor = vtk.vtkActor()
        self.orbit_actor.GetProperty().RenderLinesAsTubesOff()

    def update(self, et: float):
        p_clipper = self.intp(et)
        p_launch = self.earth.get_position(et)
        # myprint(verbose=verbose, text='Clipper at {}, launch site at {}, launch pos at {}'.format(p_clipper, p_launch, self.launch_pos))
        self.position_actor.SetPosition(p_clipper)

        if not self.has_launched and et >= self.launch:
            self.has_launched = True
            self.orbit_actor = VTKUtils.make_curve_actor([p_clipper, p_clipper], [0.5, 0.5, 0])

        if self.has_launched:
            myprint(verbose=verbose, text='adding one point to curve')
            self.orbit_actor = VTKUtils.add_point_to_curve(self.orbit_actor, p_clipper)
            self.orbit_actor.Modified()

    def get_position(self, et: float):
        return self.intp(et)

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
        spice.furnsh(os.path.join(naif_path, '19F23_VEEGA_L230511_A290930_LP01_V2_scpse.bsp')) # clipper
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
                'Mercury': 88*Units.day,
                'Venus': 225*Units.day,
                'Earth': 365.256*Units.day,
                'Mars': 686.971*Units.day,
                'Jupiter': 4332.59*Units.day,
                'Saturn': 10759.22*Units.day,
                'Uranus': 30688.5*Units.day,
                'Neptune': 60182*Units.day,
                'Io': 4332.59*Units.day,
                'Europa': 4332.59*Units.day,
                'Ganymede': 4332.59*Units.day,
                'Callisto': 4332.59*Units.day,
                'Moon': 365.256*Units.day,
                 }

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
        spice.spkcov(os.path.join(naif_path, '19F23_VEEGA_L230511_A290930_LP01_V2_scpse.bsp'), self.spice_id, etb)
        # arrival time
        arrival_time = spice.str2et('2029 SEP 27 18:26:02.4221 TDB')
        init_time = etb[0]
        final_time = etb[1]
        selected_interval = [init_time, final_time]
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

        # def __init__(self, name, texture_name, meridian, radii, color)
        myprint(verbose=verbose, text='Initializing body objects...')
        self.bodies = {}
        self.tags = {}
        for name in self.all_body_names:
            self.bodies[name] = Body(name, self.texture_files[name], self.meridian_rotation[name], self.radii[name], self.planet_colors[name])
            self.bodies[name].set_states(self.full_paths[name], self.schedule.full_timeline, default_frame)
        
        self.clipper = Clipper(self)

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
    return (et-data.schedule.initial_time)/Units.day

'''
VTK helper functions
'''
class VTKUtils:
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

    def make_curve_actor(pts: List[List[float]], color=[0.5, 0.5, 0], size=3):
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
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetLineWidth(size)
        actor.GetProperty().SetColor(color[0], color[1], color[2])
        return actor

    def add_point_to_curve(curve_actor: vtk.vtkActor, point: List[float]) -> vtk.vtkActor:
        polydata = curve_actor.GetMapper().GetInput()
        points = polydata.GetPoints()
        n = points.GetNumberOfPoints()
        points.InsertNextPoint(point[0], point[1], point[2])
        points.GetData().Modified()
        lines = polydata.GetLines()
        lines.InsertNextCell(2)
        lines.InsertCellPoint(n-1)
        lines.InsertCellPoint(n)
        lines.GetData().Modified()
        polydata.SetPoints(points)
        polydata.SetLines(lines)
        curve_actor.GetMapper().SetInputData(polydata)
        curve_actor.GetMapper().Update()
        return curve_actor

    def copy_transformation(actor_from: vtk.vtkActor, actor_to: vtk.vtkActor):
        actor_to.SetUserMatrix(actor_from.GetUserMatrix())

    def vantage_point(target_name: str, data: DataLoader, et: int) -> Dict[str, List[float]]:
        if target_name == 'Clipper':
            myprint(verbose=verbose, text='handling special case of Clipper')
            p_now = data.full_paths['clipper'].get_position(et)
            p_old = data.full_paths['clipper'].get_position(et-2*Units.hour)
            return { 'pos': p_old, 'focal': p_now, 'up': [0,0,1] }
        elif target_name == 'Clipper-Venus':
            return { 'pos': data.full_paths['clipper'].get_position(et), 'focal': data.bodies['Venus'].get_position(et), 'up': [0, 0, 1] }
        elif target_name == 'Clipper-Earth':
            return { 'pos': data.full_paths['clipper'].get_position(et), 'focal': data.bodies['Earth'].get_position(et), 'up': [0, 0, 1] }
        elif target_name == 'Clipper-Mars':
            return { 'pos': data.full_paths['clipper'].get_position(et), 'focal': data.bodies['Mars'].get_position(et), 'up': [0, 0, 1] }
        elif target_name == 'Clipper-Jupiter':
            return { 'pos': data.full_paths['clipper'].get_position(et), 'focal': data.bodies['Jupiter'].get_position(et), 'up': [0, 0, 1] }
        elif target_name == 'Clipper-Europa':
            return { 'pos': data.full_paths['clipper'].get_position(et), 'focal': data.bodies['Europa'].get_position(et), 'up': [0, 0, 1] }
        elif target_name == 'ecliptic':
            target_pos = data.bodies['Sun'].get_position(et)
            return { 'pos': target_pos + np.array([0, 0, Units.au_to_km(10)]), 'focal': target_pos, 'up': [0, 1, 0] }

        main_body = data.reference_body(target_name)
        main_body_pos = data.bodies[main_body].get_position(et)
        target_pos = data.bodies[target_name].get_position(et)
        if target_name.upper() == 'SUN':
            myprint(verbose=verbose, text='handling special case of Sun')
            dir = target_pos - data.bodies['Earth'].get_position(et)
        else:
            dir = target_pos - main_body_pos
        dir /= np.linalg.norm(dir)
        radius = data.bodies[target_name].get_max_radius()
        return { 'pos': target_pos - 5*radius*dir, 'focal': target_pos, 'up': [0,0,1] }

    def clock_as_str(cl: int) -> str:
        return spice.timout(cl, 'AP:MN:SC AMPM Month DD, YYYY')

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
        return 'p={}, f={}, u={}'.format(pos_str, foc_str, up_str)

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
            actor.SetMapper(mapper)
            actor.SetTexture(ring_texture)
            ring_actors.append(actor)
        return ring_actors

    # Return: actor_plane, actor_perimeter
    # Ecliptic disk and surrounding perimeter
    def make_ecliptic_plane(data: DataLoader) -> List[vtk.vtkActor]:
        ecliptic_actors = []
        ecliptic = vtk.vtkRegularPolygonSource()
        ecliptic.SetCenter(0,0,0)
        ecliptic.SetNumberOfSides(500)
        ecliptic.SetRadius(data.distance('Neptune', 'Sun', data.schedule.initial_time))
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(ecliptic.GetOutputPort())
        actor = vtk.vtkActor()
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
        f_actor.SetMapper(f_mapper)
        f_actor.GetProperty().SetColor(0,0.2,1)
        f_actor.GetProperty().SetOpacity(0.1)
        f_actor.GetProperty().SetAmbient(1)
        f_actor.GetProperty().SetLineWidth(2)
        f_actor.GetProperty().RenderLinesAsTubesOn()
        return actor, f_actor

    # Return: orbit_actors
    # Dictionary associated body name to trajectory depiction actor
    def make_orbits(data: DataLoader) -> Dict[str, vtk.vtkActor]:
        orbit_actors = {}
        sched = data.schedule
        body_groups = [ [ data.planet_names, Units.day ], [ data.jupiter_moon_names, Units.hour ] ]
        for names, period in body_groups:
            for name in names:
                myprint(verbose=verbose, text=' * {}'.format(name))
                body = data.bodies[name]
                points = []
                for t in np.arange(sched.initial_time, sched.initial_time + min(1.01*data.orbital_period[name], sched.final_time-sched.initial_time), period):
                    points.append(body.get_position(t))
                actor = VTKUtils.make_curve_actor(points, data.planet_colors[name], 2)
                orbit_actors[name] = actor
        myprint(verbose=verbose, text=' * Clipper')
        intp = interpolate.make_interp_spline(x=data.schedule.mission_timeline, y=data.full_paths['clipper'])
        points = []
        for t in np.arange(sched.initial_time, sched.final_time, Units.hour):
            points.append(intp(t))
        actor = VTKUtils.make_curve_actor(points, [0.3, 0, 0], 0.25)
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

# Class to encapsulate various aspects of the simulation state
class SimulationState:
    def __init__(self, params: Dict[str, Any], time_step: float, clock: float):
        self.params = params
        self.time_step = time_step
        self.clock = clock
        self.planet_focus_changed = True
        self.saved_relative_position = False 

# Class to encapsulate various graphical elements of the visualization
class GraphicsObjects:
    def __init__(self, window: vtk.vtkRenderWindow, all_actors: Dict[str, List[vtk.vtkActor]], resize: vtk.vtkImageResize):
        self.window = window
        self.all_actors = all_actors
        self.resize = resize
        self.renderers = None
        self.frame_counter = 0

# Overarching class that contains and controls all aspects of the simulation
# and visualization
class Simulation:
    def __init__(self, state: SimulationState, data: DataLoader, graphics: GraphicsObjects):
        self.state = state
        self.data = data
        self.graphics = graphics
        self.render_window = self.graphics.window

    '''
    Series of callback functions to affect the variables of the simulation
    '''
    def key_pressed_cb(self, obj, event):
        print('key pressed')
        print(f'current Saturn position is {self.graphics.all_actors["bodies"]["Saturn"].GetPosition()}')
        print(f'current Earth position is {self.graphics.all_actors["bodies"]["Earth"].GetPosition()}')
        print(f'current ring positions are {self.graphics.all_actors["rings"][0].GetPosition()}')
        new_key = obj.GetKeySym()
        if new_key.isdigit():
            key = int(new_key)
            if key <= 8 and self.data.all_body_names[key] != self.state.params.planet_focus:
                myprint(verbose=True, text='changing focus from {} to {}'.format(self.state.params.planet_focus, self.data.all_body_names[key]))
                self.state.params.planet_focus = self.data.all_body_names[key]
                self.state.planet_focus_changed = True
                self.state.saved_relative_position = False
                self.state.params.do_tether = True
                cam_setting = VTKUtils.vantage_point(target_name=self.state.params.planet_focus, data=self.data, et=self.state.clock)
                VTKUtils.get_camera(self.graphics.window).SetPosition(cam_setting['pos'])
                VTKUtils.get_camera(self.graphics.window).SetFocalPoint(cam_setting['focal'])
                VTKUtils.get_camera(self.graphics.window).SetViewUp(cam_setting['up'])
                self.graphics.all_actors['text']['bodies'].SetInput('Current Focus: {}'.format(self.state.params.planet_focus))

            # Overrides the default vtk behavior for keypress '3'
            if int(new_key) == 3:
                self.graphics.window.StereoRenderOn()
        elif new_key == 'minus':
            self.state.params.planet_focus = 'N/A'
            self.state.planet_focus_changed = False
            self.graphics.all_actors['text']['bodies'].SetInput('Current Focus: None')
            self.state.params.do_tether = False
        elif new_key == 'o' or new_key == 'O':
            self.state.params.show_orbits = not self.state.params.show_orbits
            if self.state.params.show_orbits:
                myprint(verbose=True, text='now showing orbits')
                for name in self.graphics.all_actors['orbits']:
                    if self.state.clock < self.data.schedule.arrival_time and name in self.data.jupiter_moon_names:
                        continue
                    myprint(verbose=True, text='added orbit of {}'.format(name))
                    myprint(verbose=True, text='there are {} vertices in this orbit'.format(self.graphics.all_actors['orbits'][name].GetMapper().GetInput().GetPoints().GetNumberOfPoints()))
                    VTKUtils.get_renderer(self.graphics.window).AddActor(self.graphics.all_actors['orbits'][name])
            else:
                myprint(verbose=True, text='now hiding orbits')
                for name in self.graphics.all_actors['orbits']:
                    VTKUtils.get_renderer(self.graphics.window).Remove_actor(self.graphics.all_actors['orbits'][name])
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
                new_focus_planet = 'Clipper Backward'
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
                self.state.params.planet_focus = self.state.new_focus_planet.split()[0]
                self.state.planet_focus_changed = True
                self.state.saved_relative_position = False
                self.state.params.do_tether = True
                camera = VTKUtils.get_camera(self.graphics.window)
                cam_setting = VTKUtils.vantage_point(target_name=self.state.params.planet_focus, data=self.data, et=self.state.clock)
                camera.SetPosition(cam_setting['pos'])
                camera.SetFocalPoint(cam_setting['focal'])
                camera.SetViewUp(cam_setting['up'])
                self.graphics.all_actors['text']['title'].SetInput('Current Focus: {}'.format(self.state.params.planet_focus))
            camera = VTKUtils.get_camera(self.graphics.window)
            if self.state.params.zoom_factor != 1:
                camera.Zoom(self.state.params.zoom_factor)
                self.state.params.zoom_factor = 1
                myprint(verbose=True, text='camera view angle: {}'.format(camera.GetViewAngle()))
        elif new_key == 't':
            # tether mode
            if not self.state.params.do_tether:
                self.state.params.do_tether = True
                myprint(verbose=True, text='tethering activated')
            else:
                self.state.params.do_tether = False
                self.state.saved_relative_position = None
                myprint(verbose=True, text='tethering deactivated')
        elif new_key == 'h':
            myprint(verbose=True, text='List of keyboard commands:')
            myprint(verbose=True, text='\'0\', \'1\', \'2\', ..., \'8\': Point the camera to the Sun, Mercury, Venus, ..., Neptune and track')
            myprint(verbose=True, text='\'!\', \'@\', \'#\', \'$\':  Point the camera to Jupiter\'s Galilean moons')
            myprint(verbose=True, text='\'^\':                     Point the camera to the Moon')
            myprint(verbose=True, text='\'&\':                     Point the camera to Clipper\'s position')
            myprint(verbose=True, text='\'-\':                     Turn off tracking')
            myprint(verbose=True, text='\'f\':                     Toggle full screen on and off')
            myprint(verbose=True, text='\'a\', ..., \'d\':           Point camera from Clipper to Earth, Mars, Jupiter, Europa')
            myprint(verbose=True, text='\'Ctrl-9\':                Point camera from Clipper to the Moon')
            myprint(verbose=True, text='\'Ctrl-!\', \'Ctrl-@\':      Point camera from Clipper to Jupiter\'s moon')
            myprint(verbose=True, text='\'t\':                     Toggle tethering on and off')
            myprint(verbose=True, text='\'v\':                     Point the camera to the Sun\'s North pole')
            myprint(verbose=True, text='\'o\' (\'O\'):               Turn on/off depiction of planets\' orbits')
            myprint(verbose=True, text='\'h\':                     Print this information')
        else:
            myprint(verbose=verbose, text='unrecognized entered key is {}'.format(new_key))

    def time_step_change_cb(self, value):
        print('time step change')
        self.state.time_step = float(value*Units.hour)/self.state.params.refresh_rate
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
            actor.SetScale(self.state.params.planet_scale, self.state.params.planet_scale, self.state.params.planet_scale)
        self.graphics.window.Render()

    def window_size_change_cb(self, obj, event):
        print('window size change')
        width, height = self.graphics.window.GetSize()
        if self.state.params.resolution[0] != width or self.state.params.resolution[1] != height:
            self.state.params.resolution[0] = width
            self.state.params.resolution[1] = height
            if not self.state.params.shadow:
                self.graphics.resize.SetOutputDimensions(self.state.params.resolution[0], self.state.params.resolution[1], 0)
                self.graphics.resize.Update()
                myprint(verbose=verbose, text='texture has been resized')
            self.graphics.window.Render()

    def camera_change_cb(self, caller, event):
        print('camera change')
        if not self.state.params.in_progress:
            self.state.params.cam_changed = True

    def interaction_start_cb(self, caller, event):
        print('interaction starts')
        if self.state.params.do_tether:
            self.state.params.paused = True

    def interaction_end_cb(self, caller, event):
        print('interaction ends')
        if self.state.params.do_tether:
            self.state.params.paused = False

    def date_change_cb(self, value):
        current = self.state.clock
        et = value
        self.state.clock = et*Units.hour + self.data.schedule.initial_time

        renderer = VTKUtils.get_renderer(self.render_window)
        renderer.RemoveActor(self.graphics.all_actors['clipper orbit'])
        renderer.RemoveActor(self.graphics.all_actors['clipper position'])
        # renderer.RemoveActor(data.full_paths['clipper'].launch_orbit_actor)
        # renderer.RemoveActor(data.full_paths['clipper'].launch_site_curve_actor)
        self.data.clipper = Clipper(self.data, self.state.clock)
        renderer.AddActor(self.data.clipper.orbit_actor)
        renderer.AddActor(self.data.clipper.position_actor)
        # renderer.AddActor(self.data.full_paths['clipper'].launch_orbit_actor)
        # renderer.AddActor(data.full_paths['clipper'].launch_site_curve_actor)

    # function called at every step of the internal clock of the simulation
    def timer_cb(self, obj, value):
        camera = self.graphics.renderers['bodies'].GetActiveCamera()
        print('timer callback: clock={}'.format(self.state.clock))
        cam = self.render_window.GetRenderers().GetFirstRenderer().GetActiveCamera()
        if self.state.params.paused:
            return
        clock = self.state.clock
        timestr = VTKUtils.clock_as_str(clock)
        # print(f'time string is {timestr}')
        if clock >= self.data.schedule.final_time:
            return
        current_focal_vector = VTKUtils.focal_vector(camera)
        # print(f'current focal vector is {current_focal_vector}')
        # Update camera for currently selected planet
        # print(f'focus planet = {self.state.params.focus_planet}')
        if self.state.params.focus_planet != 'N/A':
            # tether the camera position to the moving planet to maintain the relative position
            # that it has at the beginning or after a user change to the camera
            self.state.params.show_clipper_orbit = True
            self.state.params.in_progress = True
            if self.state.params.focus_planet.split()[0] == 'Clipper':
                new_body_pos = self.data.full_paths['clipper'].get_position(clock)
            elif self.state.params.focus_planet == 'Clipper-Earth':
                new_body_pos = self.data.bodies['Earth'].get_position(clock)
            elif self.state.params.focus_planet == 'Clipper-Mars':
                new_body_pos = self.data.bodies['Mars'].get_position(clock)
            elif self.state.params.focus_planet == 'Clipper-Jupiter':
                new_body_pos = self.data.bodies['Jupiter'].get_position(clock)
            elif self.state.params.focus_planet == 'Clipper-Europa':
                new_body_pos = self.data.bodies['Europa'].get_position(clock)
            elif self.state.params.focus_planet == 'Clipper-Venus':
                new_body_pos = self.data.bodies['Venus'].get_position(clock)
            elif self.state.params.focus_planet != 'ecliptic':
                new_body_pos = self.data.bodies[self.state.params.focus_planet].get_position(
                    clock)
            else:
                new_body_pos = self.data.bodies['Sun'].get_position(clock)
            camera.SetFocalPoint(new_body_pos)
            self.state.params.in_progress = False
            if self.state.params.focus_planet == 'Clipper-Earth' or \
                self.state.params.focus_planet == 'Clipper-Mars' or \
                self.state.params.focus_planet == 'Clipper-Jupiter' or \
                self.state.params.focus_planet == 'Clipper-Europa' or \
                self.state.params.focus_planet == 'Clipper-Venus':
                self.state.params.show_data.full_paths['clipper'] = False
                planet_name = self.state.params.focus_planet[8:]
                clipper_new_pos = self.data.full_paths['clipper'].get_position(clock)
                self.state.params.do_tether = False
                self.state.params.in_progress = True
                self.graphics.camera.SetPosition(clipper_new_pos)
                self.state.params.in_progress = False
                self.text.SetInput('Current Focus: {}'.format(self.state.params.focus_planet))
            else:
                cam_pos = np.array(camera.GetPosition())
                self.graphics.all_actors['text']['bodies'].SetInput('Current Focus: ' + self.state.params.focus_planet)
                if self.state.params.do_tether:
                    if self.state.params.saved_relative_position is None or self.state.params.cam_changed:
                        self.state.params.saved_relative_position = current_focal_vector
                        self.state.params.cam_changed = False
                    self.state.params.in_progress = True
                    camera.SetPosition(new_body_pos - self.state.params.saved_relative_position)
                    self.state.params.in_progress = False

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

        # extend Clipper's orbit
        self.data.clipper.update(clock)

        self.graphics.all_actors['text']['time'].SetInput(timestr)
        obj.GetRenderWindow().GetRenderers().GetFirstRenderer().RemoveActor(self.data.clipper.orbit_actor)
        if self.state.params.show_clipper_orbit:
            obj.GetRenderWindow().GetRenderers().GetFirstRenderer().AddActor(self.data.clipper.orbit_actor)
        obj.GetRenderWindow().Render()

        if not self.state.params.paused:
            # print('state is not paused: incrementing')
            self.state.clock += self.state.params.time_step
        else:
            print('state is paused')

    def screenshot_cb(self):
        # ---------------------------------------------------------------
        # Save current contents of render window to PNG file
        # ---------------------------------------------------------------
        file_name = 'screenshot' + str(self.graphics.frame_counter).zfill(5) + ".png"
        image = vtk.vtkWindowToImageFilter()
        image.SetInput(self.graphics.window)
        png_writer = vtk.vtkPNGWriter()
        png_writer.SetInputConnection(image.GetOutputPort())
        png_writer.SetFileName(file_name)
        self.graphics.window.Render()
        png_writer.Write()
        self.graphics.frame_counter += 1

class MainWindow(QMainWindow):
    def __init__(self, parent = None):
        QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.mission_state = None
        self.mission_data = None
        self.clipper_mission = None
        print('Main window initialized')

    def timer_callback(self):
        print('in timer callback')
        self.ui.slider_date.setValue(int((self.mission_state.clock-self.mission_data.schedule.initial_time)/Units.hour))
        self.clipper_mission.timer_cb(obj=None, event=None)

    def main(self, args):
        verbose = args.verbose
        ''' Mission Data '''
        print('creating mission data')
        self.mission_data = DataLoader(media_path=args.media_path, naif_path=args.naif_path, scale=args.scale)
        print('mission data created')

        # Timer Frequency
        timer_period = 1000/args.refresh_rate  # in ms
        args.time_step = Units.second/args.refresh_rate
        args.zoom_factor = 1

        # global state variables
        args.planet_scale = 1
        args.ambient = 0.1
        args.paused = False

        # interaction state variables
        args.planet_focus_changed = False
        args.do_tether = False
        args.saved_relative_position = None
        args.in_progress = False
        args.cam_changed = True

        # Render size set to HD Display
        # args.planet_focus = 'Sun'  # Initial focus planet is the sun
        args.show_clipper_orbit = True

        ''' Mission State '''
        print('creating mission state')
        self.mission_state = SimulationState(params=args, time_step=args.time_step, clock=self.mission_data.schedule.initial_time)
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
        all_actors['rings'] = VTKUtils.make_rings(self.mission_data)
        print('after initialization, ring actor has position {}'.format(all_actors['rings'][0].GetPosition()))
        myprint(verbose=verbose, text='...done.')

        myprint(verbose=verbose, text='Creating ecliptic plane...')
        all_actors['ecliptic'] = VTKUtils.make_ecliptic_plane(self.mission_data)
        myprint(verbose=verbose, text='...done.')

        myprint(verbose=verbose, text='Creating bodies\' orbit curves...')
        all_actors['orbits'] = VTKUtils.make_orbits(self.mission_data)
        myprint(verbose=verbose, text='...done.')

        myprint(verbose=verbose, text='Creating Clipper object...')
        all_actors['clipper position'] = self.mission_data.clipper.position_actor
        all_actors['clipper orbit'] = self.mission_data.clipper.orbit_actor
        myprint(verbose=verbose, text='...done.')

        idx = 0
        myprint(verbose=verbose, text='Initializing VTK actors...')
        body_actors = {}
        body_point_actors = {}
        for name in self.mission_data.all_body_names:
            print(f'processing {name}')
            body = self.mission_data.bodies[name]

            # Create Actors and add textures, rotations, and orbit position
            body_actors[name] = body.move_actor(self.mission_state.clock)
            body_actors[name].GetProperty().SetAmbient(self.mission_state.params.ambient)
            body_actors[name].GetProperty().SetSpecular(0) #.01)
            body_actors[name].GetProperty().SetDiffuse(0.65)

            body_point_actors[name] = body.move_point_actor(self.mission_state.clock)

            if name == 'Saturn':
                myprint(verbose=verbose, text='moving rings')
                for actor in all_actors['rings']:
                    # print(f'ring actor is {actor}')
                    # print(f'Saturn\'s position is {body.get_position(self.mission_state.clock)}')
                    # print(f'Saturn\'s actor position is {body.get_actor().GetPosition()}')
                    # actor.SetPosition(body.get_position(self.mission_state.clock))
                    # print(f'ring actor moved to {actor.GetPosition()}')
                    VTKUtils.copy_transformation(body.get_actor(), actor)
                    # print(f'After transformation, ring actor is {actor}')
                    actor.GetProperty().SetAmbient(self.mission_state.params.ambient)
        myprint(verbose=verbose, text='...done.')
        all_actors['bodies'] = body_actors
        all_actors['points'] = body_point_actors
        resize_image = vtk.vtkImageResize()

        ''' Mission Graphics '''
        mission_graphics = GraphicsObjects(window=render_window, all_actors=all_actors, resize=resize_image)

        ''' Mission '''
        self.clipper_mission = Simulation(state=self.mission_state, data=self.mission_data, graphics=mission_graphics)

        # Create Lighting
        sunlight1 = vtk.vtkLight()
        sunlight1.SetPosition(0, 0, 0)
        sunlight1.SetFocalPoint(0.0000001*self.mission_state.params.scale, 0, 0)  # Right hemisphere lighting
        sunlight1.SetConeAngle(180)
        sunlight1.SetPositional(0)

        # Left hemisphere lighting
        sunlight2 = vtk.vtkLight()
        sunlight2.SetPosition(0, 0, 0)
        sunlight2.SetFocalPoint(-0.0000001*self.mission_state.params.scale, 0, 0)
        sunlight2.SetConeAngle(180)
        sunlight2.SetPositional(0)

        # Saturn specific lighting
        sunlight3 = vtk.vtkLight()
        sunlight3.SetPosition(0, 0, 0)
        sunlight3.SetFocalPoint(self.mission_data.bodies['Saturn'].get_position(et=self.mission_state.clock))
        sunlight3.SetConeAngle(90)
        sunlight3.SetPositional(1)

        # Add Background Starry Image
        if not args.show_shadows:
            myprint(verbose=verbose, text='Loading background image...')
            jpeg_reader = vtk.vtkJPEGReader()
            jpeg_reader.SetFileName(self.mission_data.texture_files['Stars'])
            jpeg_reader.Update()

            # Resize Image to Match Render Window Size
            #resize_image = vtk.vtkImageResize()
            resize_image.SetInputData(jpeg_reader.GetOutput())
            resize_image.SetResizeMethod(0)
            if args.resolution[0] > 0 and args.resolution[1] > 0:
                resize_image.SetOutputDimensions(args.resolution[0], args.resolution[1], 0)
            resize_image.Update()

            background_mapper = vtk.vtkImageMapper()
            background_mapper.SetInputData(resize_image.GetOutput())
            background_actor = vtk.vtkActor2D()
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
        # body_renderer.UseShadowsOn()
        renderers['bodies'].SetBackground(0, 0, 0)  # black background
        renderers['bodies'].AddLight(sunlight1)
        renderers['bodies'].AddLight(sunlight2)
        # renderers['bodies'].AddLight(sunlight3)
        myprint(verbose=verbose, text='Adding all actors to renderer...')
        for name in self.mission_data.all_body_names:
            renderers['bodies'].AddActor(body_actors[name])
            renderers['bodies'].AddActor(body_point_actors[name])
        for actor in all_actors['rings']:
            # myprint(verbose=verbose, text='adding one ring actor {}'.format(actor))
            renderers['bodies'].AddActor(actor)
        if args.show_ecliptic:
            for actor in self.graphics.ecliptic_actors:
                renderers['bodies'].AddActor(actor)
        if args.show_orbits:
            myprint(verbose=verbose, text='orbit_actors={}'.format(mission_graphics.all_actors['orbits']))
            for body, actor in all_actors['bodies'].items():
                myprint(verbose=verbose, text='body={}, actor={}'.format(body, actor))
                renderers['bodies'].AddActor(actor)
        if args.show_clipper_orbit:
            renderers['bodies'].AddActor(all_actors['clipper orbit'])
        renderers['bodies'].AddActor(all_actors['clipper position'])
        myprint(verbose=verbose, text='...done.')

        # Create the camera focused on the Sun
        cam1 = vtk.vtkCamera()
        cam1.SetFocalPoint(0, 0, 0)  # Sun Center
        cam1.SetPosition(-self.mission_data.distance('Mercury', 'Sun', self.mission_state.clock), 0, 0.5*self.mission_data.distance('Mercury', 'Sun', self.mission_state.clock))
        cam1.Elevation(0)  # look at off-angle system (0 is top down)
        cam1.Azimuth(0)
        cam1.SetViewUp(0,0,1)
        cam1.SetClippingRange(100, 100000000000)
        myprint(verbose=verbose, text='camera view angle is: {}'.format(cam1.GetViewAngle()))

        renderers['bodies'].SetActiveCamera(cam1)
        if not args.show_shadows:
            renderers['bodies'].SetLayer(1)  # In front of background image

        if args.show_shadows:
            myprint(verbose=verbose, text='adding shadows')
            render_window.SetMultiSamples(0)
            shadows = vtk.vtkShadowMapPass()
            # shadows.DebugOn()
            sequence = vtk.vtkSequencePass()
            # sequence.DebugOn()
            passes = vtk.vtkRenderPassCollection()
            # passes.DebugOn()
            passes.AddItem(shadows.GetShadowMapBakerPass())
            passes.AddItem(shadows)
            sequence.SetPasses(passes)
            camera_pass = vtk.vtkCameraPass()
            # camera_pass.DebugOn()
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
        text_actors['time'].SetInput('{}'.format(spice.timout(self.mission_data.schedule.initial_time, 'AP:MN:SC AMPM Month DD, YYYY')))
        text_actors['time'].GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['time'].GetPositionCoordinate().SetValue(0.01, 0.95)
        text_actors['time'].GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actors['time'].GetPosition2Coordinate().SetValue(0.25, 0.95)
        text_actors['time'].GetTextProperty().SetFontSize(10)
        text_actors['time'].GetTextProperty().SetColor(1.0, 1.0, 1.0)

        all_actors['text'] = text_actors
        myprint(verbose=verbose, text='...done')

        renderers['bodies'].AddActor2D(text_actors['time'])
        renderers['bodies'].AddActor2D(text_actors['bodies'])

        render_window.AddRenderer(renderers['bodies'])

        # Create a background renderer
        if not args.show_shadows:
            background_renderer = vtk.vtkRenderer()
            # background_renderer.DebugOn()
            background_renderer.SetLayer(0)
            background_renderer.InteractiveOff()
            background_renderer.AddActor(background_actor)
            render_window.AddRenderer(background_renderer)
            # render_window.DebugOn()

        # visualize ecliptic plane
        if args.show_ecliptic:
            for actor in all_actors['ecliptic']:
                renderers['bodies'].AddActor(actor)

        self.clipper_mission.graphics.renderers = renderers

        # Set-up interactor
        render_window_interactor = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        # render_window_interactor.DebugOn()
        style = vtk.vtkInteractorStyleTrackballCamera()
        style.AutoAdjustCameraClippingRangeOff()
        # style.DebugOn()
        render_window_interactor.SetInteractorStyle(style)
        render_window_interactor.Initialize()
        render_window.Render()

        render_window.AddObserver('WindoResizeEvent', self.clipper_mission.window_size_change_cb)

        myprint(verbose=verbose, text='Creating sliders...')

        # Setting up widgets
        self.ui.timescale_spinbox.setValue(self.mission_state.time_step)
        self.ui.planetscale_dial.setRange(1, int(args.max_planet_scale))
        self.ui.planetscale_dial.setValue(self.mission_state.params.planet_scale)
        self.ui.date_slider.setRange(0, int(what_day(data=self.mission_data, et=self.mission_data.schedule.final_time)))
        self.ui.date_slider.setOrientation(QtCore.Qt.Horizontal)
        self.ui.date_slider.setValue(int(what_day(self.mission_data, self.mission_state.clock)))
        self.ui.date_slider.setTracking(False)
        self.ui.date_slider.setTickInterval(365)
        self.ui.date_slider.setTickPosition(QSlider.TicksAbove)

        myprint(verbose=verbose, text='Initializing interactor and adding observers')
        interactor = render_window.GetInteractor()
        # interactor.AddObserver('TimerEvent', self.timer_callback)
        interactor.AddObserver('KeyPressEvent', self.clipper_mission.key_pressed_cb)
        # interactor.AddObserver('InteractionEvent', interaction_callback)
        interactor.AddObserver('StartInteractionEvent', self.clipper_mission.interaction_start_cb)
        interactor.AddObserver('EndInteractionEvent', self.clipper_mission.interaction_end_cb)
        # interactor.CreateRepeatingTimer(int(timer_period))  # ms between calls
        self.timer = QTimer(self)
        self.timer.setInterval(int(timer_period))
        print(f'interval set to {timer_period}')

        # interactor.EnableRenderOn()
        myprint(verbose=verbose, text='...done')
        render_window.Render()

        myprint(verbose=verbose, text='Connecting widget to callbacks')
        self.ui.timescale_spinbox.valueChanged.connect(self.clipper_mission.time_step_change_cb)
        self.ui.date_slider.valueChanged.connect(self.clipper_mission.date_change_cb)
        self.ui.planetscale_dial.valueChanged.connect(self.clipper_mission.planet_scale_change_cb)
        self.ui.push_screenshot.clicked.connect(self.clipper_mission.screenshot_cb)
        self.timer.timeout.connect(self.timer_callback)
        myprint(verbose=verbose, text='...done')
        sys.exit(app.exec_())

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Visualize Clipper Spacecraft and Solar System using NAIF Data',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--show_ecliptic', action='store_true', help='Show ecliptic plane')
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
    args = parser.parse_args()

    app = QApplication(sys.argv)
    window = MainWindow()
    window.ui.vtkWidget.GetRenderWindow().SetSize(args.resolution[0], args.resolution[1])
    window.setWindowState(Qt.WindowMaximized)
    window.show()
    window.main(args)
