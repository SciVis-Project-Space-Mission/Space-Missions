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

# faulthandler.enable()

# last modification, interrupted, attempted to add the ability to control
# the time range of the visualization, its time step, and the automatic
# recording of each individual frame, as well as the specification of the
# camera setting, including in terms of the bodies and Clipper's position

# load the data
print('loading the data from NAIF files...')
spice.furnsh('naif/19F23_VEEGA_L230511_A290930_LP01_V2_scpse.bsp') # clipper
spice.furnsh('naif/naif0012.tls') # bodies' dynamics
spice.furnsh('naif/pck00010.tpc') # bodies' constant values and orientation
print('...done')

# set up constants
day     = 86400
hour    = 3600
minute  = 60

# handy unit conversions
def au_to_km(au):
    return au * 149597871

def km_to_au(km):
    return km / 149597871

def deg_to_rad(deg):
    return deg/180.*math.pi

def to_unit_rgb(rgb):
    r = float(rgb[0]/255.)
    g = float(rgb[1]/255.)
    b = float(rgb[2]/255.)
    return [r, g, b]

# set clipper SPICE ID
scid    = -159
scidstr = '-159'

planet_names_spice = [ 'Mercury', 'Venus', 'Earth', 'Mars Barycenter',
                       'Jupiter', 'Saturn barycenter', 'Uranus barycenter',
                       'Neptune barycenter' ]
planet_names = [ name.split()[0] for name in planet_names_spice ]
jupiter_moon_names = [ 'Io', 'Europa', 'Ganymede', 'Callisto' ]
earth_moon_name = [ 'Moon' ]
all_body_names = ['Sun'] + planet_names + jupiter_moon_names + earth_moon_name
all_body_names_spice = ['Sun'] + planet_names_spice + jupiter_moon_names + earth_moon_name

planet_scale = 1
ambient = 0.1

saturn_rings = { 'D': (66900, 74510), 'C': (74658, 92000),
                 'B': (92000, 117580), 'Cassini': (117580, 122170),
                 'A': (122170, 136775), 'Roche': (136775, 139380),
                 'F': (140183, 140680) }

texture_files = { 'Sun':     'media/2k_sun_bright.jpg',
                  'Mercury': 'media/2k_mercury.jpg',
                  'Venus':   'media/2k_venus.jpg',
                  # 'Earth':   'media/2k_earth_daymap.jpg',
                  'Earth':   'media/earth/world.topo.200412.3x5400x2700.jpg',
                  'Mars':    'media/2k_mars.jpg',
                  'Jupiter': 'media/2k_jupiter.jpg',
                  'Saturn':  'media/2k_saturn.jpg',
                  'Uranus':  'media/2k_uranus.jpg',
                  'Neptune': 'media/2k_neptune.jpg',
                  'Stars':   'media/starry_background.jpg',
                  'Rings':   'media/ring_colors2d.jpg',
                  'Io':      'media/2k_io.jpg',
                  'Europa':  'media/2k_europa.jpg',
                  'Ganymede':'media/2k_ganymede.jpg',
                  'Callisto':'media/2k_callisto.jpg',
                  'Moon':    'media/2k_moon.jpg'
                }

meridian_rotation = { 'Sun': 0,
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

orbital_period = { 'Sun': 0,
                   'Mercury': 88*day,
                   'Venus': 225*day,
                   'Earth': 365.256*day,
                   'Mars': 686.971*day,
                   'Jupiter': 4332.59*day,
                   'Saturn': 10759.22*day,
                   'Uranus': 30688.5*day,
                   'Neptune': 60182*day,
                   'Io': 4332.59*day,
                   'Europa': 4332.59*day,
                   'Ganymede': 4332.59*day,
                   'Callisto': 4332.59*day,
                   'Moon': 365.256*day,
                 }

planet_colors = { 'Sun':      to_unit_rgb([255,167,0]),
                  'Mercury':  to_unit_rgb([112,109,113]),
                  'Venus':    to_unit_rgb([240,218,166]),
                  'Earth':    to_unit_rgb([44,49,115]),
                  'Mars':     to_unit_rgb([244,108,76]),
                  'Jupiter':  to_unit_rgb([240,210,145]),
                  'Saturn':   to_unit_rgb([200,170,100]),
                  'Uranus':   to_unit_rgb([193,238,238]),
                  'Neptune':  to_unit_rgb([81,106,161]),
                  'Io':       to_unit_rgb([252,240,113]),
                  'Europa':   to_unit_rgb([148,115,68]),
                  'Ganymede': to_unit_rgb([81,70,59]),
                  'Callisto': to_unit_rgb([86,77,53]),
                  'Moon':     to_unit_rgb([128, 128, 128])
                 }

# set frame and corrections
default_frame = 'ECLIPJ2000'
abcorr = 'NONE'

print('Setting temporal domain...')
# coverage dates for Clipper
etb = spice.cell_double(10000)
# entire duration of the mission
spice.spkcov('naif/19F23_VEEGA_L230511_A290930_LP01_V2_scpse.bsp', scid, etb)
# arrival time
etarrive= spice.str2et('2029 SEP 27 18:26:02.4221 TDB')

init_time = etb[0]
final_time = etb[1]
clock = init_time
selected_interval = [init_time, final_time]

paused = False
video_mode = False

# set up time array
# total duration of the covered portion of Clipper in days
etfull      = np.arange(etb[0], etb[1], day)
# 1.5 day from launch in minutes
etlaunch    = etb[0] + np.arange(0, day*1.5, minute)
# interplanetary travel in days
etinterpl   = np.arange(etb[0], etarrive, day)
# moon tour in hours
ettour      = np.arange(etarrive, etb[1], hour)

print(f'etfull={etfull}, etarrive={etarrive}, ettour={ettour}')
print('time to launch: {} minutes, interplanetary travel time: {} days, touring time {} hours'.format(1.5*day/minute, (etarrive-init_time)/day-1.5, (etb[1]-etarrive)/hour))

etinterpl_no_overlap = np.arange(etb[0]+day*1.5, etarrive, day)
etclipper = np.concatenate((etlaunch, etinterpl_no_overlap, ettour))
# print('etclipper={}'.format(etclipper))

# interaction state variables
planet_focus_changed = False
do_tether = False
saved_relative_position = None
in_progress = False
cam_changed = True

# query states during entire mission for clipper and relevant planets
# states are given in (km, km/s)
full_paths = {}
interpl_paths = {}
tour_paths = {}

print('Acquiring bodies\' paths across time domain...')
for id, name in enumerate(all_body_names):
    print(' * {}'.format(name))
    spice_name = all_body_names_spice[id]
    full_paths[name], not_used = spice.spkezr(spice_name, etfull, default_frame, abcorr, 'Sun')
    interpl_paths[name], not_used = spice.spkezr(spice_name, etinterpl, default_frame, abcorr, 'Sun')
    full_paths[name] = np.array(full_paths[name])[:,:3]
    interpl_paths[name] = np.array(interpl_paths[name])[:,:3]

for name in jupiter_moon_names:
    print(' * {}'.format(name))
    tour_paths[name], not_used = spice.spkezr(name, ettour, default_frame, abcorr, 'Jupiter')
    tour_paths[name] = np.array(tour_paths[name])[:,:3]

# clipper trajectory
print(' * Clipper')
clipper_full, not_used = spice.spkezr(scidstr, etclipper, default_frame, abcorr, 'Sun')
clipper_full = np.array(clipper_full)[:,:3]
print('...done')

# NAIF/SPICE helper functions
def get_radii(body_name):
    radii = spice.bodvrd( body_name.split()[0], 'RADII', 3 )
    print('radii={}'.format(radii[1]))
    return radii[1]

def coordinate_transform(coord, from_frame, to_frame, et):
    mat = spice.pxform(from_frame, to_frame, et)
    return spice.mxv(mat, coord)

def body_orientation(body_name, et, to_frame=default_frame, from_frame=None):
    if from_frame is None:
        from_frame = 'IAU_' + body_name
    return spice.pxform(from_frame, to_frame, et)

# VTK helper functions
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

def make_curve_actor(pts, color=[0.5, 0.5, 0], size=3):
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

def add_point_to_curve(curve_actor, point):
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

class Body(object):
    def __init__(self, name):
        self.name = name.split()[0]
        self.img_filename = texture_files[self.name]
        sphere = vtk.vtkTexturedSphereSource()
        sphere.SetRadius(1)
        sphere.SetThetaResolution(40)
        sphere.SetPhiResolution(40)

        shape = vtk.vtkTransform()
        shape.RotateZ(-meridian_rotation[self.name])
        self.radii = get_radii(self.name)
        shape.Scale(self.radii[0], self.radii[1], self.radii[2])
        transform = vtk.vtkTransformPolyDataFilter()
        transform.SetInputConnection(sphere.GetOutputPort())
        transform.SetTransform(shape)
        transform.Update()
        self.no_point = True
        self.ellipsoid = transform.GetOutput()
        if name != 'Moon' and name != 'Callisto' and name != 'Io' \
            and name != 'Europa' and name != 'Ganymede':
            self.no_point = False
        if not self.no_point:
            self.point_actor = make_point_actor([0,0,0], [1,1,1], 4)
            self.point_actor.GetProperty().RenderPointsAsSpheresOn()
        else:
            self.point_actor = vtk.vtkActor()
        if name == 'Sun':
            self.point_actor.GetProperty().SetAmbient(1000)
            self.point_actor.GetProperty().SetAmbientColor(planet_colors['Sun'])

        self.local_frame = 'IAU_' + self.name.upper()
        self.states = None
        self.times = None
        self.actor = None
        self.orbit = None
        self.orbit_frame = None

    def GetRadii(self):
        return self.radii

    def GetMaxRadius(self):
        return np.amax(self.radii)

    def SetStates(self, states, ets, frame):
        self.states = states
        self.times = ets
        self.orbit_frame = frame
        self.orbit = interpolate.make_interp_spline(self.times, self.states)

    def GetCover(self):
        if self.times is not None:
            return [self.times[0], self.times[-1]]
        else:
            print('Invalid request: GetCover: no time coordinates available')

    def GetLocalFrameName(self):
        return self.local_frame

    def GetGeometry(self):
        return self.ellipsoid

    def GetActor(self):
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

    def ScaleActor(self, scaling):
        self.actor.SetScale(scaling, scaling, scaling)

    def GetPosition(self, et=clock, frame=default_frame):
        if self.orbit is not None:
            if self.orbit_frame == frame:
                return self.orbit(et)
            else:
                print('doing coordinate transformation')
                return coordinate_transform(self.orbit(et), self.orbit_frame, frame, et)
        else:
            print('Invalid request: GetPosition: no states available')

    def GetRotation(self, et=clock, frame=default_frame):
        mat = body_orientation(self.name, et, frame)
        return mat

    def GetMatrix(self, et=clock, frame=default_frame):
        mat = self.GetRotation(et, frame)
        pos = self.GetPosition(et, frame)
        matrix = vtk.vtkMatrix4x4()
        matrix.Identity()
        for row in range(3):
            matrix.SetElement(row, 0, mat[row, 0])
            matrix.SetElement(row, 1, mat[row, 1])
            matrix.SetElement(row, 2, mat[row, 2])
            matrix.SetElement(row, 3, pos[row])
        return matrix;

    def GetLocalPosition(self, pos, et=clock):
        mat = self.GetRotation(et)
        center =self.GetPosition(et)
        dp = pos-center
        d = np.linalg.norm(dp)
        e1 = mat[:,0]
        e2 = mat[:,1]
        e3 = mat[:,2]
        x = np.dot(dp, e1)
        y = np.dot(dp, e2)
        z = np.dot(dp, e3)
        return [x, y, z]

    def GetGlobalPosition(self, pos, et=clock):
        mat = self.GetRotation(et)
        center = self.GetPosition(et)
        return center + pos[0]*mat[:,0] + pos[1]*mat[:,1] + pos[2]*mat[:,2]

    def GetEquatorialPlane(self, et=clock, frame=default_frame):
        mat = self.GetRotation(et, frame)
        p = np.array(self.GetPosition(et, frame))
        e1 = mat[:,0]
        e2 = mat[:,1]
        return [p, e1, e2]

    def MoveActor(self, et=clock, frame=default_frame):
        actor = self.GetActor()
        actor.SetUserMatrix(self.GetMatrix(et, frame))
        return actor

    def MovePointActor(self, et=clock, frame=default_frame):
        actor = self.point_actor
        actor.SetUserMatrix(self.GetMatrix(et, frame))
        return actor

print('Initializing body objects...')
bodies = {}
tags = {}
for name in all_body_names:
    bodies[name] = Body(name)
    bodies[name].SetStates(full_paths[name], etfull, default_frame)
    bodies[name].point_actor.GetProperty().SetColor(planet_colors[name])

class ClipperOrbit:
    global etclipper
    global clipper_full
    global init_time

    def launch_time(self):
        for t in np.arange(init_time, final_time, hour):
            p_clipper = self.GetPosition(t)
            p_earth = bodies['Earth'].GetPosition(t)
            dist = np.linalg.norm(p_earth - p_clipper)
            if dist >= 6375:
                return t
        return final_time

    def __init__(self):
        self.has_launched = False
        self.intp = interpolate.make_interp_spline(x=etclipper, y=clipper_full)
        self.launch = self.launch_time()
        p_clipper = self.intp(clock)
        p_clipper[2] += 10
        self.position_actor = make_point_actor(p_clipper, [1,1,0], 10)
        self.position_actor.GetProperty().RenderPointsAsSpheresOff()

        if clock >= self.launch:
            self.has_launched = True
            points = []
            for t in np.arange(self.launch, clock, hour):
                p_clipper = self.GetPosition(t)
                p_clipper[2] += 10
                points.append(p_clipper)
            self.orbit_actor = make_curve_actor(points, [0.5, 0.5, 0])
        else:
            self.orbit_actor = vtk.vtkActor()
        self.orbit_actor.GetProperty().RenderLinesAsTubesOff()

    def Update(self, et=clock):
        p_clipper = self.intp(et)
        p_launch = bodies['Earth'].GetPosition(et)
        # print('Clipper at {}, launch site at {}, launch pos at {}'.format(p_clipper, p_launch, self.launch_pos))
        self.position_actor.SetPosition(p_clipper)

        if not self.has_launched and et >= self.launch:
            self.has_launched = True
            self.orbit_actor = make_curve_actor([p_clipper, p_clipper], [0.5, 0.5, 0])

        if self.has_launched:
            # print('adding one point to curve')
            self.orbit_actor = add_point_to_curve(self.orbit_actor, p_clipper)
            self.orbit_actor.Modified()

    def GetPosition(self, et=clock):
        return self.intp(et)

def copy_transformation(actor_from, actor_to):
    actor_to.SetUserMatrix(actor_from.GetUserMatrix())

def distance(body1, body2='Sun', et=clock, frame=default_frame):
    p1 = bodies[body1].GetPosition(et, frame)
    p2 = bodies[body2].GetPosition(et, frame)
    return np.linalg.norm(p1-p2)

def reference_body(name):
    if name == 'Sun' or name in planet_names:
        return 'Sun'
    elif name in jupiter_moon_names:
        return 'Jupiter'
    elif name == 'Moon':
        return 'Earth'
    elif name == 'Clipper':
        return 'Earth'
    else:
        print('unrecognized body name: {}'.format(name))

def vantage_point(target_name, et=clock):
    global clipper_orbit

    if target_name == 'Clipper':
        print('handling special case of Clipper')
        p_now = clipper_orbit.GetPosition()
        p_old = clipper_orbit.GetPosition(et-2*hour)
        return { 'pos': p_old, 'focal': p_now, 'up': [0,0,1] }
    elif target_name == 'Clipper-Venus':
        return { 'pos': clipper_orbit.GetPosition(), 'focal': bodies['Venus'].GetPosition(), 'up': [0, 0, 1] }
    elif target_name == 'Clipper-Earth':
        return { 'pos': clipper_orbit.GetPosition(), 'focal': bodies['Earth'].GetPosition(), 'up': [0, 0, 1] }
    elif target_name == 'Clipper-Mars':
        return { 'pos': clipper_orbit.GetPosition(), 'focal': bodies['Mars'].GetPosition(), 'up': [0, 0, 1] }
    elif target_name == 'Clipper-Jupiter':
        return { 'pos': clipper_orbit.GetPosition(), 'focal': bodies['Jupiter'].GetPosition(), 'up': [0, 0, 1] }
    elif target_name == 'Clipper-Europa':
        return { 'pos': clipper_orbit.GetPosition(), 'focal': bodies['Europa'].GetPosition(), 'up': [0, 0, 1] }
    elif target_name == 'ecliptic':
        target_pos = bodies['Sun'].GetPosition()
        return { 'pos': target_pos + np.array([0, 0, au_to_km(10)]), 'focal': target_pos, 'up': [0, 1, 0] }

    main_body = reference_body(target_name)
    main_body_pos = bodies[main_body].GetPosition()
    target_pos = bodies[target_name].GetPosition()
    if target_name.upper() == 'SUN':
        print('handling special case of Sun')
        dir = target_pos - bodies['Earth'].GetPosition()
    else:
        dir = target_pos - main_body_pos
    dir /= np.linalg.norm(dir)
    radius = bodies[target_name].GetMaxRadius()
    return { 'pos': target_pos - 5*radius*dir, 'focal': target_pos, 'up': [0,0,1] }

def clock_as_str(cl=clock):
    return spice.timout(cl, 'AP:MN:SC AMPM Month DD, YYYY')

sun_mercury_distance = distance('Mercury', 'Sun', init_time)
sun_radius = bodies['Sun'].GetMaxRadius()

def get_renderer():
    return render_window.GetRenderers().GetFirstRenderer()

def get_camera():
    return get_renderer().GetActiveCamera()

def camera_description(camera):
    pos = camera.GetPosition()
    pos_str = '({:06.2f}, {:06.2f}, {:06.2f})'.format(pos[0], pos[1], pos[2])
    foc = camera.GetFocalPoint()
    foc_str = '({:06.2f}, {:06.2f}, {:06.2f})'.format(foc[0], foc[1], foc[2])
    up = camera.GetViewUp()
    up_str = '({:06.2f}, {:06.2f}, {:06.2f})'.format(up[0], up[1], up[2])
    return 'p={}, f={}, u={}'.format(pos_str, foc_str, up_str)

def make_rings():
    global ring_actors

    ring_reader = vtk.vtkJPEGReader()
    ring_reader.SetFileName(texture_files['Rings'])
    ring_texture = vtk.vtkTexture()
    ring_texture.SetInputConnection(ring_reader.GetOutputPort())

    minR = saturn_rings['C'][0]
    maxR = saturn_rings['A'][1]
    for ring_name in ['C', 'B', 'A']:
        inner, outer = saturn_rings[ring_name]
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

def make_ecliptic_plane():
    global all_actors
    global ecliptic_actors

    ecliptic = vtk.vtkRegularPolygonSource()
    ecliptic.SetCenter(0,0,0)
    ecliptic.SetNumberOfSides(500)
    ecliptic.SetRadius(distance('Neptune', 'Sun', init_time))
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(ecliptic.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0,0.2,1)
    actor.GetProperty().SetOpacity(0.05)
    actor.GetProperty().SetAmbient(1)
    all_actors['ecliptic plane'] = actor
    ecliptic_actors.append(actor)
    frame = vtk.vtkRegularPolygonSource()
    frame.SetCenter(0,0,0)
    frame.SetNumberOfSides(500)
    frame.SetRadius(distance('Neptune', 'Sun', init_time))
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
    all_actors['ecliptic contour'] = f_actor
    ecliptic_actors.append(f_actor)

def make_orbits():
    global orbit_actors
    global all_actors

    body_groups = [ [ planet_names, day ], [ jupiter_moon_names, hour ] ]
    for names, period in body_groups:
        for name in names:
            print(' * {}'.format(name))
            body = bodies[name]
            points = []
            for t in np.arange(init_time, init_time + min(1.01*orbital_period[name], final_time-init_time), period):
                points.append(body.GetPosition(t))
            actor = make_curve_actor(points, planet_colors[name], 2)
            orbit_actors[name] = actor
            all_actors['orbit of ' + name] = actor
    print(' * Clipper')
    intp = interpolate.make_interp_spline(x=etclipper, y=clipper_full)
    points = []
    for t in np.arange(init_time, final_time, hour):
        points.append(intp(t))
    actor = make_curve_actor(points, [0.3, 0, 0], 0.25)
    orbit_actors['Clipper'] = actor
    all_actors['orbit of Clipper'] = actor

def focal_vector():
    cam = get_camera()
    p = np.array(cam.GetPosition())
    f = np.array(cam.GetFocalPoint())
    return f-p

def print_camera():
    cam = get_camera()
    settings = {
        'from' : cam.GetPosition(),
        'to'   : cam.GetFocalPoint(),
        'up'   : cam.GetViewUp(),
        'clip' : cam.GetClippingRange(),
        'view' : cam.GetViewAngle()
    }
    txt = json.dump(settings)
    print(txt)

def read_camera(str):
    cam = get_camera()
    settings = json.load(str)
    if 'from' in settings:
        cam.SetPosition(settings['from'])
    if 'to' in settings:
        cam.SetFocalPoint(settings['to'])
    if 'up' in settings:
        cam.SetViewUp(settings['up'])
    if 'clip' in settings:
        cam.SetClippingRange(settings['clip'])
    if 'view' in settings:
        cam.SetViewAngle(settings['view'])

def focal_distance():
    return np.linalg.norm(focal_vector())

# Callback function for the GUI slider
def slider_callback(obj, event):
    global time_step
    value = obj.GetRepresentation().GetValue()
    # number of days per second
    # refresh_rate * step = value*day
    time_step = float(value*hour)/refresh_rate
    print('time_step is now {}'.format(time_step))

def planet_slider_callback(obj, event):
    global ring_actors
    planet_scale = obj.GetRepresentation().GetValue()
    for name in bodies:
        body = bodies[name]
        if name == 'Sun' and planet_scale > sun_max:
            body.ScaleActor(sun_max)
        else:
             body.ScaleActor(planet_scale)
    # adjust Saturn rings accordingly
    for actor in ring_actors:
        actor.SetScale(planet_scale, planet_scale, planet_scale)

def time_slider_callback(obj, event):
    global clock
    global clipper_orbit
    global render_window
    global video_mode

    current = clock
    et = obj.GetRepresentation().GetValue()
    clock = et*hour + etb[0]

    if video_mode:
        video_mode = False
        print('Video mode turned OFF after date change')
    

    renderer = get_renderer()
    renderer.RemoveActor(clipper_orbit.orbit_actor)
    renderer.RemoveActor(clipper_orbit.position_actor)
    # renderer.RemoveActor(clipper_orbit.launch_orbit_actor)
    # renderer.RemoveActor(clipper_orbit.launch_site_curve_actor)
    clipper_orbit = ClipperOrbit()
    renderer.AddActor(clipper_orbit.orbit_actor)
    renderer.AddActor(clipper_orbit.position_actor)
    # renderer.AddActor(clipper_orbit.launch_orbit_actor)
    # renderer.AddActor(clipper_orbit.launch_site_curve_actor)

# Handles the key press to change camera focal point
def key_pressed_callback(obj, event):
    global focus_planet
    global show_ecliptic
    global show_orbits
    global orbit_actors
    global ecliptic_actors
    global resolution
    global prev_resolution
    global do_tether
    global planet_focus_changed
    global zoom_factor
    global video_mode

    new_key = obj.GetKeySym()
    print(f'newkey is {new_key}')
    if new_key.isdigit():
        key = int(new_key)
        if key <= 8 and all_body_names[key] != focus_planet:
            print('changing focus from {} to {}'.format(focus_planet, all_body_names[key]))
            focus_planet = all_body_names[key]
            planet_focus_changed = True
            saved_relative_position = False
            do_tether = True
            cam_setting = vantage_point(focus_planet)
            get_camera().SetPosition(cam_setting['pos'])
            get_camera().SetFocalPoint(cam_setting['focal'])
            get_camera().SetViewUp(cam_setting['up'])
            text_actor_body.SetInput('Current Focus: {}'.format(focus_planet))

        # Overrides the default vtk behavior for keypress '3'
        if int(new_key) == 3:
            render_window.StereoRenderOn()
    elif new_key == 'minus':
        focus_planet = 'N/A'
        planet_focus_changed = False
        text_actor_body.SetInput('Current Focus: None')
        do_tether = False
    elif new_key == 'o' or new_key == 'O':
        show_orbits = not show_orbits
        if show_orbits:
            print('now showing orbits')
            for name in orbit_actors:
                if clock < etarrive and name in jupiter_moon_names:
                    continue
                print('added orbit of {}'.format(name))
                print('there are {} vertices in this orbit'.format(orbit_actors[name].GetMapper().GetInput().GetPoints().GetNumberOfPoints()))
                get_renderer().AddActor(orbit_actors[name])
        else:
            print('now hiding orbits')
            for name in orbit_actors:
                get_renderer().RemoveActor(orbit_actors[name])
        render_window.Render()
    elif new_key == 'f':
        if render_window.GetFullScreen():
            render_window.FullScreenOff()
        render_window.Render()
    elif new_key == 'v':
        video_mode = not video_mode
        if video_mode:
            print('Video mode turned ON')
        else:
            print('Video mode turned OFF')
        return
    elif new_key == 'exclam' or new_key == 'at' or new_key == 'numbersign' \
         or new_key == 'dollar' or new_key == 'asciicircum' \
         or new_key == 'plus' or new_key == 'equal' or new_key == 'ampersand' \
         or new_key == 'n' or new_key == 'a' or new_key == 'b' or new_key == 'c' \
         or new_key == 'd' or new_key == 'g' or new_key == 'semicolon' \
         or new_key == 'quoteright' or new_key == 'k':
        new_focus_planet = focus_planet
        if new_key == 'exclam':
            new_focus_planet = 'Io'
        elif new_key == 'at':
            new_focus_planet = 'Europa'
        elif new_key == 'numbersign':
            new_focus_planet = 'Ganymede'
        elif new_key == 'dollar':
            new_focus_planet = 'Callisto'
        elif new_key == 'asciicircum':
            new_focus_planet = 'Moon'
        elif new_key == 'plus':
            new_focus_planet = 'Clipper Forward'
        elif new_key == 'equal':
            new_focus_planet = 'Clipper Backward'
        elif new_key == 'ampersand':
            new_focus_planet = 'Clipper'
        elif new_key == 'n':
            new_focus_planet = 'ecliptic'
        elif new_key == 'a':
            new_focus_planet = 'Clipper-Earth'
        elif new_key == 'b':
            new_focus_planet = 'Clipper-Mars'
        elif new_key == 'c':
            new_focus_planet = 'Clipper-Jupiter'
        elif new_key == 'd':
            new_focus_planet = 'Clipper-Europa'
        elif new_key == 'g':
            new_focus_planet = 'Clipper-Venus'
        elif new_key == 'semicolon':
            zoom_factor /= 2.0
        elif new_key == 'quoteright':
            zoom_factor *= 2.0
        elif new_key == 'k':
            print_camera()
            return
        if new_focus_planet != focus_planet:
            focus_planet = new_focus_planet.split()[0]
            planet_focus_changed = True
            saved_relative_position = False
            do_tether = True
            camera = get_camera()
            cam_setting = vantage_point(focus_planet)
            camera.SetPosition(cam_setting['pos'])
            camera.SetFocalPoint(cam_setting['focal'])
            camera.SetViewUp(cam_setting['up'])
            text_actor_body.SetInput('Current Focus: {}'.format(focus_planet))
        camera = get_camera()
        if zoom_factor != 1:
            camera.Zoom(zoom_factor)
            zoom_factor = 1
            print('camera view angle: {}'.format(camera.GetViewAngle()))
    elif new_key == 't':
        # tether mode
        if not do_tether:
            do_tether = True
            print('tethering activated')
        else:
            do_tether = False
            saved_relative_position = None
            print('tethering deactivated')
    elif new_key == 'h':
        print('List of keyboard commands:')
        print('\'0\', \'1\', \'2\', ..., \'8\': Point the camera to the Sun, Mercury, Venus, ..., Neptune and track')
        print('\'!\', \'@\', \'#\', \'$\':  Point the camera to Jupiter\'s Galilean moons')
        print('\'^\':                     Point the camera to the Moon')
        print('\'&\':                     Point the camera to Clipper\'s position')
        print('\'-\':                     Turn off tracking')
        print('\'f\':                     Toggle full screen on and off')
        print('\'a\', ..., \'d\':           Point camera from Clipper to Earth, Mars, Jupiter, Europa')
        print('\'Ctrl-9\':                Point camera from Clipper to the Moon')
        print('\'Ctrl-!\', \'Ctrl-@\':      Point camera from Clipper to Jupiter\'s moon')
        print('\'t\':                     Toggle tethering on and off')
        print('\'n\':                     Point the camera to the Sun\'s North pole')
        print('\'v\':                     Toggle video recording mode')
        print('\'o\' (\'O\'):               Turn on/off depiction of planets\' orbits')
        print('\'h\':                     Print this information')
    else:
        print('unrecognized entered key is {}'.format(new_key))

def window_resized_callback(obj, event):
    global render_window
    global resize_image
    global render_x, render_y
    global resolution
    width, height = render_window.GetSize()
    if resolution[0] != width or resolution[1] != height:
        resolution[0] = width
        resolution[1] = height
        resize_image.SetOutputDimensions(resolution[0], resolution[1], 0)
        resize_image.Update()
        print('texture has been resized')
        render_window.Render()

def camera_modif_callback(caller, event):
    global in_progress
    global cam_changed
    if not in_progress:
        cam_changed = True

def interaction_start_callback(caller, event):
    global paused
    global do_tether

    if do_tether:
        paused = True

def interaction_end_callback(caller, event):
    global paused
    global do_tether

    if do_tether:
        paused = False

def interaction_callback(caller, event):
    global paused

    cam = render_window.GetRenderers().GetFirstRenderer().GetActiveCamera()

def save_frame(window, clock):
    global video_basename
    global args
    # ---------------------------------------------------------------
    # Save current contents of render window to PNG file
    # ---------------------------------------------------------------
    file_name = video_basename + str(int(clock)) + ".png"
    image = vtk.vtkWindowToImageFilter()
    image.SetInput(window)
    png_writer = vtk.vtkPNGWriter()
    png_writer.SetInputConnection(image.GetOutputPort())
    png_writer.SetFileName(file_name)
    window.Render()
    png_writer.Write()
    if args.verbose:
        print(file_name + " has been successfully exported")

# Handles timers from the interactive render window and updates planet positions
class TimerCallback:
    def __init__(self):
        global scoped_lock
        global clock
        self.ring_actors = []
        self.arrow_actor = {}
        self.body_actors = {}
        self.body_point_actors = {}
        self.body_tag_actors = {}
        self.camera = None
        self.text = None
        self.time_display = None

    def execute(self, obj, event):
        global focus_planet
        global time_step
        global render_window
        global clock
        global planet_focus_changed
        global do_tether
        global saved_relative_position
        global cam_changed
        global in_progress
        global paused
        global clipper_orbit
        global launch_orbit
        global time_widget
        global show_clipper_orbit
        global zoom_factor
        global video_mode

        cam = render_window.GetRenderers().GetFirstRenderer().GetActiveCamera()
        if paused:
            return

        timestr = clock_as_str(clock)
        if video_mode:
            save_frame(render_window, clock)
        if clock >= final_time:
            return
        current_focal_vector = focal_vector()
        # Update camera for currently selected planet
        if focus_planet != 'N/A':
            # tether the camera position to the moving planet to maintain the relative position
            # that it has at the beginning or after a user change to the camera
            show_clipper_orbit = True
            in_progress = True
            if focus_planet.split()[0] == 'Clipper':
                new_body_pos = clipper_orbit.GetPosition(clock)
            elif focus_planet == 'Clipper-Earth':
                new_body_pos = bodies['Earth'].GetPosition(clock)
            elif focus_planet == 'Clipper-Mars':
                new_body_pos = bodies['Mars'].GetPosition(clock)
            elif focus_planet == 'Clipper-Jupiter':
                new_body_pos = bodies['Jupiter'].GetPosition(clock)
            elif focus_planet == 'Clipper-Europa':
                new_body_pos = bodies['Europa'].GetPosition(clock)
            elif focus_planet == 'Clipper-Venus':
                new_body_pos = bodies['Venus'].GetPosition(clock)
            elif focus_planet != 'ecliptic':
                new_body_pos = bodies[focus_planet].GetPosition(clock)
            else:
                new_body_pos = bodies['Sun'].GetPosition(clock)
            self.camera.SetFocalPoint(new_body_pos)
            in_progress = False
            if focus_planet == 'Clipper-Earth' or \
                focus_planet == 'Clipper-Mars' or \
                focus_planet == 'Clipper-Jupiter' or \
                focus_planet == 'Clipper-Europa' or \
                focus_planet == 'Clipper-Venus':
                show_clipper_orbit = False
                planet_name = focus_planet[8:]
                clipper_new_pos = clipper_orbit.GetPosition(clock)
                do_tether = False
                in_progress = True
                self.camera.SetPosition(clipper_new_pos)
                in_progress = False
                self.text.SetInput('Current Focus: {}'.format(focus_planet))
            else:
                cam_pos = np.array(self.camera.GetPosition())
                self.text.SetInput('Current Focus: ' + focus_planet)
                if do_tether:
                    if saved_relative_position is None or cam_changed:
                        saved_relative_position = current_focal_vector
                        cam_changed = False
                    in_progress = True
                    self.camera.SetPosition(new_body_pos - saved_relative_position)
                    in_progress = False

        for name in all_body_names:
            body = bodies[name]
            # Update position and rotation
            self.body_actors[name] = body.MoveActor(clock)
            self.body_point_actors[name] = body.MovePointActor(clock)
            if name == 'Saturn':
                for ra in self.ring_actors:
                    ra.SetPosition(self.body_actors['Saturn'].GetPosition())
                    copy_transformation(body.GetActor(), ra)

        # extend Clipper's orbit
        clipper_orbit.Update(clock)


        self.time_display.SetInput(timestr)
        time_widget.GetRepresentation().SetValue((clock-etb[0])/hour)
        obj.GetRenderWindow().GetRenderers().GetFirstRenderer().RemoveActor(clipper_orbit.orbit_actor)
        if show_clipper_orbit:
            obj.GetRenderWindow().GetRenderers().GetFirstRenderer().AddActor(clipper_orbit.orbit_actor)
        obj.GetRenderWindow().Render()

        if not paused:
            clock += time_step

def make_slider(minval, maxval, init_val, title, xmin, xmax, y=0.1, format="%.2f"):
    global render_window

    slider = vtk.vtkSliderRepresentation2D()
    slider.DebugOn()
    slider.SetMinimumValue(minval)
    slider.SetMaximumValue(maxval)
    slider.GetSliderProperty().SetLineWidth(0.5)
    slider.GetTubeProperty().SetLineWidth(0.5)
    slider.GetLabelProperty().SetFontSize(1)
    slider.SetLabelHeight(0.015)
    slider.GetLabelProperty().BoldOff()
    slider.SetValue(init_val)
    slider.SetTitleText(title)
    slider.GetTitleProperty().SetFontSize(1)
    slider.GetTitleProperty().BoldOff()
    slider.SetTitleHeight(0.02)
    slider.SetTubeWidth(0.005)
    slider.SetLabelFormat(format)  # float format
    slider.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider.GetPoint1Coordinate().SetValue(xmin, y)
    slider.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider.GetPoint2Coordinate().SetValue(xmax, y)

    widget = vtk.vtkSliderWidget()
    widget.DebugOn()
    widget.SetInteractor(render_window.GetInteractor())
    widget.SetRepresentation(slider)
    widget.SetAnimationModeToJump()
    widget.EnabledOn()

    return [slider, widget]

def main():
    global render_window
    global focus_planet
    global refresh_rate
    global time_step
    global resize_image
    global text_actor_body
    global show_ecliptic
    global show_orbits
    global show_clipper_orbit
    global resolution
    global planet_min_max
    global bodies
    global orbit_actors
    global ecliptic_actors
    global ring_actors
    global body_actors
    global body_point_actors
    global body_tag_actors
    global all_actors
    global clipper_orbit
    global time_widget
    global zoom_factor
    global do_shadows

    # Timer Frequency
    refresh_rate = 5  # Hz
    timer_period = 1000/refresh_rate  # in ms
    time_step = minute/refresh_rate

    zoom_factor = 1

    show_clipper_orbit = True

    print('time step is {}'.format(time_step))

    # Render size set to HD Display
    focus_planet = 'Sun'  # Initial focus planet is the sun

    body_mappers = {}
    body_actors = {}
    body_point_actors = {}
    body_tag_actors = {}
    arrows_actor = {}
    ring_actors = []
    orbit_actors = {}
    ecliptic_actors = []
    all_actors = {}

    # Create render window
    render_window = vtk.vtkRenderWindow()
    if not do_shadows:
        render_window.SetNumberOfLayers(2)
    if resolution[0] > 0 and resolution[1] > 0:
        render_window.SetSize(resolution[0], resolution[1])  # Set the window size
    else:
        resolution = render_window.GetScreenSize()
        render_window.SetSize(resolution)
    render_window.SetWindowName("Solar System")
    render_window.StereoRenderOff()

    print('Creating Jupiter\'s rings...')
    make_rings()
    print('...done.')
    print('Creating ecliptic plane...')
    make_ecliptic_plane()
    print('...done.')
    print('Creating bodies\' orbit curves...')
    make_orbits()
    print('...done.')
    print('Creating Clipper object...')
    clipper_orbit = ClipperOrbit()
    print('...done.')

    idx = 0
    print('Initializing VTK actors...')
    for name in all_body_names:
        body = bodies[name]

        # Create Actors and add textures, rotations, and orbit position
        body_actors[name] = body.MoveActor(clock)
        body_actors[name].GetProperty().SetAmbient(ambient)
        body_actors[name].GetProperty().SetSpecular(0.01)
        body_actors[name].GetProperty().SetDiffuse(0.65)
        all_actors['Body ' + name] = body_actors[name]

        body_point_actors[name] = body.MovePointActor(clock)
        all_actors['Point ' + name] = body_point_actors[name]

        if name == 'Saturn':
            print('moving rings')
            for actor in ring_actors:
                actor.SetPosition(body.GetPosition())
                copy_transformation(body.GetActor(), actor)
                actor.GetProperty().SetAmbient(ambient)
    print('...done.')

    # Create Lighting
    sunlight1 = vtk.vtkLight()
    sunlight1.SetPosition(0, 0, 0)
    sunlight1.SetFocalPoint(0.0000001*scale, 0, 0)  # Right hemisphere lighting
    sunlight1.SetConeAngle(180)
    sunlight1.SetPositional(1)

    sunlight2 = vtk.vtkLight()
    sunlight2.SetPosition(0, 0, 0)
    sunlight2.SetFocalPoint(-0.0000001*scale, 0, 0)  # Left hemisphere lighting
    sunlight2.SetConeAngle(180)
    sunlight2.SetPositional(1)

    # Add Background Starry Image
    if not do_shadows:
        print('Loading background image...')
        jpeg_reader = vtk.vtkJPEGReader()
        jpeg_reader.SetFileName(texture_files['Stars'])
        jpeg_reader.Update()

        # Resize Image to Match Render Window Size
        resize_image = vtk.vtkImageResize()
        resize_image.SetInputData(jpeg_reader.GetOutput())
        resize_image.SetResizeMethod(0)
        if resolution[0] > 0 and resolution[1] > 0:
            resize_image.SetOutputDimensions(resolution[0], resolution[1], 0)
        resize_image.Update()

        background_mapper = vtk.vtkImageMapper()
        background_mapper.SetInputData(resize_image.GetOutput())
        background_actor = vtk.vtkActor2D()
        background_actor.SetMapper(background_mapper)
        all_actors['background'] = background_actor
        background_actor.DebugOn()
        print('...done')

    # Light up the sun no matter what
    sun_props = vtk.vtkProperty()
    sun_props.SetLighting(0)
    body_actors['Sun'].SetProperty(sun_props)

    # Create a renderer and add the actors to the scene
    print('Rendering Planets...')
    body_renderer = vtk.vtkRenderer()
    # body_renderer.UseShadowsOn()
    body_renderer.SetBackground(0, 0, 0)  # black background
    body_renderer.AddLight(sunlight1)
    body_renderer.AddLight(sunlight2)
    print('Adding all actors to renderer...')
    for name in all_body_names:
        body_renderer.AddActor(body_actors[name])
        body_renderer.AddActor(body_point_actors[name])
    for actor in ring_actors:
        # print('adding one ring actor {}'.format(actor))
        body_renderer.AddActor(actor)
    if show_ecliptic:
        for actor in ecliptic_actors:
            body_renderer.AddActor(actor)
    if show_orbits:
        print('orbit_actors={}'.format(orbit_actors))
        for body, actor in orbit_actors.items():
            print('body={}, actor={}'.format(body, actor))
            body_renderer.AddActor(actor)
    if show_clipper_orbit:
        body_renderer.AddActor(clipper_orbit.orbit_actor)
    body_renderer.AddActor(clipper_orbit.position_actor)
    print('...done.')

    # Create the camera focused on the Sun
    cam1 = vtk.vtkCamera()
    cam1.SetFocalPoint(0, 0, 0)  # Sun Center
    cam1.SetPosition(-distance('Mercury', 'Sun'), 0, 0.5*distance('Mercury', 'Sun'))
    cam1.Elevation(0)  # look at off-angle system (0 is top down)
    cam1.Azimuth(0)
    cam1.SetViewUp(0,0,1)
    cam1.SetClippingRange(100, 100000000000)
    print('camera view angle is: {}'.format(cam1.GetViewAngle()))

    body_renderer.SetActiveCamera(cam1)
    if not do_shadows:
        body_renderer.SetLayer(1)  # In front of background image

    if do_shadows:
        print('adding shadows')
        render_window.SetMultiSamples(0)
        shadows = vtk.vtkShadowMapPass()
        sequence = vtk.vtkSequencePass()
        passes = vtk.vtkRenderPassCollection()
        passes.AddItem(shadows.GetShadowMapBakerPass())
        passes.AddItem(shadows)
        sequence.SetPasses(passes)
        camera_pass = vtk.vtkCameraPass()
        camera_pass.SetDelegatePass(sequence)

        gl_renderer = body_renderer
        gl_renderer.SetPass(camera_pass)


    # Add useful text to the renderer
    # Render Title
    print('Setting text actors...')
    text_actor_title = vtk.vtkTextActor()
    text_actor_title.SetInput("The Solar System")
    text_actor_title.SetPosition(40, resolution[1] - 90)
    text_actor_title.GetTextProperty().SetFontSize(30)
    text_actor_title.GetTextProperty().SetColor(1.0, 1.0, 1.0)
    all_actors['title'] = text_actor_title

    # Current planet of focus
    text_actor_body = vtk.vtkTextActor()
    text_actor_body.SetInput('Current Focus: Sun')
    text_actor_body.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    text_actor_body.GetPositionCoordinate().SetValue(.80, .95)
    text_actor_body.GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    text_actor_body.GetPosition2Coordinate().SetValue(.99, .95)
    text_actor_body.GetTextProperty().SetFontSize(30)
    text_actor_body.GetTextProperty().SetColor(1.0, 1.0, 1.0)
    all_actors['focus body'] = text_actor_body

    time_text = vtk.vtkTextActor()
    time_text.SetInput('{}'.format(spice.timout(init_time, 'AP:MN:SC AMPM Month DD, YYYY')))
    time_text.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    time_text.GetPositionCoordinate().SetValue(0.01, 0.95)
    time_text.GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    time_text.GetPosition2Coordinate().SetValue(0.25, 0.95)
    time_text.GetTextProperty().SetFontSize(30)
    time_text.GetTextProperty().SetColor(1.0, 1.0, 1.0)
    all_actors['time display'] = time_text
    print('...done')

    body_renderer.AddActor2D(text_actor_body)
    body_renderer.AddActor2D(time_text)

    all_actors = body_renderer.GetActors()
    all_actors.InitTraversal()
    while True:
        actor = all_actors.GetNextActor()
        if actor:
            actor.DebugOn()
        else:
            break

    render_window.AddRenderer(body_renderer)

    # Create a background renderer
    if not do_shadows:
        background_renderer = vtk.vtkRenderer()
        background_renderer.DebugOn()
        background_renderer.SetLayer(0)
        background_renderer.InteractiveOff()
        background_renderer.AddActor(background_actor)
        render_window.AddRenderer(background_renderer)
        render_window.DebugOn()

    # visualize ecliptic plane
    if show_ecliptic:
        for actor in ecliptic_actors:
            body_renderer.AddActor(actor)

    # Set-up interactor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.DebugOn()
    style = vtk.vtkInteractorStyleTrackballCamera()
    style.AutoAdjustCameraClippingRangeOff()
    style.DebugOn()
    render_window_interactor.SetInteractorStyle(style)
    render_window_interactor.SetRenderWindow(render_window)
    render_window_interactor.Initialize()
    render_window.Render()

    render_window.AddObserver('ModifiedEvent', window_resized_callback)

    '''
    def make_slider(minval, maxval, init_val, title, xmin, xmax, y=0.1, format="%.2f"):
    '''

    print('Creating sliders...')
    slider, day_widget = make_slider(minval=0, maxval=168, init_val=0.1, title='Hours per second', xmin=0.02, xmax=0.22)
    day_widget.AddObserver("InteractionEvent", slider_callback)

    # Widget to control the scale of the planets
    slider2, planet_widget = make_slider(planet_min_max[0], planet_min_max[1], 1, 'Bodies scaling factor', .24, .44)
    planet_widget.AddObserver("InteractionEvent", planet_slider_callback)

    time_slider, time_widget = make_slider(0, (etb[1]-etb[0])/hour, 0, 'ET in hours since {}'.format(clock_as_str(etb[0])), 0.46, 0.96)
    time_widget.AddObserver('InteractionEvent', time_slider_callback)

    # Sign up to receive TimerEvent
    cb = TimerCallback()
    cb.body_actors = body_actors
    cb.ring_actors = ring_actors
    cb.camera = cam1
    cb.camera.AddObserver(vtk.vtkCommand.ModifiedEvent, camera_modif_callback)
    cb.text = text_actor_body
    cb.time_display = time_text

    render_window_interactor.AddObserver('TimerEvent', cb.execute)
    render_window_interactor.AddObserver('KeyPressEvent', key_pressed_callback)
    render_window_interactor.AddObserver('InteractionEvent', interaction_callback)
    render_window_interactor.AddObserver('StartInteractionEvent', interaction_start_callback)
    render_window_interactor.AddObserver('EndInteractionEvent', interaction_end_callback)
    render_window_interactor.CreateRepeatingTimer(int(timer_period))  # ms between calls
    render_window_interactor.EnableRenderOn()
    print('...done')
    print('Starting interactive loop...')
    render_window_interactor.Initialize()
    render_window.Render()
    render_window_interactor.Start()

if __name__ == '__main__':
    global show_ecliptic
    global show_orbits
    global planet_min_max
    global resolution
    global do_record
    global filebase
    global do_shadows
    global video_basename

    parser = argparse.ArgumentParser(
            description='Visualize Clipper Spacecraft and Solar System using NAIF Data',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-e', '--ecliptic', action='store_true', help='Show ecliptic plane')
    parser.add_argument('-o', '--orbits', action='store_true', help='Show planets\' orbits')
    parser.add_argument('-r', '--resolution', metavar='int', default=[2560, 1440], type=int, nargs=2, help='window resolution (-1: full screen)')
    parser.add_argument('-s', '--sun', metavar='float', default='40', type=float, help='Max sun radius (in multiples of Mercury semi-major axis)')
    parser.add_argument('-p', '--planet', metavar='float', default=[1, 2000], type=float, nargs=2, help='Range of scaling factors for planets')
    parser.add_argument('-x', '--scale', metavar='int', default=10000, type=int, help='Uniform space scaling factor')
    parser.add_argument('-t', '--time', metavar='float', default=[-1, -1], type=float, nargs=2, help='Temporal range to consider')
    parser.add_argument('-v', '--video', metavar='filename', default='ephemeris', type=str, help='Filename to store each individual frame of the animation')
    parser.add_argument('-c', '--camera', metavar='json', type=str, help='Json string containing camera setting')
    parser.add_argument('--shadow', action='store_true', help='Use shadow map pass')
    parser.add_argument('--verbose', action='store_true', help='Toggle verbose mode')
    args = parser.parse_args()

    show_ecliptic = args.ecliptic
    show_orbits = args.orbits
    scale = args.scale
    sun_max = args.sun
    planet_min_max = args.planet
    resolution = args.resolution
    if args.video is not None:
        do_record = True
        filebase = args.video
    if args.camera is not None:
        load_camera(args.camera)
    if args.time is not None:
        selected_interval = args.time
        print('selected interval = {}'.format(selected_interval))
    if args.shadow:
        do_shadows = args.shadow
    else:
        do_shadows = False

    video_basename = args.video

    main()
