#########################################################################
#
# CS530 - Final Project
#
# Author: Drew Brandsen
# due: 04/29/2018
#
# Textures courtesy of:
#   https://www.solarsystemscope.com/textures/
# Background image courtesy of:
#   http://www.wallpixa.com/download-space-hd-picture-wallpaper-background-starry-sky-stars-night-116961-1920x1080-3548/
#
# Function: solarsystem.py
#     input:
#           Textures for each planet (Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune) and backdrop
#     output:
#           Solar System starting from current time and moving forward 0.15 days per second
#     interaction:
#           Pan and zoom, use keys 1-9 to change the camera focal point to each planet (+ sun)
#           Use the slider bar to change animation speed from stopped to 10 earth-days per second
#
#########################################################################
#! python2

import sys
import vtk
import math
from datetime import datetime
import argparse

# spatial scaling
scale=10000
planet_scale = 100
ambient = 0.1
bodies = {}
show_arrows = False
sun_max = 2
prev_resolution = [-1, -1]

all_actors = {}


files = {   'Sun': '2k_sun.jpg', 
            'Mercury': '2k_mercury.jpg', 
            'Venus': '2k_venus.jpg', 
            'Earth': '2k_earth_daymap.jpg', 
            'Mars': '2k_mars.jpg', 
            'Jupiter': '2k_jupiter.jpg', 
            'Saturn': '2k_saturn.jpg', 
            'Uranus': '2k_uranus.jpg', 
            'Neptune': '2k_neptune.jpg', 
            'Stars': 'starry_background.jpg',
            'Rings': 'ring_colors2d.jpg' }

planet_radius = [ 0.0046504673, 0.000016308387, 0.000040453784, 0.000042587505, 0.000022657408, 0.00046732617, 0.00038925688, 0.00016953450, 0.00016458790]
day_length = [24.0, 58.25, 116.25, 1.0, 1.5, 10.0, 10.75, 17.25, 16.0]
n_arr = [(0, 0), (48.3313, 3.24587E-5), (76.6799, 2.46590E-5), (0, 0), (49.5574, 2.11081E-5),
         (100.4542, 2.76854E-5), (113.6634, 2.38980E-5), (74.0005, 1.3978E-5), (131.7806, 3.0173E-5)]
i_arr = [(0, 0), (7.0047, 5.00E-8), (3.3946, 2.75E-8), (0, 0), (1.8497, -1.78E-8),
         (1.3030, -1.557E-7), (2.4886, -1.081E-7), (0.7733, 1.9E-8), (1.7700, -2.55E-7)]
w_arr = [(0, 0), (29.1241, 1.01444E-5), (54.8910, 1.38374E-5), (102.9, 4.70935E-5), (286.5016, 2.92961E-5),
         (273.8777, 1.64505E-5), (339.3939, 2.97661E-5), (96.6612, 3.0565E-5), (272.8461, -6.027E-6)]
a_arr = [(0, 0), (0.387098, 0), (0.72333, 0), (1, 0), (1.523688, 0),
         (5.20256, 0), (9.55475, 0), (19.18171, -1.55E-8), (30.05826, 3.313E-8)]
e_arr = [(0, 0), (0.205635, 5.59E-10), (0.006773, -1.302E-9), (0.016709, - 1.151E-9), (0.093405, 2.516E-9),
         (0.048498, 4.469E-9), (0.055546, -9.499E-9), (0.047318, 7.45E-9), (0.008606, 2.15E-9)]
m_arr = [(0, 0), (168.6562, 4.092334436), (48.0052, 1.602130224), (356.0470, 0.985600258), (18.6021, 0.524020776),
             (19.8950, 0.0830853001), (316.9670, 0.0334442282), (142.5905, 0.011725806), (260.2471, 0.005995147)]
axial_tilt = [-7.25, -0.03, -2.64, -23.44, -25.19, -3.13, -26.73, -82.23, -28.32]
years = [0, 87.969, 224.7, 365.2564, 687.011, 4332.59, 10765.866, 30688.5, 60182]
body_names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
saturn_rings = { 'D': (66900, 74510), 'C': (74658, 92000), 'B': (92000, 117580), 'Cassini': (117580, 122170), \
                 'A': (122170, 136775), 'Roche': (136775, 139380), 'F': (140183, 140680) }
          
def au_to_miles(au):
    return au * 92955807
    
def miles_to_au(miles):
    return miles / 92955807
    
def au_to_km(au):
    return au * 149597871
    
def km_to_au(km):
    return km / 149597871
    
def deg_to_rad(deg):
    return deg/180.*math.pi
    
def print_camera_info(camera):
    print('Camera info:\n*position:{0}\n*focal point:{1}\n*clipping range:{2}\n*view up:{3}\n*scale:{4}'.format(camera.GetPosition(), camera.GetFocalPoint(), camera.GetClippingRange(), camera.GetViewUp(), camera.GetParallelScale()))
    
def which_actor(actor):
    for name in all_actors:
        if actor == all_actors[name]:
            return name
    return 'unknown actor'
    
def print_all_actors(renderer):
    actors = renderer.GetActors()
    actors.InitTraversal()
    print('all actors in renderer:')
    for i in range(actors.GetNumberOfItems()):
        actor = actors.GetNextActor()
        print('Actor {}:\n{}\n'.format(which_actor(actor), actor))
        
def make_rings():
    global ring_transforms
    
    ring_transforms = []
    minR = km_to_au(saturn_rings['C'][0])
    maxR = km_to_au(saturn_rings['A'][1])
    for ring_name in ['C', 'B', 'A']:
        inner, outer = saturn_rings[ring_name]
        inner = km_to_au(inner)
        outer = km_to_au(outer)
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
        transform = vtk.vtkTransform()
        transform.Scale(scale*planet_scale, scale*planet_scale, scale*planet_scale)
        myfilter = vtk.vtkTransformPolyDataFilter()
        myfilter.SetTransform(transform)
        myfilter.SetInputData(aring)
        ring_transforms.append(myfilter)

# Callback function for the GUI slider
def slider_callback(obj, event):
    global time_step
    global refresh_rate
    value = obj.GetRepresentation().GetValue()
    time_step = value / refresh_rate  # Update the time_step
    
def planet_slider_callback(obj, event):
    global planet_sphere
    global render_window
    global sun_max
    planet_scale = obj.GetRepresentation().GetValue()
    # scale planets according to new value
    for name in bodies:
        body = bodies[name]
        if name == 'Sun':
            if planet_scale*body.radius < sun_max*bodies['Mercury'].orbit()['distance']:
                actual_scale = planet_scale
            else:
                actual_scale = sun_max*bodies['Mercury'].orbit()['distance']/body.radius
            planet_sphere[name].SetRadius(actual_scale*body.radius)
        else:
            planet_sphere[name].SetRadius(planet_scale*body.radius)
        planet_sphere[name].Update()
    # adjust Saturn rings accordingly
    for rt in ring_transforms:
        rt.GetTransform().Identity()
        rt.GetTransform().Scale(scale*planet_scale, scale*planet_scale, scale*planet_scale)
        rt.Update()
    render_window.Render()


# Handles the key press to change camera focal point
def key_pressed_callback(obj, event):
    global render_window
    global focus_planet
    global text_actor_planet
    global show_ecliptic
    global show_orbits
    global orbit_actors
    global ecliptic_actors
    global resolution
    global prev_resolution
    
    new_key = obj.GetKeySym()
    # print('key pressed:', new_key)
    if new_key.isdigit():
        if int(new_key) <= 9:
            focus_planet = body_names[int(new_key)]
            text_actor_planet.SetInput('Current Focus: {}'.format(focus_planet))

            # Overrides the default vtk behavior for keypress '3'
            if int(new_key) == 3:
                render_window.StereoRenderOn()
        elif int(new_key) == 9:
            focus_planet = 'N/A'
            text_actor_planet.SetInput('Current Focus: None')
    elif new_key == 'o' or new_key == 'O':
        show_orbits = not show_orbits
        if show_orbits:
            for name in orbit_actors:
                render_window.GetRenderers().GetFirstRenderer().AddActor(orbit_actors[name])
        else:
            for name in orbit_actors:
                render_window.GetRenderers().GetFirstRenderer().RemoveActor(orbit_actors[name])
        render_window.Render()
    elif new_key == 'f':
        if render_window.GetFullScreen():
            render_window.FullScreenOff()
        render_window.Render()
    elif new_key == 'h':
        print('List of keyboard commands:')
        print('\'0\', \'1\', \'2\', ..., \'8\': Point the camera to the Sun, Mercury, Venus, ..., Neptune and track')
        print('\'9\':                     Turn off tracking')
        print('\'o\' (\'O\'):               Turn on/off depiction of planets\' orbits')
        print('\'h\':                     Print this information')
    # else:
        # print('invalid key entry')

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

class Transient(object):
    def __init__(self, v=0, dv=0):
        self.v = v
        self.dv = dv 
        
    def __call__(self, day):
        return self.v + self.dv * day
 
class Orbit(object):
    def __init__(self, longitude_of_ascending_node, 
                 inclination, ascending_node_perihelion_angle, 
                 mean_distance, eccentricity, mean_anomaly):
        self.long = longitude_of_ascending_node
        self.incl = inclination
        self.peri = ascending_node_perihelion_angle
        self.dist = mean_distance
        self.ecc  = eccentricity
        self.ano  = mean_anomaly
        
    def __call__(self, day=0):
        return { 'longitude': self.long(day),
                 'inclination': self.incl(day),
                 'perihelion': self.peri(day),
                 'distance': self.dist(day),
                 'eccentricity': self.ecc(day),
                 'anomaly': self.ano(day) }
    
    def long_asc_node(self, day=0):
        return self.long(day)
        
    def inclination(self, day=0):
        return self.incl(day)
        
    def perihelion(self, day=0):
        return self.peri(day)
        
    def distance(self, day=0):
        return self.dist(day)
        
    def eccentricity(self, day=0):
        return self.ecc(day)
        
    def mean_anomaly(self, day=0):
        return self.ano(day)
     
class Body(object):
    def __init__(self, radius, orbit, day, tilt, year):
        self.orbit = orbit
        self.radius = radius
        self.day = day
        self.tilt = tilt
        self.year = year

def initialize_bodies():
    for i in range(9):
        bodies[body_names[i]] = Body(
            radius=scale*planet_radius[i], 
            orbit=Orbit(longitude_of_ascending_node = Transient(n_arr[i][0], n_arr[i][1]),
                        inclination = Transient(i_arr[i][0], i_arr[i][1]), 
                        ascending_node_perihelion_angle = Transient(w_arr[i][0], w_arr[i][1]), 
                        mean_distance = Transient(scale*a_arr[i][0], scale*a_arr[i][1]), 
                        eccentricity = Transient(e_arr[i][0], e_arr[i][1]),
                        mean_anomaly = Transient(m_arr[i][0], m_arr[i][1])),
            day=day_length[i], tilt=axial_tilt[i], year=years[i])
                                     
    # print('Bodies: ', bodies)
                                 
# Calculates Planet Rotations
def rotation_calculator(body_name, step):
    body = bodies[body_name]
    
    # Spin Resolution is based on time_step
    spin_degrees = (360 / body.day) * step

    # Venus and Uranus spin in retrograde
    if body_name=='Venus' or body_name=='Uranus':
        return spin_degrees  # positive for sun rising in the west
    else:
        return -spin_degrees  # negative for sun rising in the east


# Calculates planet positions using http://www.stjarnhimlen.se/comp/tutorial.html
def orbit_calculator(body_name, current_time, step):
    # Sun does not 'orbit'
    if body_name == 'Sun':
        return (0, 0, 0)
    
    body = bodies[body_name]

    # Covert to days since January 1, 2000 0:00
    j2000 = datetime(2000, 1, 1, 0, 0)
    delta = current_time - j2000
    d = delta.days + step
    
    orbital_info = body.orbit(d)
    anomaly = orbital_info['anomaly']
    perihelion = orbital_info['perihelion']
    eccentricity = orbital_info['eccentricity']
    distance = orbital_info['distance']
    inclination = orbital_info['inclination']
    longitude = orbital_info['longitude']

    # Step 3
    anomaly = (anomaly % 360)  # Modulus M so that it's between 0 and 360 degrees
    e_star = (180 / math.pi) * eccentricity
    big_e = anomaly + (e_star * math.sin(math.radians(anomaly)))

    # Do a third order approximation
    for order_idx in range(0, 3):
        delta_m = anomaly - (big_e - (e_star * math.sin(math.radians(big_e))))
        delta_e = delta_m / (1 - (eccentricity * math.cos(math.radians(big_e))))
        big_e += delta_e

    # Calculate rectangular coordinates
    x = distance * (math.cos(math.radians(big_e)) - eccentricity)
    y = distance * math.sqrt(1 - eccentricity * eccentricity) * math.sin(math.radians(big_e))

    r_orig = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
    r = r_orig
    v = math.atan2(y, x) * (180 / math.pi)

    # Compute ecliptic coordinates
    x_eclip = r * (math.cos(math.radians(longitude)) * math.cos(math.radians(v + perihelion)) -
                   math.sin(math.radians(longitude)) * math.sin(math.radians(v + perihelion)) * math.cos(math.radians(inclination)))
    y_eclip = r * (math.sin(math.radians(longitude)) * math.cos(math.radians(v + perihelion)) +
                   math.cos(math.radians(longitude)) * math.sin(math.radians(v + perihelion)) * math.cos(math.radians(inclination)))
    z_eclip = r * math.sin(math.radians(v + perihelion)) * math.sin(math.radians(inclination))

    return (x_eclip, y_eclip, z_eclip)


# Handles timers from the interactive render window and updates planet positions
class VtkTimerCallback:
    def __init__(self):
        self.timer_count = 0
        self.current_time = 0
        self.planet_actor = {}
        self.ring_actors = []
        self.arrow_actor = {}
        self.camera = vtk.vtkCamera()
        self.text = vtk.vtkTextActor()

    def execute(self, obj, event):
        global focus_planet
        global time_step
        global render_window
        
        for name in bodies:
            body = bodies[name]
            # Get new orbital position and rotation
            spin_degrees = rotation_calculator(name, time_step)
            x, y, z = orbit_calculator(name, self.current_time, self.timer_count)

            # Update position and rotation
            self.planet_actor[name].SetPosition(x, y, z)
            self.planet_actor[name].RotateZ(spin_degrees)  # Negative for prograde rotation, positive for retrograde            
            if name == 'Saturn':
                for ra in self.ring_actors: 
                    ra.SetPosition(x, y, z)
                    ra.RotateZ(spin_degrees)

            # Update camera for currently selected planet
            if name == focus_planet:
                self.camera.SetFocalPoint(x, y, z)
                self.text.SetInput('Current Focus: ' + focus_planet)
                
        obj.GetRenderWindow().Render()
        self.timer_count += time_step  # Increment time_step (in earth days)


def main():
    global render_window
    global focus_planet
    global refresh_rate
    global time_step
    global resize_image
    global planet_sphere
    global ring_transforms
    global text_actor_planet
    global orbit_actors
    global ecliptic_actors
    global show_ecliptic
    global show_orbits
    global resolution
    global planet_min_max

    # Timer Frequency
    refresh_rate = 30  # Hz
    timer_freq = 1000/refresh_rate  # ms
    time_step = 0.005  # in earth days

    # Render size set to HD Display
    focus_planet = 'Sun'  # Initial focus planet is the sun
    
    planet_pic = {}
    planet_reader = {}
    flip_planet_image = {}
    planet_texture = {}
    planet_sphere = {}
    texture_to_planet = {}
    planet_mapper = {}
    planet_actor = {}
    arrows_actor = {}
    ring_actors = []
    orbit_actors = {}
    ecliptic_actors = []
    
    initialize_bodies()
    
    # Create render window
    render_window = vtk.vtkRenderWindow()
    render_window.SetNumberOfLayers(2)
    if resolution[0] > 0 and resolution[1] > 0:
        render_window.SetSize(resolution[0], resolution[1])  # Set the window size
    else:
        # render_window.FullScreenOn()
        resolution = render_window.GetScreenSize()
        print('resolution: {}'.format(resolution))
        render_window.SetSize(resolution)
    render_window.SetWindowName("Solar System")
    render_window.StereoRenderOff()
    
    if show_arrows:
        arrow_source = vtk.vtkArrowSource()
        arrow_source.InvertOn()
        transform = vtk.vtkTransform()
        transform.RotateY(-90)
        transform.Scale(bodies['Sun'].radius, bodies['Sun'].radius, bodies['Sun'].radius)
        for name in bodies:
            body = bodies[name]
            p = body.orbit()['distance']
            transform_poly = vtk.vtkTransformPolyDataFilter()
            transform_poly.SetTransform(transform)
            transform_poly.SetInputConnection(arrow_source.GetOutputPort())
            arrow_mapper = vtk.vtkPolyDataMapper()
            arrow_mapper.SetInputConnection(transform_poly.GetOutputPort())
            arrows_actor[name] = vtk.vtkActor()
            arrows_actor[name].SetMapper(arrow_mapper)
            arrows_actor[name].SetPosition(p, 0, body.radius)
            arrows_actor[name].GetProperty().SetColor(1,0,0)
            arrows_actor[name].GetProperty().SetAmbient(1)
            all_actors['arrow(' + name + ')'] = arrows_actor[name]
        
    make_rings()

    idx = 0
    for name in bodies:
        body = bodies[name]
        # Read Planet Pictures
        planet_pic[name] = files[name]
        planet_reader[name] = vtk.vtkJPEGReader()
        planet_reader[name].SetFileName(planet_pic[name])

        # Flip Images
        flip_planet_image[name] = vtk.vtkImageFlip()
        flip_planet_image[name].SetFilteredAxis(0)  # flip x axis
        flip_planet_image[name].SetInputConnection(planet_reader[name].GetOutputPort())
        flip_planet_image[name].Update()

        # Create the Textures
        planet_texture[name] = vtk.vtkTexture()
        planet_texture[name].SetInputConnection(flip_planet_image[name].GetOutputPort())

        # Create the Planets
        planet_sphere[name] = vtk.vtkTexturedSphereSource()
        # planet_sphere[-1].SetCenter(0, 0, 0)
        if name == 'Sun':
            planet_sphere[name].SetRadius(planet_scale*body.radius)
        else:
            planet_sphere[name].SetRadius(planet_scale*body.radius)
        planet_sphere[name].SetThetaResolution(40)
        planet_sphere[name].SetPhiResolution(40)

        # Map Textures to Planets
        texture_to_planet[name] = vtk.vtkTextureMapToSphere()
        texture_to_planet[name].SetInputConnection(planet_sphere[name].GetOutputPort())
        texture_to_planet[name].PreventSeamOff()

        # Create Mappers
        planet_mapper[name] = vtk.vtkPolyDataMapper()
        planet_mapper[name].SetInputConnection(texture_to_planet[name].GetOutputPort())

        # Create Actors and add textures, rotations, and orbit position
        planet_actor[name] = vtk.vtkActor()
        planet_actor[name].SetMapper(planet_mapper[name])
        planet_actor[name].SetTexture(planet_texture[name])
        planet_actor[name].RotateY(body.tilt)
        planet_actor[name].RotateX(180)
        planet_actor[name].SetPosition(body.orbit.distance(), 0, 0)
        planet_actor[name].GetProperty().SetAmbient(ambient)
        planet_actor[name].GetProperty().SetSpecular(0.01)
        planet_actor[name].GetProperty().SetDiffuse(0.65)
        all_actors['planet ' + name] = planet_actor[name]
        
        j=0
        if name == 'Saturn':
            ring_reader = vtk.vtkJPEGReader()
            ring_reader.SetFileName(files['Rings'])
            ring_texture = vtk.vtkTexture()
            ring_texture.SetInputConnection(ring_reader.GetOutputPort())
            for rt in ring_transforms:
                mapper = vtk.vtkPolyDataMapper()
                mapper.SetInputConnection(rt.GetOutputPort())
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                actor.SetTexture(ring_texture)
                actor.RotateY(body.tilt)
                actor.RotateX(180)
                actor.SetPosition(body.orbit.distance(), 0, 0)
                actor.GetProperty().SetAmbient(ambient)
                ring_actors.append(actor)
                all_actors['ring #{}'.format(j)] = actor
                j = j+1

    # Create Lighting
    planet_light1 = vtk.vtkLight()
    planet_light1.SetPosition(0, 0, 0)
    planet_light1.SetFocalPoint(0.0000001*scale, 0, 0)  # Right hemisphere lighting
    planet_light1.SetConeAngle(180)
    planet_light1.SetPositional(1)

    planet_light2 = vtk.vtkLight()
    planet_light2.SetPosition(0, 0, 0)
    planet_light2.SetFocalPoint(-0.0000001*scale, 0, 0)  # Left hemisphere lighting
    planet_light2.SetConeAngle(180)
    planet_light2.SetPositional(1)

    # Add Background Starry Image
    jpeg_reader = vtk.vtkJPEGReader()
    jpeg_reader.SetFileName(files['Stars'])
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

    # Light up the sun no matter what
    sun_props = vtk.vtkProperty()
    sun_props.SetLighting(0)
    planet_actor['Sun'].SetProperty(sun_props)

    # Create a renderer and add the actors to the scene
    # print('Rendering Planets...')
    planet_renderer = vtk.vtkRenderer()
    planet_renderer.SetBackground(0, 0, 0)  # black background
    planet_renderer.AddLight(planet_light1)
    planet_renderer.AddLight(planet_light2)
    for name in bodies:
        planet_renderer.AddActor(planet_actor[name])
        # planet_renderer.AddActor(arrows_actor[name])
    for act in ring_actors:
        planet_renderer.AddActor(act)

    # Create the camera focused on the Sun
    cam1 = vtk.vtkCamera()
    cam1.SetFocalPoint(0, 0, 0)  # Sun Center
    cam1.SetPosition(-bodies['Jupiter'].orbit.distance(), 0, 0)
    cam1.Elevation(0)  # look at off-angle system (0 is top down)
    cam1.Azimuth(0)
    cam1.SetViewUp(0,0,1)

    planet_renderer.SetActiveCamera(cam1)
    planet_renderer.SetLayer(1)  # In front of background image

    # Add useful text to the renderer
    # Render Title
    text_actor_title = vtk.vtkTextActor()
    text_actor_title.SetInput("The Solar System")
    text_actor_title.SetPosition(40, resolution[1] - 90)
    text_actor_title.GetTextProperty().SetFontSize(60)
    text_actor_title.GetTextProperty().SetColor(1.0, 1.0, 1.0)
    all_actors['text title'] = text_actor_title

    # Current planet of focus
    text_actor_planet = vtk.vtkTextActor()
    text_actor_planet.SetInput('Current Focus: Sun')
    text_actor_planet.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    text_actor_planet.GetPositionCoordinate().SetValue(.85, .95)
    text_actor_planet.GetPosition2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    text_actor_planet.GetPosition2Coordinate().SetValue(.99, .95)
    # text_actor_planet.SetPosition(50, 40)
    text_actor_planet.GetTextProperty().SetFontSize(40)
    text_actor_planet.GetTextProperty().SetColor(1.0, 1.0, 1.0)
    all_actors['text planet'] = text_actor_planet

    # planet_renderer.AddActor2D(text_actor_title)
    # planet_renderer.AddActor2D(text_actor_subtitle)
    planet_renderer.AddActor2D(text_actor_planet)

    # Create a background renderer
    background_renderer = vtk.vtkRenderer()
    background_renderer.SetLayer(0)
    background_renderer.InteractiveOff()
    background_renderer.AddActor(background_actor)
    
    render_window.AddRenderer(planet_renderer)
    render_window.AddRenderer(background_renderer)
    
    # visualize ecliptic plane
    if show_ecliptic:
        ecliptic = vtk.vtkRegularPolygonSource()
        ecliptic.SetCenter(0,0,0)
        ecliptic.SetNumberOfSides(500)
        ecliptic.SetRadius(40*scale)
        ecmapper = vtk.vtkPolyDataMapper()
        ecmapper.SetInputConnection(ecliptic.GetOutputPort())
        e_actor = vtk.vtkActor()
        e_actor.SetMapper(ecmapper)
        e_actor.GetProperty().SetColor(1,1,1)
        e_actor.GetProperty().SetOpacity(0.2)
        e_actor.GetProperty().SetAmbient(1)
        all_actors['ecliptic plane'] = e_actor
        ecliptic_actors.append(e_actor)
        ecliptic_frame = vtk.vtkRegularPolygonSource()
        ecliptic_frame.SetCenter(0,0,0)
        ecliptic_frame.SetNumberOfSides(500)
        ecliptic_frame.SetRadius(40*scale)
        ecliptic_frame.GeneratePolygonOff()
        ecfmapper = vtk.vtkPolyDataMapper()
        ecfmapper.SetInputConnection(ecliptic_frame.GetOutputPort())
        ecf_actor = vtk.vtkActor()
        ecf_actor.SetMapper(ecfmapper)
        ecf_actor.GetProperty().SetColor(1,1,1)
        ecf_actor.GetProperty().SetOpacity(1)
        ecf_actor.GetProperty().SetAmbient(1)
        ecf_actor.GetProperty().SetLineWidth(5)
        ecf_actor.GetProperty().RenderLinesAsTubesOn()
        all_actors['ecliptic contour'] = ecf_actor
        ecliptic_actors.append(ecf_actor)
        planet_renderer.AddActor(e_actor)
        planet_renderer.AddActor(ecf_actor)
    
    # visualize projection of planets orbits on the ecliptic plane
    if show_orbits:
        for name in bodies:
            minz = 1000000
            maxz = -1000000
            if name == 'Sun':
                continue
            body = bodies[name]
            y = body.year
            now = datetime.now()
            
            points = vtk.vtkPoints()
            zvals = vtk.vtkFloatArray()
            zvals.SetNumberOfComponents(1)
            for d in range(int(y)):
                (x, y, z) = orbit_calculator(name, now, d)
                points.InsertNextPoint((x, y, z))
                if z<minz:
                    minz = z
                if z>maxz:
                    maxz = z
                zvals.InsertNextTuple([z])
                
            polyline = vtk.vtkPolyLine()
            polyline.GetPointIds().SetNumberOfIds(points.GetNumberOfPoints()+1)
            for i in range(points.GetNumberOfPoints()):
                polyline.GetPointIds().SetId(i,i)
            polyline.GetPointIds().SetId(points.GetNumberOfPoints(), 0)
            cells = vtk.vtkCellArray()
            cells.InsertNextCell(polyline)
            orbit = vtk.vtkPolyData()
            orbit.SetPoints(points)
            orbit.SetLines(cells)
            orbit.GetPointData().SetScalars(zvals)
            omapper = vtk.vtkPolyDataMapper()
            omapper.SetInputData(orbit)
            if name != 'Earth':
                ctf = vtk.vtkColorTransferFunction()
                ctf.AddRGBPoint(minz, 0, 0, 1)
                ctf.AddRGBPoint(0, 1, 1, 1)
                ctf.AddRGBPoint(maxz, 1, 0, 0) 
                omapper.SetLookupTable(ctf)
            else:
                omapper.ScalarVisibilityOff()
            oactor = vtk.vtkActor()
            oactor.SetMapper(omapper)
            oactor.GetProperty().SetColor(1,1,1)
            oactor.GetProperty().SetAmbient(1)
            oactor.GetProperty().RenderLinesAsTubesOn()
            orbit_actors[name] = oactor
            all_actors['orbit of ' + name] = oactor
            planet_renderer.AddActor(oactor)

    # Set-up interactor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    style = vtk.vtkInteractorStyleTrackballCamera()
    render_window_interactor.SetInteractorStyle(style)
    render_window_interactor.SetRenderWindow(render_window)
    render_window_interactor.Initialize()
    render_window.Render()
    
    render_window.AddObserver('ModifiedEvent', window_resized_callback)

    # Widget to control the animation speed (in earth-days per second)
    slider = vtk.vtkSliderRepresentation2D()
    slider.SetMinimumValue(0)
    slider.SetMaximumValue(30)
    slider.GetSliderProperty().SetLineWidth(1)
    slider.GetTubeProperty().SetLineWidth(1)
    slider.GetLabelProperty().SetFontSize(12)
    slider.GetLabelProperty().BoldOff()
    slider.SetValue(time_step * timer_freq)
    slider.SetTitleText("Days per Second")
    slider.SetLabelFormat("%.2f")  # float format
    slider.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider.GetPoint1Coordinate().SetValue(.02, .1)
    slider.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider.GetPoint2Coordinate().SetValue(.2, .1)

    day_widget = vtk.vtkSliderWidget()
    day_widget.SetInteractor(render_window_interactor)
    day_widget.SetRepresentation(slider)
    day_widget.SetAnimationModeToAnimate()
    day_widget.EnabledOn()
    day_widget.AddObserver("InteractionEvent", slider_callback)

    # Widget to control the scale of the planets
    slider2 = vtk.vtkSliderRepresentation2D()
    slider2.SetMinimumValue(planet_min_max[0])
    slider2.SetMaximumValue(planet_min_max[1])
    slider2.GetLabelProperty().SetFontSize(12)
    slider2.GetSliderProperty().SetLineWidth(2)
    slider2.GetLabelProperty().BoldOff()
    slider2.SetValue(0.1*sun_max)
    slider2.SetTitleText("Bodies scaling factor")
    slider2.SetLabelFormat("%.2f")  # integer format
    slider2.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider2.GetPoint1Coordinate().SetValue(.30, .1)
    slider2.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider2.GetPoint2Coordinate().SetValue(.98, .1)

    planet_widget = vtk.vtkSliderWidget()
    planet_widget.SetInteractor(render_window_interactor)
    planet_widget.SetRepresentation(slider2)
    planet_widget.SetAnimationModeToAnimate()
    planet_widget.EnabledOn()
    planet_widget.AddObserver("InteractionEvent", planet_slider_callback)

    # Sign up to receive TimerEvent
    cb = VtkTimerCallback()
    cb.planet_actor = planet_actor
    cb.arrow_actor = arrows_actor
    cb.ring_actors = ring_actors
    cb.camera = cam1
    cb.text = text_actor_planet
    cb.current_time = datetime.now()
    render_window_interactor.AddObserver('TimerEvent', cb.execute)
    render_window_interactor.AddObserver('KeyPressEvent', key_pressed_callback)
    render_window_interactor.CreateRepeatingTimer(int(timer_freq))  # ms between calls
    render_window_interactor.EnableRenderOn()
    render_window_interactor.Initialize()
    
    # planet_renderer.ResetCamera()
    # planet_renderer.GetActiveCamera().SetClippingRange([0.000001, 100])
    render_window.Render()
    
    
    render_window_interactor.Start()
    


if __name__ == '__main__':
    global show_ecliptic
    global show_orbits
    global planet_min_max
    global resolution
    parser = argparse.ArgumentParser(
            description='Visualize Solar System (based on code by Drew Brandsen, Purdue CS530, 2018)',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-e', '--ecliptic', action='store_true', help='Show ecliptic plane')
    parser.add_argument('-o', '--orbits', action='store_true', help='Show planets\' orbits')
    parser.add_argument('-r', '--resolution', metavar='int', default=[2560, 1440], type=int, nargs=2, help='window resolution (-1: full screen)')
    parser.add_argument('-s', '--sun', metavar='float', default='10', type=float, help='Max sun radius (in multiples of Mercury semi-major axis)')
    parser.add_argument('-p', '--planet', metavar='float', default=[1, 2000], type=float, nargs=2, help='Range of scaling factors for planets')
    parser.add_argument('-x', '--scale', metavar='int', default=10000, type=int, help='Uniform space scaling factor')
    args = parser.parse_args()
    
    print('args: {}'.format(args))
    
    show_ecliptic = args.ecliptic
    show_orbits = args.orbits
    scale = args.scale
    sun_max = args.sun
    planet_min_max = args.planet
    resolution = args.resolution
        
    main()
