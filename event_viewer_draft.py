def date_change(self, et):

    current = self.state.clock
    # self.state.clock = Units.QDateTime2time(datetime)
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


def vantage_point(camera_name: str, data: DataLoader, be=False) -> vtk.vtkCamera:
        
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
            if be:
                camera.SetPosition(0, 0, 5*r)
                camera.SetViewUp(0, 1, 0)
            else:
                camera.SetPosition(5*r, 0, 0)
                camera.SetViewUp(0, 0, 1)
            camera.SetFocalPoint(0, 0, 0)
            camera.SetViewAngle(45)
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


def change_planet_focus(target):
    if target != self.state.params.planet_focus:
        myprint(verbose=True, text='changing focus from {} to {}'.format(
            self.state.params.planet_focus, target))
        self.state.params.planet_focus = target
        self.state.planet_focus_changed = True
        self.state.params.do_tether = True
        cam_setting = VTKUtils.vantage_point(
            target_name=self.state.params.planet_focus, data=self.data, et=self.state.cloc, be=True) # TODO: what is target_name, et?
        VTKUtils.get_camera(self.graphics.window).SetPosition(
            cam_setting['pos'])
        VTKUtils.get_camera(self.graphics.window).SetFocalPoint(
            cam_setting['focal'])
        VTKUtils.get_camera(self.graphics.window).SetViewUp(
            cam_setting['up'])
        self.graphics.all_actors['text']['bodies'].SetInput(
            'Current Focus: {}'.format(self.state.params.planet_focus))


def play_event(self, event_name):

    events_time_st = {
        'launch': 'year month, day hour:minute:second'
    }
    events_time_ed = {
        'launch': 'year month, day hour:minute:second'
    }
    events_time_scale = {
        'launch': time_step_changed_1minute_cb
    }
    events_target = {
        'launch': (None, 'Earth')
    }

    time_st = events_time_st[event_name]
    et_st = spice.utc2et(time_st)
    date_change(et_st)

    events_time_scale[event_name]()

    target = events_target[event_name]
    change_planet_focus(target)

    self.time_end = events_time_ed[event_name]
    # TODO: need to poll for when self.state.clock >= self.time_end
