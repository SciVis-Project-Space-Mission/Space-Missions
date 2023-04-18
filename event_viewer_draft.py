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


# TODO: number keys are bugged, use change_planet_focus() instead of vantage_point()
def change_planet_focus_key(self, target):
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


def play_event_cb(self, event_name):

    events_time_st = {
        'launch': 'year month, day hour:minute:second'
    }
    events_time_ed = {
        'launch': 'year month, day hour:minute:second'
    }
    events_time_scale = {
        'launch': self.time_step_changed_1minute_cb
    }
    events_frame = {
        'launch': (None, 'Earth')
    }

    time_st = events_time_st[event_name]
    et_st = spice.utc2et(time_st)
    date_change(et_st)

    events_time_scale[event_name]()

    frame = events_frame[event_name]
    anchor = frame[0]
    target = frame[1]
    event_camera = VTKUtils.vantage_point(camera_name=frame, data=self.data, be=True)
    self.change_planet_focus(anchor, target, ref_camera=event_camera)

    self.time_end = events_time_ed[event_name]
    # TODO: need to poll for when self.state.clock >= self.time_end
