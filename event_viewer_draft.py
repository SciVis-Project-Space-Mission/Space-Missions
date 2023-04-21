# TODO: number keys are bugged, use change_planet_focus() instead of vantage_point()
def change_planet_focus_key(self, target):
    if target != self.state.params.planet_focus:
        myprint(verbose=True, text='changing focus from {} to {}'.format(
            self.state.params.planet_focus, target))
        self.state.params.planet_focus = target
        self.state.planet_focus_changed = True
        self.state.params.do_tether = True
        cam_setting = VTKUtils.vantage_point(
            target_name=self.state.params.planet_focus, data=self.data, et=self.state.cloc, be=True) # params are wrong
        VTKUtils.get_camera(self.graphics.window).SetPosition(
            cam_setting['pos'])
        VTKUtils.get_camera(self.graphics.window).SetFocalPoint(
            cam_setting['focal'])
        VTKUtils.get_camera(self.graphics.window).SetViewUp(
            cam_setting['up'])
        self.graphics.all_actors['text']['bodies'].SetInput(
            'Current Focus: {}'.format(self.state.params.planet_focus))
