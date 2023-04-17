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


def play_event(self, event_name):

    events_time_st = {
        'launch': 'year month, day hour:minute:second'
    }
    events_time_ed = {
        'launch': 'year month, day hour:minute:second'
    }
    events_ref = {
        'launch': ('Earth', None)
    }

    # TODO: start paused?

    time_st = events_time_st[event_name]
    et_st = spice.utc2et(time_st)
    date_change(et_st)

    # TODO: orient camera birds eye for event_ref

    self.time_end = events_time_ed[event_name]
    # TODO: need to poll for when self.state.clock >= self.time_end
