#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

from ztf_plan_too.plan import ObservationPlan


class ObsBot:

    # The constructor for the class. It takes the channel name as the a
    # parameter and then sets it as an instance variable
    def __init__(self, channel, name):
        self.channel = channel
        self.name = name

    # Craft and return the entire message payload as a dictionary.
    def create_plot(self):
        plan = ObservationPlan(name=self.name)
        return
