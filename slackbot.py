#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

from ztf_plan_too.plan import ObservationPlan


class ObsBot:

    # The constructor for the class. It takes the channel name as the a
    # parameter and then sets it as an instance variable
    def __init__(self, channel, name, ra=None, dec=None, date=None):
        self.channel = channel
        self.name = name
        self.ra = ra
        self.dec = dec
        self.date = date

    # Craft and return the entire message payload as a dictionary.
    def create_plot(self):
        plan = ObservationPlan(name=self.name, ra=self.ra, dec=self.dec, date=self.date)
        plan.plot_target()
        self.summary = plan.get_summary()
