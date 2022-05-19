#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

from ztf_plan_obs.plan import PlanObservation, AirmassError, ParsingError, is_ztf_name
from ztf_plan_obs.multiday_plan import MultiDayObservation
import matplotlib.pyplot as plt


class Slackbot:
    def __init__(
        self,
        channel,
        name,
        ra=None,
        dec=None,
        date=None,
        multiday=False,
        alertsource=None,
    ):
        self.channel = channel
        self.name = name
        self.ra = ra
        self.dec = dec
        self.date = date
        self.multiday = multiday
        self.alertsource = alertsource

    # Craft and return the entire message payload as a dictionary.
    def create_plot(self):
        try:
            plan = PlanObservation(
                name=self.name,
                ra=self.ra,
                dec=self.dec,
                date=self.date,
                multiday=self.multiday,
                alertsource=self.alertsource,
            )
            plan.plot_target()
            plt.close()
            self.summary = plan.get_summary()
            if plan.observable is True:
                self.fields = plan.request_ztf_fields()
            else:
                self.summary = "Not observable!"
                self.fields = None

            if self.multiday:
                multiday_plan = MultiDayObservation(
                    name=self.name, ra=self.ra, dec=self.dec, startdate=self.date
                )
                if plan.observable is True:
                    self.multiday_summary = multiday_plan.summarytext
                else:
                    self.multiday_summary = "Not observable!"

        except ParsingError:
            self.summary = "GCN parsing error"
            self.fields = None
