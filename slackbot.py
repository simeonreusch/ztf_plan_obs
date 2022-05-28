#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

from ztf_plan_obs.plan import PlanObservation, AirmassError, ParsingError
from ztf_plan_obs.utils import is_ztf_name
from ztf_plan_obs.multiday_plan import MultiDayObservation
from ztf_plan_obs.api import Queue
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
        submit_trigger=False,
        alertsource=None,
        site=None,
    ):
        self.channel = channel
        self.name = name
        self.ra = ra
        self.dec = dec
        self.date = date
        self.multiday = multiday
        self.submit_trigger = submit_trigger
        self.alertsource = alertsource
        self.site = site
        self.fields = None
        self.recommended_field = None
        self.coverage = None

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
                site=self.site,
            )
            plan.plot_target()
            plt.close()

            self.summary = plan.get_summary()
            if plan.observable is True:
                if self.site == "Palomar":
                    self.fields = plan.request_ztf_fields()
                    if plan.ra_err:
                        self.recommended_field = plan.recommended_field
                        self.coverage = plan.coverage
                else:
                    self.summary = "No fields available (select 'Palomar' as site)"
            else:
                self.summary = "Not observable!"

            if self.multiday:
                multiday_plan = MultiDayObservation(
                    name=self.name, ra=self.ra, dec=self.dec, startdate=self.date
                )
                if plan.observable is True:
                    self.multiday_summary = multiday_plan.summarytext
                else:
                    self.multiday_summary = "Not observable!"

                if self.submit_trigger:
                    triggers = multiday_plan.triggers
                    q = Queue(user="DESY")

                    for i, trigger in enumerate(triggers):
                        q.add_trigger_to_queue(
                            trigger_name=f"ToO_{self.name}",
                            validity_window_start_mjd=trigger["mjd_start"],
                            field_id=trigger["field_id"],
                            filter_id=trigger["filter_id"],
                            exposure_time=trigger["exposure_time"],
                        )
                    q.submit_queue()

                    self.multiday_summary += f"\nYOU HAVE TRIGGERED ALL OBSERVATIONS ({len(q.queue)} in total)!\nCheck with 'Queue -get' if they have been added successfully.\nYour triggers:\n"

                    triggertext = multiday_plan.print_triggers()
                    self.multiday_summary += triggertext

        except ParsingError:
            self.summary = "GCN parsing error"
