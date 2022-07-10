#!/usr/bin/env python3
import matplotlib.pyplot as plt

from ztf_plan_obs.plan import PlanObservation
from ztf_plan_obs.multiday_plan import MultiDayObservation
from ztf_plan_obs.api import Queue

name = "IC220624A"  # Name of the alert object
# name = "IC201007A"
date = "2022-06-24"  # This is optional, defaults to today
# ra = 242.58
# dec = 11.61
# Now no ra and dec values are given, but alertsource is set to 'icecube'. This enables GCN archive parsing for the alert name. If it is not found, it will use the latest GCN notice (these are automated).

plan = PlanObservation(name=name, date=None, alertsource="icecube")
plan.plot_target()  # Plots the observing conditions
plan.request_ztf_fields()  # Checks in which ZTF fields this object is observable
# plan.plot_finding_chart()
plt.close()


observationplan = MultiDayObservation(name=name, startdate="2022-06-23")
observationplan.print_plan()
summary = observationplan.summarytext

triggers = observationplan.triggers

q = Queue(user="DESY")

for i, trigger in enumerate(triggers):
    q.add_trigger_to_queue(
        trigger_name=f"ToO_{name}",
        validity_window_start_mjd=trigger["mjd_start"],
        field_id=trigger["field_id"],
        filter_id=trigger["filter_id"],
        exposure_time=trigger["exposure_time"],
    )

q.print()
