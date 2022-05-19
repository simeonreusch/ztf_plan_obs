#!/usr/bin/env python3
from ztf_plan_obs.plan import PlanObservation
from ztf_plan_obs.multiday_plan import MultiDayObservation
import matplotlib.pyplot as plt

name = "IC220513A"  # Name of the alert object
# name = "IC201007A"
date = "2022-05-19"  # This is optional, defaults to today
# ra = 242.58
# dec = 11.61
# Now no ra and dec values are given, but alertsource is set to 'icecube'. This enables GCN archive parsing for the alert name. If it is not found, it will use the latest GCN notice (these are automated).

plan = PlanObservation(name=name, date=date, alertsource="icecube")
plan.plot_target()  # Plots the observing conditions
plan.request_ztf_fields()  # Checks in which ZTF fields this object is observable
# # plan.plot_finding_chart()
plt.close()

name = "IC220513A"
date = "2022-05-19"

observationplan = MultiDayObservation(name=name, startdate=date)
observationplan.print_plan()
summary = observationplan.summarytext
