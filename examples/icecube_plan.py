#!/usr/bin/env python3
from ztf_plan_obs.plan import PlanObservation

name = "IC201007A"  # Name of the alert object
date = "2020-10-08"  # This is optional, defaults to today
# Now no ra and dec values are given, but alertsource is set to 'icecube'. This enables GCN archive parsing for the alert name. If it is not found, it will use the latest GCN notice (these are automated).

plan = PlanObservation(name=name, date=date, alertsource="icecube")
plan.plot_target()  # Plots the observing conditions
plan.request_ztf_fields()  # Checks in which ZTF fields this object is observable
