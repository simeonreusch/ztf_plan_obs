#!/usr/bin/env python3
from ztf_plan_obs.plan import PlanObservation

# name = "testalert"  # Name of the alert object
name = "ZTF20actqnhg"
date = "2020-11-29"  # This is optional, defaults to today
# ra = 133.7
# dec = 13.37

plan = PlanObservation(name=name, date=date)
plan.plot_target()  # Plots the observing conditions
plan.request_ztf_fields()  # Checks in which ZTF fields this object is observable
