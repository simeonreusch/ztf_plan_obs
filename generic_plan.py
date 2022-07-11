#!/usr/bin/env python3
from ztf_plan_obs.plan import PlanObservation

name = "ZTF19accdntg"  # Name of the alert object
date = "2021-07-22"  # This is optional, defaults to today
# ra = 133.7
# dec = 13.37

plan = PlanObservation(name=name, date=date, alertsource="ZTF")
plan.plot_target()  # Plots the observing conditions
plan.request_ztf_fields()  # Checks in which ZTF fields this object is observable
# plan.plot_finding_chart()
print(plan.recommended_field)
