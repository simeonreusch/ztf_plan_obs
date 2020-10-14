# ztf_plan_obs
Toolset for planning observations with ZTF. GCN parsing is currently only implemented for IceCube alerts.

# Installation
Using Pip: ```pip3 install ztf_plan_obs```

Otherwise, you can clone the repository: ```git clone https://github.com/simeonreusch/ztf_plan_obs```

# General usage
```python
from ztf_plan_obs.plan import ObservationPlan

name = "testalert" # Name of the alert object
date = "2020-05-05" #This is optional, defaults to today
ra = 133.7
dec = 13.37

plan = ObservationPlan(name=name, date=date, ra=ra, dec=dec)
plan.plot_target() # Plots the observing conditions
plan.request_ztf_fields() # Checks in which ZTF fields this object is observable
```

Note: Checking if fields have references requires ztfquery, which needs IRSA credentials

# Usage for IceCube alerts
```python
from ztf_plan_obs.plan import ObservationPlan

name = "IC201007A" # Name of the alert object
date = "2020-10-08" #This is optional, defaults to today
# Now no ra and dec values are given, but alertsource is set to 'icecube'. This enables GCN archive parsing for the alert name. If it is not found, it will use the latest GCN notice (these are automated).

plan = ObservationPlan(name=name, date=date, alertsource="icecube")
plan.plot_target() # Plots the observing conditions
plan.request_ztf_fields() # Checks in which ZTF fields this object is observable
```
