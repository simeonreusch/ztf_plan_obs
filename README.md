# ztf_plan_obs
Toolset for planning observations with ZTF. Currently only parsing for IceCube neutrino alerts is implemented.

# Installation
Using Pip: ```pip3 install ztf_plan_obs```

Otherwise, you can clone the repository: ```git clone https://github.com/simeonreusch/ztf_plan_obs```

# Usage
```python
from ztf_plan_obs.plan import ObservationPlan

NAME = "IC200929A" # Name of the alert object
date = "2020-10-05" #This is optional, defaults to today
# you can also pass ra and dec values. If no values are given, it checks
# the GCN archive

plan = ObservationPlan(name=NAME, date=date)
plan.plot_target() # Plots the observing conditions
plan.request_ztf_fields() # Checks in which ZTF fields this object is observable
```

Note: Checking if fields have references requires ztfquery, which needs IRSA credentials
