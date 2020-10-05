# ztf_plan_too
Toolset for planning target of opportunity observations with ZTF

# Installation
Using Pip: ```pip3 install git+https://github.com/simeonreusch/ztf_plan_too```

Otherwise, you can clone the repository: ```git clone https://github.com/simeonreusch/ztf_plan_too```

# Usage
```python
from ztf_plan_too.plan import ObservationPlan

NAME = "IC200929A" # Name of the alert object
date = "2020-10-05" #This is optional, defaults to today

plan = ObservationPlan(name=NAME, date=date)
plan.plot_target() # Plots the observing conditions
plan.request_ztf_fields() # Checks in which ZTF fields this object is observable
```

Note: Checking if fields have references requires ztfquery, which needs IRSA credentials
