# ztf_plan_obs
Toolset for planning observations with ZTF. GCN parsing is currently only implemented for IceCube alerts.

It checks if the object is observable with a maximum airmass on a given date, plots the airmass vs. time, computes two optimal (minimal airmass at night) observations of 300s in g- and r and downloads the ZTF field plots from [here](http://yupana.caltech.edu/cgi-bin/ptf/tb/zoc?begin) for all fields having a reference.

The output is designed so it can be directly pasted into the [GROWTH ToO Marshal](http://skipper.caltech.edu:8081/login?next=%2Fplan_manual) to schedule an observation.

# Requirements
[ztfquery](https://github.com/mickaelrigault/ztfquery) for checking if fields have a reference.

ztf_plan_obs requires at least Python 3.7

# Installation
Using Pip: ```pip3 install ztf_plan_obs```.

Otherwise, you can clone the repository: ```git clone https://github.com/simeonreusch/ztf_plan_obs```. This also gives you access to the Slackbot.

# General usage
```python
from ztf_plan_obs.plan import PlanObservation

name = "testalert" # Name of the alert object
date = "2020-05-05" #This is optional, defaults to today
ra = 133.7
dec = 13.37

plan = ObservationPlan(name=name, date=date, ra=ra, dec=dec)
plan.plot_target() # Plots the observing conditions
plan.request_ztf_fields() # Checks in which ZTF fields this object is observable and download plot for them from http://yupana.caltech.edu
```
The observation plot and the ZTF field plots will be located in the current directory/[name]
![](examples/figures/observation_plot_generic.png)

Note: Checking if fields have references requires ztfquery, which needs IRSA credentials.

# Usage for IceCube alerts
```python
from ztf_plan_obs.plan import PlanObservation

name = "IC201007A" # Name of the alert object
date = "2020-10-08" #This is optional, defaults to today
# Now no ra and dec values are given, but alertsource is set to 'icecube'. This enables GCN archive parsing for the alert name. If it is not found, it will use the latest GCN notice (these are automated).

plan = ObservationPlan(name=name, date=date, alertsource="icecube")
plan.plot_target() # Plots the observing conditions
plan.request_ztf_fields() # Checks in which ZTF fields this object is observable and download plot for them from http://yupana.caltech.edu
```
![](examples/figures/observation_plot_icecube.png)
