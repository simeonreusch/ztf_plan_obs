# ztf_plan_obs
Toolset for planning observations with ZTF. GCN parsing is currently only implemented for IceCube alerts.

It checks if the object is observable with a maximum airmass on a given date, plots the airmass vs. time, computes two optimal (minimal airmass at night) observations of 300s in g- and r and generate the ZTF field plots for all fields having a reference.

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
plan.request_ztf_fields() # Checks in which ZTF fields this 
# object is observable and generate plot for them.
```
The observation plot and the ZTF field plots will be located in the current directory/[name]
![](examples/figures/observation_plot_generic.png)

Note: Checking if fields have references requires ztfquery, which needs IPAC credentials.

# Usage for IceCube alerts
```python
from ztf_plan_obs.plan import PlanObservation

name = "IC201007A" # Name of the alert object
date = "2020-10-08" #This is optional, defaults to today
# Now no ra and dec values are given, but alertsource 
# is set to 'icecube'. This enables GCN archive parsing 
# for the alert name. If it is not found, it will use 
#the latest GCN notice (these are automated).

plan = ObservationPlan(name=name, date=date, alertsource="icecube")
plan.plot_target() # Plots the observing conditions
plan.request_ztf_fields() # Checks in which ZTF fields 
ä# this object is observable and generate plot for them.
```
![](examples/figures/observation_plot_icecube.png)

# Triggering ZTF

`ztf_plan_obs` can be used for directly scheduling ToO observations with ZTF. 
This is done through API calls to the `Kowalski` system, managed by the kowalski python manager [penquins](https://github.com/dmitryduev/penquins).

To use this functionality, you must first configure the connection details. You need both an API token, and to know the address of the Kowalski host address.
You can then set these as environment variables:

```bash
export KOWALSKI_HOST=something
export KOWALSKI_API_TOKEN=somethingelse
```

You can then import helper functions for querying, submitting and deleting ToO schedules:

## Querying

```python
from ztf_plan_obs.api import get_too_queues
existing_too_requests = get_too_queues()
print([x["name"] for x in existing_too_requests["data"]])
```

## Submitting

```python
from ztf_plan_obs.api import build_request, submit_request, get_too_queues

queue_name = "ToO_IC220513A_test"

request = build_request(
    user="yourname",
    queue_name=queue_name,
    validity_window_start_mjd=59719.309333333334,
    validity_window_end_mjd=59719.30988055556,
    subprogram_name="ToO_Neutrino",
    field_ids=[427],
    filter_ids=[1],
    exposure_times=[600.]
) 
print(request)
submit_request(request)

queue = get_too_queues()

entries = [x["name"] for x in queue["data"]]
print(entries)
assert queue_name in entries
```

## Deleting
```python
from ztf_plan_obs.api import delete_request

queue_name = "ToO_IC220513A_test"

res = delete_request(user="yourname", queue_name=queue_name)
print(res)
```