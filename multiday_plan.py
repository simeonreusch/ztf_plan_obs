#!/usr/bin/env python3
import matplotlib.pyplot as plt
from datetime import date, datetime
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
from astropy.time import Time
from astropy import units as u
from ztf_plan_obs.plan import PlanObservation
from ztf_plan_obs.plan import time_shortener as ts

NIGHTS = [1, 2, 3, 5, 7, 9]
SHORT_NIGHTS = NIGHTS[1:]

name = "testalert"
today = date.today()
now = datetime.now()
now_astropy = Time(str(now), format="iso", scale="utc", out_subfmt="date")
next_days = [(now_astropy + i - 1).value for i in NIGHTS]

ra = 180
dec = 50

plan_initial = PlanObservation(
    name=name, date=str(today), ra=ra, dec=dec
)  # alertsource="icecube")
ra = plan_initial.ra
dec = plan_initial.dec


observable = []
g_band_start = []
g_band_end = []
r_band_start = []
r_band_end = []

with PdfPages(f"{name}_multiday.pdf") as pdf:
    for i, day in enumerate(tqdm(next_days)):
        if i + 1 not in SHORT_NIGHTS:
            plan = PlanObservation(name=name, date=day, ra=ra, dec=dec, verbose=False)
        else:
            plan = PlanObservation(
                name=name, date=day, ra=ra, dec=dec, observationlength=30, verbose=False
            )
        observable.append(plan.observable)
        if observable:
            g_band_start.append(ts(plan.g_band_recommended_time_start))
            g_band_end.append(ts(plan.g_band_recommended_time_end))
            r_band_start.append(ts(plan.r_band_recommended_time_start))
            r_band_end.append(ts(plan.r_band_recommended_time_end))
        else:
            g_band_start.append(None)
            g_band_end.append(None)
            r_band_start.append(None)
            r_band_end.append(None)
        ax = plan.plot_target()
        plt.tight_layout()
        pdf.savefig()
        plt.close()

print(g_band_start)
print(g_band_end)
