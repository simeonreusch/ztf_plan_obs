#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import time, os
import astropy
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
import os, time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astroplan as ap
from astroplan import Observer
from datetime import datetime
from astroplan.plots import plot_airmass, plot_altitude
import requests
from bs4 import BeautifulSoup


class ObservationPlan:
    def __init__(
        self,
        ra: float,
        dec: float,
        name: str,
        arrivaltime: str = None,
        **kwargs,
    ):

        self.ra = ra
        self.dec = dec
        self.name = name
        self.arrivaltime = arrivaltime

        self.coordinates = SkyCoord(self.ra * u.deg, self.dec * u.deg, frame="icrs")

        self.target = ap.FixedTarget(name=self.name, coord=self.coordinates)

        self.palomar = Observer.at_site("Palomar", timezone="US/Pacific")

        self.now = Time(datetime.utcnow())

        constraints = [
            ap.AltitudeConstraint(20 * u.deg, 90 * u.deg),
            ap.AirmassConstraint(5),
            ap.AtNightConstraint.twilight_astronomical(),
        ]

    def plot_target(self):
        ax = plot_airmass(
            self.target,
            self.palomar,
            self.now,
            brightness_shading=True,
            altitude_yaxis=False,
            max_airmass=2.5,
        )
        ax.axvline(Time(self.now).plot_date, color="black", label="now", ls="dotted")
        if self.arrivaltime is not None:
            ax.axvline(
                Time(self.arrivaltime).plot_date,
                color="red",
                label="neutrino arrival",
                ls="dotted",
            )
        plt.tight_layout()
        plt.savefig(f"{self.name}_airmass.png")

    def request_ztf_fields(self):
        URL = "http://yupana.caltech.edu/cgi-bin/ptf/tb//zoc"
        IMAGE_URL_1 = "http://yupana.caltech.edu/marshals/tb//igmo_0_0.png"
        IMAGE_URL_2 = "http://yupana.caltech.edu/marshals/tb//igmo_0_1.png"
        IMAGE_URLS = [IMAGE_URL_1, IMAGE_URL_2]

        objra = self.ra
        objdec = self.dec
        radius = 60

        fieldids_total = []

        for grid in [1, 2]:

            request_data = {
                "showobject": 1,
                "objra": objra,
                "objdec": objdec,
                "grid": grid,
                "objname": "unknown",
                "radam": radius,
                "submitshowobject": "SUBMIT (Show Object)",
            }

            # Post the request
            # response = requests.post(url=URL, data=request_data, timeout=30)
            response = requests.get(URL, params=request_data)
            img_data = requests.get(IMAGE_URL_1).content
            with open("test.png", "wb") as handler:
                handler.write(img_data)

            soup = BeautifulSoup(response.text, "html5lib")

            # try:
            pre = soup.find_all("pre")[-1]
            results = pre.text.split("\n")[1:-3]
            fieldids = []
            for result in results:
                if len(result) > 10:
                    fieldid = int(result[0:7])
                    fieldids.append(fieldid)
                    fieldids_total.append(fieldid)

            for index, fieldid in enumerate(fieldids):
                img_data = requests.get(IMAGE_URLS[index]).content
                with open(f"grid_{fieldid}.png", "wb") as handler:
                    handler.write(img_data)

        print(f"Fields that are possible: {fieldids_total}")


# NEED TO INCLUDE ERROR CIRCLE CALCULATION
# RA = 96.46
# DEC = -4.33
RA = 40
DEC = 45
NAME = "IC200926A"
ARRIVALTIME = "2020-09-26 07:54:11.621"

plan = ObservationPlan(ra=RA, dec=DEC, name=NAME, arrivaltime=ARRIVALTIME)

# plan.plot_target()
plan.request_ztf_fields()
