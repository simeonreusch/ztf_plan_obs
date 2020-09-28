#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import time, os, warnings
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
from ztfquery import fields


class ObservationPlan:
    """ """

    def __init__(
        self,
        ra: float,
        dec: float,
        name: str,
        arrivaltime: str = None,
        date: str = None,
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
        self.date = date

        constraints = [
            ap.AltitudeConstraint(20 * u.deg, 90 * u.deg),
            ap.AirmassConstraint(5),
            ap.AtNightConstraint.twilight_astronomical(),
        ]

        if not os.path.exists(self.name):
            os.makedirs(self.name)

    def plot_target(self):
        """ """
        if self.date is not None:
            _date = self.date + " 12:00:00.000000"
            time = _date
        else:
            time = self.now

        ax = plot_airmass(
            self.target,
            self.palomar,
            time,
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

        if self.date is not None:
            ax.set_xlabel(f"{self.date} [UTC]")
        else:
            ax.set_xlabel(f"{self.now.datetime.date()} [UTC]")
        outpath = os.path.join(self.name, f"{self.name}_airmass.png")
        plt.grid(True, color="green", linestyle="dotted", which="both")
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(int(start), int(end), 1 / 12))
        plt.savefig(outpath)

        # NOTE: INCLUDE MOON AND SUN IN THIS!

    def request_ztf_fields(self):
        """
        This looks at yupana.caltech.edu for the fields matching
        your location and downloads the camera grid plots for these
        """
        URL = "http://yupana.caltech.edu/cgi-bin/ptf/tb//zoc"
        IMAGE_URL_1 = "http://yupana.caltech.edu/marshals/tb//igmo_0_0.png"
        IMAGE_URL_2 = "http://yupana.caltech.edu/marshals/tb//igmo_0_1.png"
        IMAGE_URL_3 = "http://yupana.caltech.edu/marshals/tb//igmo_0_2.png"
        IMAGE_URL_4 = "http://yupana.caltech.edu/marshals/tb//igmo_0_3.png"
        IMAGE_URLS = [IMAGE_URL_1, IMAGE_URL_2, IMAGE_URL_3, IMAGE_URL_4]

        objra = self.ra
        objdec = self.dec
        radius = 60

        fieldids_total = []
        fieldids_total_ref = []

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

            # Get information on fields from yupana.caltech.edu
            response = requests.get(URL, params=request_data)

            # Parse the HTML response
            soup = BeautifulSoup(response.text, "html5lib")

            pre = soup.find_all("pre")[-1]
            results = pre.text.split("\n")[1:-3]
            fieldids = []
            fieldids_ref = []

            for result in results:
                if len(result) > 10:
                    fieldid = int(result[0:7])
                    fieldids.append(fieldid)
                    fieldids_total.append(fieldid)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        refs = fields.has_field_reference(fieldid)
                    if refs["zg"] == True and refs["zr"] == True:
                        fieldids_ref.append(fieldid)
                        fieldids_total_ref.append(fieldid)

            # Download the camera images (URLS are static, images
            # seem to be regenerated after each request)
            for index, fieldid in enumerate(fieldids_ref):
                img_data = requests.get(IMAGE_URLS[index]).content
                outpath = os.path.join(self.name, f"{self.name}_grid_{fieldid}.png")
                with open(outpath, "wb") as handler:
                    handler.write(img_data)

        print(f"Fields that are possible: {fieldids_total}")
        print(f"Of these have a reference: {fieldids_total_ref}")

    def check_galactic_latitude(self):
        """ """
        print("Not implemented yet")

    def get_best_obstime(self):
        """ """
        print("Not implemented yet")


# NEED TO INCLUDE ERROR CIRCLE CALCULATION
RA = 96.46
DEC = -4.33
NAME = "IC200926A"
ARRIVALTIME = "2020-09-26 07:54:11.621"
date = "2020-09-26"

plan = ObservationPlan(ra=RA, dec=DEC, name=NAME, arrivaltime=ARRIVALTIME)  # date=date)

plan.plot_target()
# plan.request_ztf_fields()
