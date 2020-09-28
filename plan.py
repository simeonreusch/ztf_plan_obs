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
from astroplan import Observer, is_observable
from datetime import datetime
from astroplan.plots import plot_airmass, plot_altitude
import requests
from bs4 import BeautifulSoup
from ztfquery import fields


class ObservationPlan:
    """ """

    def __init__(
        self,
        name: str,
        ra: float = None,
        dec: float = None,
        arrivaltime: str = None,
        date: str = None,
        **kwargs,
    ):

        # if ra is None:

        self.name = name
        self.arrivaltime = arrivaltime

        if ra is None:
            self.parse_gcn()
        else:
            self.ra = ra
            self.dec = dec

        self.coordinates = SkyCoord(self.ra * u.deg, self.dec * u.deg, frame="icrs")
        self.target = ap.FixedTarget(name=self.name, coord=self.coordinates)
        self.palomar = Observer.at_site("Palomar", timezone="US/Pacific")
        self.now = Time(datetime.utcnow())
        self.date = date

        if self.date is not None:
            self.start_obswindow = Time(self.date + " 00:00:00.000000")

        else:
            self.start_obswindow = Time(
                str(self.now.datetime.date()) + " 00:00:00.000000"
            )

        self.end_obswindow = Time(self.start_obswindow.mjd + 1, format="mjd").iso

        constraints = [
            ap.AltitudeConstraint(20 * u.deg, 90 * u.deg),
            ap.AirmassConstraint(2.5),
            ap.AtNightConstraint.twilight_astronomical(),
        ]

        times = Time(self.start_obswindow + np.linspace(0, 24, 1000) * u.hour)

        airmass = self.palomar.altaz(times, self.target).secz
        airmass = np.ma.array(airmass, mask=airmass < 1)
        airmass = airmass.filled(fill_value=99)
        airmass = [x.value for x in airmass]

        twilight_evening = self.palomar.twilight_evening_astronomical(
            Time(self.start_obswindow), which="next"
        )
        twilight_morning = self.palomar.twilight_morning_astronomical(
            Time(self.start_obswindow), which="next"
        )

        indices_included = []
        airmasses_included = []
        times_included = []

        for index, t_mjd in enumerate(times.mjd):
            if (
                t_mjd > twilight_evening.mjd + 0.01
                and t_mjd < twilight_morning.mjd - 0.01
            ):
                if airmass[index] < 2.5:
                    indices_included.append(index)
                    airmasses_included.append(airmass[index])
                    times_included.append(times[index])

        if len(airmasses_included) == 0:
            raise Exception("No observation possible!")

        min_airmass = np.min(airmasses_included)
        min_airmass_index = np.argmin(airmasses_included)
        min_airmass_time = times_included[min_airmass_index]

        distance_to_evening = min_airmass_time.mjd - twilight_evening.mjd
        distance_to_morning = twilight_morning.mjd - min_airmass_time.mjd

        print(f"Minimal airmass ({min_airmass:.2f}) at {min_airmass_time}")

        if distance_to_morning < distance_to_evening:
            self.g_band_recommended_time_start = (
                min_airmass_time - 300 * u.s - 0.5 * u.hour
            )
            self.g_band_recommended_time_end = (
                self.g_band_recommended_time_start + 300 * u.s
            )
            self.r_band_recommended_time_start = min_airmass_time - 300 * u.s
            self.r_band_recommended_time_end = (
                self.r_band_recommended_time_start + 300 * u.s
            )

        else:
            self.g_band_recommended_time_start = (
                min_airmass_time + 300 * u.s + 0.5 * u.hour
            )
            self.g_band_recommended_time_end = (
                self.g_band_recommended_time_start + 300 * u.s
            )
            self.r_band_recommended_time_start = min_airmass_time + 300 * u.s
            self.r_band_recommended_time_end = (
                self.r_band_recommended_time_start + 300 * u.s
            )

        print("Recommended observation times:")
        print(
            f"g-band: {self.time_shortener(self.g_band_recommended_time_start)} --- {self.time_shortener(self.g_band_recommended_time_end)}"
        )
        print(
            f"r-band: {self.time_shortener(self.r_band_recommended_time_start)} --- {self.time_shortener(self.r_band_recommended_time_end)}"
        )

        if not os.path.exists(self.name):
            os.makedirs(self.name)

    def parse_gcn(self):
        """ """
        URL = "https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html"
        response = requests.get(URL)
        table = pd.read_html(response.text)[0]
        latest = table.head(1)
        date = latest["EVENT"]["Date"][0].replace("/", "-")
        obstime = latest["EVENT"]["Time UT"][0]
        ra = latest["OBSERVATION"]["RA [deg]"][0]
        dec = latest["OBSERVATION"]["Dec [deg]"][0]
        arrivaltime = Time(f"20{date} {obstime}")
        self.arrivaltime = arrivaltime
        self.ra = ra
        self.dec = dec

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
                color="indigo",
                label="neutrino arrival",
                ls="dotted",
            )

        start, end = ax.get_xlim()

        plt.text(
            start,
            0.8,
            f"Recommended observation times:\n\ng-band: {self.time_shortener(self.g_band_recommended_time_start)} --- {self.time_shortener(self.g_band_recommended_time_end)}\nr-band: {self.time_shortener(self.r_band_recommended_time_start)} --- {self.time_shortener(self.r_band_recommended_time_end)}",
        )

        if self.date is not None:
            ax.set_xlabel(f"{self.date} [UTC]")
        else:
            ax.set_xlabel(f"{self.now.datetime.date()} [UTC]")
        plt.grid(True, color="brown", linestyle="dotted", which="both")

        ax.axvspan(
            self.g_band_recommended_time_start.plot_date,
            self.g_band_recommended_time_end.plot_date,
            alpha=0.5,
            color="green",
        )

        ax.axvspan(
            self.r_band_recommended_time_start.plot_date,
            self.r_band_recommended_time_end.plot_date,
            alpha=0.5,
            color="red",
        )

        x = np.linspace(start + 0.03, end + 0.03, 9)

        # Add recommended upper limit for airmass
        y = np.full((len(x), 1), 2)
        ax.errorbar(x, y, 0.05, color="red", uplims=True, fmt=" ")

        plt.tight_layout()
        outpath_png = os.path.join(self.name, f"{self.name}_airmass.png")
        outpath_pdf = os.path.join(self.name, f"{self.name}_airmass.pdf")
        plt.savefig(outpath_png)
        plt.savefig(outpath_pdf)

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

    @staticmethod
    def time_shortener(time):
        time_short = str(time)[:-7] + ":00"
        return time_short


# NEED TO INCLUDE ERROR CIRCLE CALCULATION
NAME = "IC200926A"
RA = 90.46
DEC = -4.33
ARRIVALTIME = "2020-09-26 07:54:11.621"
# date = "2020-09-26"

plan = ObservationPlan(
    name=NAME, ra=RA, dec=DEC, arrivaltime=ARRIVALTIME
)  # , date=date)

plan.plot_target()
plan.request_ztf_fields()
