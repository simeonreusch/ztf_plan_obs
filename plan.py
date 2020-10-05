#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import time, os, warnings
import astropy
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, AltAz
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
        self.coordinates_galactic = self.coordinates.galactic
        self.target = ap.FixedTarget(name=self.name, coord=self.coordinates)
        self.palomar = Observer.at_site("Palomar", timezone="US/Pacific")
        self.now = Time(datetime.utcnow())
        self.date = date
        # if self.date is not None:
        #     moon_coords = astropy.coordinates.get_moon(time=Time(self.date + " 00:00:00.000000"), location=self.palomar.location)
        # else:
        #     moon_coords = astropy.coordinates.get_moon(time=self.now, location=self.palomar.location)
        # self.moon = ap.FixedTarget(name="moon", coord=moon_coords)

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

        # Obtain moon coordinates at Palomar for the full time window
        times = Time(self.start_obswindow + np.linspace(0, 24, 1000) * u.hour)
        moon_times = Time(self.start_obswindow + np.linspace(0, 24, 50) * u.hour)
        moon_coords = []

        for time in moon_times:
            moon_coord = astropy.coordinates.get_moon(
                time=time, location=self.palomar.location
            )
            moon_coords.append(moon_coord)
        self.moon = moon_coords

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
        print(f"RADEC = {self.coordinates.ra.deg} {self.coordinates.dec.deg}")
        print(f"Minimal airmass ({min_airmass:.2f}) at {min_airmass_time}")
        print("Separation from galactic plane:")
        print(f"{self.coordinates_galactic.b.deg:.2f} deg")

        print("Recommended observation times:")
        print(
            f"g-band: {self.time_shortener(self.g_band_recommended_time_start)} - {self.time_shortener(self.g_band_recommended_time_end)} [UTC]"
        )
        print(
            f"r-band: {self.time_shortener(self.r_band_recommended_time_start)} - {self.time_shortener(self.r_band_recommended_time_end)} [UTC]"
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

        ax = plot_altitude(
            self.target,
            self.palomar,
            time,
            brightness_shading=True,
            min_altitude=10,
        )

        # Plot a vertical line for the current time
        ax.axvline(Time(self.now).plot_date, color="black", label="now", ls="dotted")

        # Plot a vertical line for the neutrino arrival time if available
        if self.arrivaltime is not None:
            ax.axvline(
                Time(self.arrivaltime).plot_date,
                color="indigo",
                label="neutrino arrival",
                ls="dashed",
            )

        start, end = ax.get_xlim()

        plt.text(
            start,
            100,
            f"RADEC = {self.coordinates.ra.deg} {self.coordinates.dec.deg}\nSeparation from Galactic Plane: {self.coordinates_galactic.b.deg:.2f} deg\nRecommended observation times:\ng-band: {self.time_shortener(self.g_band_recommended_time_start)} --- {self.time_shortener(self.g_band_recommended_time_end)} [UTC]\n r-band: {self.time_shortener(self.r_band_recommended_time_start)} --- {self.time_shortener(self.r_band_recommended_time_end)} [UTC]",
        )

        if self.date is not None:
            ax.set_xlabel(f"{self.date} [UTC]")
        else:
            ax.set_xlabel(f"{self.now.datetime.date()} [UTC]")
        plt.grid(True, color="gray", linestyle="dotted", which="both", alpha=0.5)

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

        # Now we plot the moon altitudes and separation
        moon_altitudes = []
        moon_times = []
        moon_separations = []
        for moon in self.moon:
            moonalt = moon.transform_to(
                AltAz(obstime=moon.obstime, location=self.palomar.location)
            ).alt.deg
            moon_altitudes.append(moonalt)
            moon_times.append(moon.obstime.plot_date)
            separation = moon.separation(self.coordinates).deg
            moon_separations.append(separation)
        ax.plot(
            moon_times,
            moon_altitudes,
            color="orange",
            linestyle=(0, (1, 2)),
            label="moon",
        )

        # And we annotate the separations

        for i, moonalt in enumerate(moon_altitudes):
            if moonalt > 20 and i % 3 == 0:
                if moon_separations[i] < 20:
                    color = "red"
                else:
                    color = "green"
                ax.annotate(
                    f"{moon_separations[i]:.0f}",
                    xy=(moon_times[i], moonalt),
                    textcoords="data",
                    fontsize=6,
                    color=color,
                )

        # for i, xy in enumerate(zip(moon_times, moon_altitudes)):
        #     ax.annotate(f"{moon_separations[i]:.0f}", xy=xy, textcoords='data')

        x = np.linspace(start + 0.03, end + 0.03, 9)

        # Add recommended upper limit for airmass
        y = np.full((len(x), 1), 30)
        ax.errorbar(x, y, 2, color="red", lolims=True, fmt=" ")

        # Plot an airmass scale
        ax2 = ax.secondary_yaxis(
            "right", functions=(self.altitude_to_airmass, self.airmass_to_altitude)
        )
        altitude_ticks = np.linspace(10, 90, 9)
        airmass_ticks = np.round(self.altitude_to_airmass(altitude_ticks), 2)
        ax2.set_yticks(airmass_ticks)
        ax2.set_ylabel("Airmass")

        plt.tight_layout()
        plt.legend()
        outpath_png = os.path.join(self.name, f"{self.name}_airmass.png")
        outpath_pdf = os.path.join(self.name, f"{self.name}_airmass.pdf")
        plt.savefig(outpath_png, dpi=300, bbox_inches="tight")
        plt.savefig(outpath_pdf, bbox_inches="tight")

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

    @staticmethod
    def airmass_to_altitude(altitude):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            airmass = 90 - np.degrees(np.arccos(1 / altitude))
        return airmass

    @staticmethod
    def altitude_to_airmass(airmass):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            altitude = 1.0 / np.cos(np.radians(90 - airmass))
        return altitude


# NEED TO INCLUDE ERROR CIRCLE CALCULATION
NAME = "IC200929A"
# RA = 90.46
# DEC = -4.33
# ARRIVALTIME = "2020-09-26 07:54:11.621"
# date = "2020-09-29"

plan = ObservationPlan(name=NAME)  # , date=date
# , ra=RA, dec=DEC, arrivaltime=ARRIVALTIME

# )  # , date=date)

plan.plot_target()
# plan.request_ztf_fields()
