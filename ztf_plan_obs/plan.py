#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# GCN parsing code partially by Rober Stein (robert.stein@desy.de)
# License: BSD-3-Clause

import time, os, warnings
import astropy
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, AltAz
import os, time, re
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
from ztf_plan_obs import gcn_parser

icecube = ["IceCube", "IC", "icecube", "ICECUBE", "Icecube"]


class ObservationPlan:
    """ """

    def __init__(
        self,
        name: str,
        ra: float = None,
        dec: float = None,
        arrivaltime: str = None,
        date: str = None,
        alertsource: str = None,
        **kwargs,
    ):

        self.name = name
        self.arrivaltime = arrivaltime
        self.alertsource = alertsource
        self.ra_err = None
        self.dec_err = None
        self.warning = None
        self.observable = True

        use_archival = False

        if ra is None and self.alertsource in icecube:
            print("Parsing an IceCube alert")

            # Check if request is archival:
            archive = gcn_parser.get_gcn_circulars_archive()
            for archival_name, archival_number in archive:
                if name == archival_name:
                    gnc_nr = archival_number
                    use_archival = True
                    print("Archival data found, using these.")
                    break

            if use_archival:
                (
                    self.ra,
                    self.ra_err,
                    self.dec,
                    self.dec_err,
                ) = gcn_parser.parse_gcn_circular(archival_number)

            else:
                print(
                    "No archival data found (archive only stretches ~6 months back). Using newest notice!"
                )
                (
                    ra_notice,
                    dec_notice,
                    self.arrivaltime,
                ) = gcn_parser.parse_latest_gcn_notice()
                gcn_nr_latest = archive[0][1]
                (
                    ra_circ,
                    ra_err_circ,
                    dec_circ,
                    dec_err_circ,
                ) = gcn_parser.parse_gcn_circular(gcn_nr_latest)
                coords_notice = SkyCoord(
                    ra_notice * u.deg, dec_notice * u.deg, frame="icrs"
                )
                coords_circular = SkyCoord(
                    ra_circ * u.deg, dec_circ * u.deg, frame="icrs"
                )
                separation = coords_notice.separation(coords_circular).deg
                if separation < 1:
                    self.ra = ra_circ
                    self.dec = dec_circ
                    self.ra_err = ra_err_circ
                    self.dec_err = dec_err_circ
                else:
                    self.ra = ra_notice
                    self.dec = dec_notice

        else:
            self.ra = ra
            self.dec = dec

        self.coordinates = SkyCoord(self.ra * u.deg, self.dec * u.deg, frame="icrs")
        self.coordinates_galactic = self.coordinates.galactic
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

        self.twilight_evening = self.palomar.twilight_evening_astronomical(
            Time(self.start_obswindow), which="next"
        )
        self.twilight_morning = self.palomar.twilight_morning_astronomical(
            Time(self.start_obswindow), which="next"
        )

        indices_included = []
        airmasses_included = []
        times_included = []

        for index, t_mjd in enumerate(times.mjd):
            if (
                t_mjd > self.twilight_evening.mjd + 0.01
                and t_mjd < self.twilight_morning.mjd - 0.01
            ):
                if airmass[index] < 2.0:
                    indices_included.append(index)
                    airmasses_included.append(airmass[index])
                    times_included.append(times[index])

        if len(airmasses_included) == 0:
            self.observable = False
            # raise AirmassError("No observation due to airmass constraint!")

        if self.observable:
            min_airmass = np.min(airmasses_included)
            min_airmass_index = np.argmin(airmasses_included)
            min_airmass_time = times_included[min_airmass_index]

            distance_to_evening = min_airmass_time.mjd - self.twilight_evening.mjd
            distance_to_morning = self.twilight_morning.mjd - min_airmass_time.mjd

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
        if self.alertsource in icecube:
            summarytext = f"Name = IceCube-{self.name[2:]}\n"
        else:
            summarytext = f"Name = {self.name}\n"

        if self.ra_err is not None:
            summarytext += f"RA = {self.coordinates.ra.deg} + {self.ra_err[0]} - {self.ra_err[1]}\nDec = {self.coordinates.dec.deg} + {self.dec_err[0]} - {self.dec_err[1]}\n"
        else:
            summarytext += (
                f"RADEC = {self.coordinates.ra.deg} {self.coordinates.dec.deg}\n"
            )

        if self.observable:
            summarytext += (
                f"Minimal airmass ({min_airmass:.2f}) at {min_airmass_time}\n"
            )
        summarytext += f"Separation from galactic plane: = {self.coordinates_galactic.b.deg:.2f} deg\n"

        if self.observable:
            summarytext += "Recommended observation times:\n"
            summarytext += f"g-band: {self.time_shortener(self.g_band_recommended_time_start)} - {self.time_shortener(self.g_band_recommended_time_end)} [UTC]\n"
            summarytext += f"r-band: {self.time_shortener(self.r_band_recommended_time_start)} - {self.time_shortener(self.r_band_recommended_time_end)} [UTC]"

        print(summarytext)

        if not os.path.exists(self.name):
            os.makedirs(self.name)

        self.summarytext = summarytext

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
            # brightness_shading=True,
            min_altitude=10,
        )

        ax.axvspan(
            self.twilight_evening.plot_date,
            self.twilight_morning.plot_date,
            alpha=0.2,
            color="gray",
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
            # f"RADEC = {self.coordinates.ra.deg} {self.coordinates.dec.deg}\nSeparation from Galactic Plane: {self.coordinates_galactic.b.deg:.2f} deg\nRecommended observation times:\ng-band: {self.time_shortener(self.g_band_recommended_time_start)} --- {self.time_shortener(self.g_band_recommended_time_end)} [UTC]\n r-band: {self.time_shortener(self.r_band_recommended_time_start)} --- {self.time_shortener(self.r_band_recommended_time_end)} [UTC]",
            self.summarytext,
            fontsize=8,
        )

        if self.date is not None:
            ax.set_xlabel(f"{self.date} [UTC]")
        else:
            ax.set_xlabel(f"{self.now.datetime.date()} [UTC]")
        plt.grid(True, color="gray", linestyle="dotted", which="both", alpha=0.5)

        if self.observable:
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

        if self.observable is False:
            plt.text(
                0.5,
                0.5,
                "NOT OBSERVABLE",
                size=20,
                rotation=30.0,
                ha="center",
                va="center",
                bbox=dict(
                    boxstyle="round",
                    ec=(1.0, 0.5, 0.5),
                    fc=(1.0, 0.8, 0.8),
                ),
                transform=ax.transAxes,
            )

        plt.tight_layout()
        plt.legend()
        outpath_png = os.path.join(self.name, f"{self.name}_airmass.png")
        outpath_pdf = os.path.join(self.name, f"{self.name}_airmass.pdf")
        plt.savefig(outpath_png, dpi=300, bbox_inches="tight")
        plt.savefig(outpath_pdf, bbox_inches="tight")
        plt.close()

    def request_ztf_fields(self):
        """
        This looks at yupana.caltech.edu for the fields matching
        your location and downloads the camera grid plots for these
        """
        URL = "http://yupana.caltech.edu/cgi-bin/ptf/tb//zoc"
        image_url = "http://yupana.caltech.edu/marshals/tb//igmo_0_"
        image_urls = [image_url + f"{x}.png" for x in [0, 1, 2, 3]]

        objra = self.ra
        objdec = self.dec
        radius = 60

        fieldids_total = []
        fieldids_total_ref = []

        if self.ra_err is not None:
            radec_err = max(max(self.ra_err), max(self.dec_err))
            radius = 60 * radec_err

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
                img_data = requests.get(image_urls[index]).content
                outpath = os.path.join(self.name, f"{self.name}_grid_{fieldid}.png")
                with open(outpath, "wb") as handler:
                    handler.write(img_data)

        print(f"Fields that are possible: {fieldids_total}")
        print(f"Of these have a reference: {fieldids_total_ref}")

        return fieldids_total_ref

    def get_summary(self):
        return self.summarytext

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


class ParsingError(Exception):
    """Base class for parsing error"""

    pass


class AirmassError(Exception):
    """Base class for parsing error"""

    pass
