#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# GCN parsing code partially by Rober Stein (robert.stein@desy.de)
# License: BSD-3-Clause

import time, os, warnings, typing
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
from astroplan.plots import plot_finder_image
from datetime import datetime
from astroplan.plots import plot_airmass, plot_altitude
from ztfquery import fields, query
from ztf_plan_obs import gcn_parser
from shapely.geometry import Polygon

icecube = ["IceCube", "IC", "icecube", "ICECUBE", "Icecube"]
ztf = ["ZTF", "ztf"]


class PlanObservation:
    """ """

    def __init__(
        self,
        name: str,
        ra: float = None,
        dec: float = None,
        arrivaltime: str = None,
        date: str = None,
        max_airmass=2.0,
        observationlength: float = 300,
        bands: list = ["g", "r"],
        multiday: bool = False,
        alertsource: str = None,
        site: str = "Palomar",
        verbose: bool = True,
        **kwargs,
    ):

        self.name = name
        self.arrivaltime = arrivaltime
        self.alertsource = alertsource
        self.site = site
        self.max_airmass = max_airmass
        self.observationlength = observationlength
        self.bands = bands
        self.multiday = multiday
        self.ra_err = None
        self.dec_err = None
        self.warning = None
        self.observable = True
        self.rejection_reason = None
        self.datasource = None
        self.found_in_archive = False
        self.search_full_archive = False
        self.coverage = None
        self.recommended_field = None

        if ra is None and self.alertsource in icecube:
            if verbose:
                print("Parsing an IceCube alert")

            # Check if request is archival:
            archive, latest_archive_no = gcn_parser.get_gcn_circulars_archive()

            # Check if the alert is younger than latest archive entry
            archival_names = [entry[0] for entry in archive]
            archival_dates = [int(entry[2:-1]) for entry in archival_names]
            latest_archival = max(archival_dates)
            this_alert_date = int(self.name[2:-1])
            if this_alert_date > latest_archival:
                if verbose:
                    print(
                        "Alert too new, no GCN circular available yet. Using latest GCN notice"
                    )
            else:
                if verbose:
                    print("Alert info should be in GCN circular archive")
                self.search_full_archive = True

            if self.search_full_archive:
                self.search_match_in_archive(archive)

                # Well, if it's not in the latest archive, use the full
                # backwards search
                while self.found_in_archive is False:
                    archive, _ = gcn_parser.get_gcn_circulars_archive(latest_archive_no)
                    self.search_match_in_archive(archive)
                    latest_archive_no -= 1

            if self.found_in_archive:
                gcn_info = gcn_parser.parse_gcn_circular(self.gcn_nr)
                self.ra = gcn_info["ra"]
                self.ra_err = gcn_info["ra_err"]
                self.dec = gcn_info["dec"]
                self.dec_err = gcn_info["dec_err"]
                self.arrivaltime = gcn_info["time"]

            else:
                if verbose:
                    print("No archival GCN circular found. Using newest notice!")
                (
                    ra_notice,
                    dec_notice,
                    self.arrivaltime,
                    revision,
                ) = gcn_parser.parse_latest_gcn_notice()
                gcn_nr_latest = archive[0][1]
                gcn_info = gcn_parser.parse_gcn_circular(gcn_nr_latest)
                ra_circ = gcn_info["ra"]
                ra_err_circ = gcn_info["ra_err"]
                dec_circ = gcn_info["dec"]
                dec_err_circ = gcn_info["dec_err"]
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
                    self.datasource = f"GCN Circular {gcn_nr_latest}\n"
                else:
                    self.ra = ra_notice
                    self.dec = dec_notice
                    self.datasource = f"GCN Notice (Rev. {revision})\n"

        elif ra is None and self.alertsource in ztf:
            if is_ztf_name(name):
                print(f"{name} is a ZTF name. Looking in Fritz database for ra/dec")
                from ztf_plan_obs.fritzconnector import FritzInfo

                fritz = FritzInfo([name])

                self.ra = fritz.queryresult["ra"]
                self.dec = fritz.queryresult["dec"]

                self.datasource = "Fritz\n"

                if np.isnan(self.ra):
                    raise ValueError("Object apparently not found on Fritz")

                print("\nFound ZTF object information on Fritz")
        elif ra is None:
            raise ValueError("Please enter ra and dec")

        else:
            self.ra = ra
            self.dec = dec

        self.coordinates = SkyCoord(self.ra * u.deg, self.dec * u.deg, frame="icrs")
        self.coordinates_galactic = self.coordinates.galactic
        self.target = ap.FixedTarget(name=self.name, coord=self.coordinates)
        self.site = Observer.at_site(self.site, timezone="US/Pacific")
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
            ap.AirmassConstraint(max_airmass),
            ap.AtNightConstraint.twilight_astronomical(),
        ]

        # Obtain moon coordinates at Palomar for the full time window
        times = Time(self.start_obswindow + np.linspace(0, 24, 1000) * u.hour)
        moon_times = Time(self.start_obswindow + np.linspace(0, 24, 50) * u.hour)
        moon_coords = []

        for time in moon_times:
            moon_coord = astropy.coordinates.get_moon(
                time=time, location=self.site.location
            )
            moon_coords.append(moon_coord)
        self.moon = moon_coords

        airmass = self.site.altaz(times, self.target).secz
        airmass = np.ma.array(airmass, mask=airmass < 1)
        airmass = airmass.filled(fill_value=99)
        airmass = [x.value for x in airmass]

        self.twilight_evening = self.site.twilight_evening_astronomical(
            Time(self.start_obswindow), which="next"
        )
        self.twilight_morning = self.site.twilight_morning_astronomical(
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
            self.rejection_reason = "airmass"

        if np.abs(self.coordinates_galactic.b.deg) < 10:
            self.observable = False
            self.rejection_reason = "proximity to gal. plane"

        self.g_band_recommended_time_start = None
        self.g_band_recommended_time_end = None
        self.r_band_recommended_time_start = None
        self.r_band_recommended_time_end = None

        if self.observable:
            min_airmass = np.min(airmasses_included)
            min_airmass_index = np.argmin(airmasses_included)
            min_airmass_time = times_included[min_airmass_index]

            distance_to_evening = min_airmass_time.mjd - self.twilight_evening.mjd
            distance_to_morning = self.twilight_morning.mjd - min_airmass_time.mjd

            if distance_to_morning < distance_to_evening:
                if "g" in self.bands:
                    self.g_band_recommended_time_start = round_time(
                        min_airmass_time - self.observationlength * u.s - 0.5 * u.hour
                    )
                    self.g_band_recommended_time_end = (
                        self.g_band_recommended_time_start
                        + self.observationlength * u.s
                    )
                if "r" in self.bands:
                    self.r_band_recommended_time_start = round_time(
                        min_airmass_time - self.observationlength * u.s
                    )
                    self.r_band_recommended_time_end = (
                        self.r_band_recommended_time_start
                        + self.observationlength * u.s
                    )

            else:
                if "g" in self.bands:
                    self.g_band_recommended_time_start = round_time(
                        min_airmass_time + self.observationlength * u.s + 0.5 * u.hour
                    )
                    self.g_band_recommended_time_end = (
                        self.g_band_recommended_time_start
                        + self.observationlength * u.s
                    )
                if "r" in self.bands:
                    self.r_band_recommended_time_start = round_time(
                        min_airmass_time + self.observationlength * u.s
                    )
                    self.r_band_recommended_time_end = (
                        self.r_band_recommended_time_start
                        + self.observationlength * u.s
                    )
        if self.alertsource in icecube:
            summarytext = f"Name = IceCube-{self.name[2:]}\n"
        else:
            summarytext = f"Name = {self.name}\n"

        if self.ra_err:
            if self.ra_err[0]:
                summarytext += f"RA = {self.coordinates.ra.deg} + {self.ra_err[0]} - {self.ra_err[1]*-1}\nDec = {self.coordinates.dec.deg} + {self.dec_err[0]} - {self.dec_err[1]*-1}\n"
        else:
            summarytext += f"RADEC = {self.coordinates.ra.deg:.8f} {self.coordinates.dec.deg:.8f}\n"

        if self.datasource is not None:
            summarytext += f"Data source: {self.datasource}"

        if self.observable:
            summarytext += (
                f"Minimal airmass ({min_airmass:.2f}) at {min_airmass_time}\n"
            )
        summarytext += f"Separation from galactic plane: {self.coordinates_galactic.b.deg:.2f} deg\n"

        if self.site.name != "Palomar":
            summarytext += f"Site: {self.site.name}"

        if self.site.name == "Palomar":

            if self.observable and not self.multiday:
                summarytext += "Recommended observation times:\n"
                if "g" in self.bands:
                    gbandtext = f"g-band: {short_time(self.g_band_recommended_time_start)} - {short_time(self.g_band_recommended_time_end)} [UTC]"
                if "r" in self.bands:
                    rbandtext = f"r-band: {short_time(self.r_band_recommended_time_start)} - {short_time(self.r_band_recommended_time_end)} [UTC]"

                if (
                    "g" in bands
                    and "r" in bands
                    and self.g_band_recommended_time_start
                    < self.r_band_recommended_time_start
                ):
                    bandtexts = [gbandtext + "\n", rbandtext]
                elif (
                    "g" in bands
                    and "r" in bands
                    and self.g_band_recommended_time_start
                    > self.r_band_recommended_time_start
                ):
                    bandtexts = [rbandtext + "\n", gbandtext]
                elif "g" in bands and "r" not in bands:
                    bandtexts = [gbandtext]
                else:
                    bandtexts = [rbandtext]

                for item in bandtexts:
                    summarytext += item

        if verbose:
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
            self.site,
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

        midnight = min(self.twilight_evening, self.twilight_morning) + 0.5 * (
            max(self.twilight_evening, self.twilight_morning)
            - min(self.twilight_evening, self.twilight_morning)
        )

        ax.annotate(
            "Night",
            xy=[midnight.plot_date, 85],
            color="dimgray",
            ha="center",
            fontsize=12,
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
            self.summarytext,
            fontsize=8,
        )

        if self.date is not None:
            ax.set_xlabel(f"{self.date} [UTC]")
        else:
            ax.set_xlabel(f"{self.now.datetime.date()} [UTC]")
        plt.grid(True, color="gray", linestyle="dotted", which="both", alpha=0.5)

        if self.site.name == "Palomar":
            if self.observable:
                if "g" in self.bands:
                    ax.axvspan(
                        self.g_band_recommended_time_start.plot_date,
                        self.g_band_recommended_time_end.plot_date,
                        alpha=0.5,
                        color="green",
                    )
                if "r" in self.bands:
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
                AltAz(obstime=moon.obstime, location=self.site.location)
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
        y = np.full((len(x), 0), 30)
        y = np.ones(len(x)) * 30

        ax.errorbar(x, y, 2, color="red", lolims=True, fmt=" ")

        # Plot an airmass scale
        ax2 = ax.secondary_yaxis(
            "right", functions=(self.altitude_to_airmass, self.airmass_to_altitude)
        )
        altitude_ticks = np.linspace(10, 90, 9)
        airmass_ticks = np.round(self.altitude_to_airmass(altitude_ticks), 2)
        ax2.set_yticks(airmass_ticks)
        ax2.set_ylabel("Airmass")

        if self.observable:
            plt.legend()

        if self.observable is False:
            plt.text(
                0.5,
                0.5,
                f"NOT OBSERVABLE\ndue to {self.rejection_reason}",
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

        if self.site.name == "Palomar":
            outpath_png = os.path.join(self.name, f"{self.name}_airmass.png")
            outpath_pdf = os.path.join(self.name, f"{self.name}_airmass.pdf")
        else:
            outpath_png = os.path.join(
                self.name, f"{self.name}_airmass_{self.site.name}.png"
            )
            outpath_pdf = os.path.join(
                self.name, f"{self.name}_airmass_{self.site.name}.pdf"
            )
        plt.savefig(outpath_png, dpi=300, bbox_inches="tight")
        plt.savefig(outpath_pdf, bbox_inches="tight")

        return ax

    def search_match_in_archive(self, archive):
        """ """
        for archival_name, archival_number in archive:
            if self.name == archival_name:
                self.gcn_nr = archival_number
                self.found_in_archive = True
                self.datasource = f"GCN Circular {self.gcn_nr}\n"
                print("Archival data found, using these.")

    def request_ztf_fields(self, plot=True):
        """
        This looks at yupana.caltech.edu for the fields matching
        your location and downloads the camera grid plots for these
        """

        # URL = "http://yupana.caltech.edu/cgi-bin/ptf/tb//zoc"
        # image_url = "http://yupana.caltech.edu/marshals/tb//igmo_0_"
        # image_urls = [image_url + f"{x}.png" for x in [0, 1, 2, 3]]

        objra = self.ra
        objdec = self.dec

        radius = 0

        fieldids = list(fields.get_fields_containing_target(ra=self.ra, dec=self.dec))
        fieldids_ref = []

        zq = query.ZTFQuery()
        querystring = f"field={fieldids[0]}"

        if len(fieldids) > 1:
            for f in fieldids[1:]:
                querystring += f" OR field={f}"

        print(
            f"Checking IPAC if references are available in g- and r-band for fields {fieldids}"
        )

        zq.load_metadata(kind="ref", sql_query=querystring)
        mt = zq.metatable

        for f in mt.field.unique():
            d = {k: k in mt["filtercode"].values for k in ["zg", "zr", "zi"]}
            if d["zg"] == True and d["zr"] == True:
                fieldids_ref.append(f)

        print(f"Fields that contain target: {fieldids}")
        print(f"Of these have a reference: {fieldids_ref}")

        self.fieldids_ref = fieldids_ref

        if plot:
            self.plot_fields()

        return fieldids_ref

    def plot_fields(self):
        """
        Plot the ZTF field(s) with the target
        """
        ccds = fields._CCD_COORDS

        coverage = {}

        for f in self.fieldids_ref:
            centroid = fields.get_field_centroid(f)

            fig, ax = plt.subplots(dpi=300)

            ax.set_aspect("equal")

            ccd_polygons = []
            covered_area = 0

            for c in ccds.CCD.unique():
                ccd = ccds[ccds.CCD == c][["EW", "NS"]].values
                ccd_draw = Polygon(ccd + centroid)
                ccd_polygons.append(ccd_draw)
                x, y = ccd_draw.exterior.xy
                ax.plot(x, y, color="black")

            if self.ra_err:
                # Create errorbox
                ul = [self.ra + self.ra_err[1], self.dec + self.dec_err[0]]
                ur = [self.ra + self.ra_err[0], self.dec + self.dec_err[1]]
                ll = [self.ra + self.ra_err[1], self.dec + self.dec_err[1]]
                lr = [self.ra + self.ra_err[0], self.dec + self.dec_err[0]]

                errorbox = Polygon([ul, ll, ur, lr])
                x, y = errorbox.exterior.xy

                ax.plot(x, y, color="red")

                for ccd in ccd_polygons:
                    covered_area += errorbox.intersection(ccd).area

                cov = covered_area / errorbox.area * 100

                coverage.update({f: cov})

            ax.scatter([self.ra], [self.dec], color="red")

            ax.set_xlabel("RA")
            ax.set_ylabel("Dec")
            if self.ra_err:
                ax.set_title(f"Field {f} (Coverage: {cov:.2f}%)")
            else:
                ax.set_title(f"Field {f}")
            plt.tight_layout()

            outpath_png = os.path.join(self.name, f"{self.name}_grid_{f}.png")

            fig.savefig(outpath_png, dpi=300)

        self.coverage = coverage

        max_coverage_field = max(coverage, key=coverage.get)

        self.recommended_field = max_coverage_field

    def plot_finding_chart(self):
        """ """
        ax, hdu = plot_finder_image(
            self.target,
            fov_radius=2 * u.arcmin,
            survey="DSS2 Blue",
            grid=True,
            reticle=False,
        )
        outpath_png = os.path.join(self.name, f"{self.name}_finding_chart.png")
        plt.savefig(outpath_png, dpi=300)
        plt.close()

    def get_summary(self):
        return self.summarytext

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


def is_ztf_name(name):
    """
    Checks if a string adheres to the ZTF naming scheme
    """
    return re.match("^ZTF[1-2]\d[a-z]{7}$", name)


def is_icecube_name(name):
    """
    Checks if a string adheres to the IceCube naming scheme
    (e.g. IC201021B)
    """
    return re.match(
        "^IC((\d{2}((0[13578]|1[02])(0[1-9]|[12]\d|3[01])|(0[13456789]|1[012])(0[1-9]|[12]\d|30)|02(0[1-9]|1\d|2[0-8])))|([02468][048]|[13579][26])0229)[a-zA-Z]$",
        name,
    )


def round_time(time):
    """
    Better readable time - round to next minute
    """
    secs = float(str(time)[-6:])
    if secs < 30:
        time_rounded = time - secs * u.s
    else:
        time_rounded = time + (60 - secs) * u.s
    return time_rounded


def short_time(time):
    """
    Better readable time - remove subseconds
    """
    return str(time)[:-4]


class ParsingError(Exception):
    """Base class for parsing error"""

    pass


class AirmassError(Exception):
    """Base class for parsing error"""

    pass
