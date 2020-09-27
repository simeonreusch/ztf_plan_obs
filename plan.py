#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os, time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astroplan as ap
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroplan import Observer
from datetime import datetime
from astroplan.plots import plot_airmass, plot_altitude


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
            altitude_yaxis=True,
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


RA = 96.46
DEC = -4.33
NAME = "IC200926A"
ARRIVALTIME = "2020-09-26 07:54:11.621"

plan = ObservationPlan(ra=RA, dec=DEC, name=NAME, arrivaltime=ARRIVALTIME)

plan.plot_target()
