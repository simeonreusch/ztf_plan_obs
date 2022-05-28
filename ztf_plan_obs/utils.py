#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import re
from astropy.time import Time
from astropy import units as u


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


def mjd_delta_to_seconds(mjd_start, mjd_end):
    """
    Convert t_end - t_start (duration of obs)
    given in mjd into a time delta in seconds
    """
    return round((mjd_end - mjd_start) * 86400)


def isotime_delta_to_seconds(isotime_start, isotime_end):
    """
    Convert t_end - t_start (duration of obs) given in iso-time
    into a time delta in seconds
    """

    mjd_start = isotime_to_mjd(isotime_start)
    mjd_end = isotime_to_mjd(isotime_end)

    return round((mjd_end - mjd_start) * 86400)


def isotime_to_mjd(isotime: str):
    """
    Convert time in iso-format to mjd
    """
    return float(Time(isotime, format="iso", scale="utc").mjd)


def mjd_to_isotime(mjd: float):
    """
    Convert time in mjd to iso-format
    """
    return Time(mjd, format="mjd", scale="utc").iso
