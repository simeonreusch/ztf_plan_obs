#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# GCN parsing code partially by Robert Stein (robert.stein@desy.de)
# License: BSD-3-Clause

import os, time, re
import numpy as np
import pandas as pd
from astropy.time import Time
import requests


def get_gcn_circulars_archive(archive_no=None):
    if archive_no is None:
        response = requests.get("https://gcn.gsfc.nasa.gov/gcn3_archive.html")
    else:
        response = requests.get(
            f"https://gcn.gsfc.nasa.gov/gcn3_arch_old{archive_no}.html"
        )

    gcns = []
    _archive_numbers = []
    for line in response.text.splitlines():
        if "IceCube observation of a high-energy neutrino" in line:
            res = line.split(">")
            gcn_no = "".join([x for x in res[2] if x.isdigit()])
            long_name = re.findall(
                r"(IceCube-[12][0-9][0-9][0-9][0-3][0-9][A-Z])", line
            )[0]
            short_name = "IC" + long_name[8:]
            gcns.append((short_name, gcn_no))
        elif "gcn3_arch_old" in line:
            url = line.split('"')[1]
            _archive_no = int(url[13:].split(".")[0])
            _archive_numbers.append(_archive_no)

    if archive_no is not None:
        print(f"Processed archive number {archive_no}")

    return gcns, max(_archive_numbers)


def parse_gcn_circular(gcn_number):
    url = f"https://gcn.gsfc.nasa.gov/gcn3/{gcn_number}.gcn3"
    response = requests.get(url)
    returndict = {}
    mainbody_starts_here = 999
    splittext = response.text.splitlines()
    splittext = list(filter(None, splittext))
    for i, line in enumerate(splittext):
        if "SUBJECT" in line:
            name = line.split(" - ")[0].split(": ")[1]
            returndict.update({"name": name})
        elif "FROM" in line:
            base = line.split("at")[0].split(": ")[1].split(" ")
            author = [x for x in base if x != ""][1]
            returndict.update({"author": author})
        elif (
            ("RA" in line or "Ra" in line)
            and ("DEC" in splittext[i + 1] or "Dec" in splittext[i + 1])
            and i < mainbody_starts_here
        ):
            ra, ra_upper, ra_lower = parse_radec(line)
            dec, dec_upper, dec_lower = parse_radec(splittext[i + 1])
            ra_err = [ra_upper, -ra_lower]
            dec_err = [dec_upper, -dec_lower]
            returndict.update(
                {"ra": ra, "ra_err": ra_err, "dec": dec, "dec_err": dec_err}
            )
            mainbody_starts_here = i + 2
        elif ("Time" in line or "TIME" in line) and i < mainbody_starts_here:
            raw_time = [
                x for x in line.split(" ") if x not in ["Time", "", "UT", "UTC"]
            ][1]
            raw_time = "".join(
                [x for x in raw_time if np.logical_or(x.isdigit(), x in [":", "."])]
            )
            raw_date = name.split("-")[1][:6]
            ut_time = f"20{raw_date[0:2]}-{raw_date[2:4]}-{raw_date[4:6]}T{raw_time}"
            time = Time(ut_time, format="isot", scale="utc")
            returndict.update({"time": time})

    return returndict


def parse_radec(str):
    regex_findall = re.findall(r"[-+]?\d*\.\d+|\d+", str)
    if len(regex_findall) == 4:
        pos = float(regex_findall[0])
        pos_upper = float(regex_findall[1])
        pos_lower = float(regex_findall[1])
    elif len(regex_findall) == 5:
        pos, pos_upper, pos_lower = regex_findall[0:3]
        pos = float(pos)
        pos_upper = float(pos_upper.replace("+", ""))
        pos_lower = float(pos_lower.replace("-", ""))
    else:
        raise ParsingError(f"Could not parse GCN ra and dec")

    return pos, pos_upper, pos_lower


def parse_latest_gcn_notice():
    """ """
    url = "https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html"
    response = requests.get(url)
    table = pd.read_html(response.text)[0]
    latest = table.head(1)
    revision = latest["EVENT"]["Rev"][0]
    date = latest["EVENT"]["Date"][0].replace("/", "-")
    obstime = latest["EVENT"]["Time UT"][0]
    ra = latest["OBSERVATION"]["RA [deg]"][0]
    dec = latest["OBSERVATION"]["Dec [deg]"][0]
    arrivaltime = Time(f"20{date} {obstime}")
    return ra, dec, arrivaltime, revision


class ParsingError(Exception):
    pass
