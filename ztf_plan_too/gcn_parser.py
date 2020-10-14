#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# GCN parsing code partially by Rober Stein (robert.stein@desy.de)
# License: BSD-3-Clause

import os, time, re
import pandas as pd
from astropy.time import Time
import requests


def get_gcn_circulars_archive():
    response = requests.get("https://gcn.gsfc.nasa.gov/gcn3_archive.html")

    gcns = []
    for line in response.text.splitlines():
        if "IceCube observation of a high-energy neutrino" in line:
            res = line.split(">")
            gcn_no = "".join([x for x in res[2] if x.isdigit()])
            long_name = re.findall(
                r"(IceCube-[12][0-9][0-9][0-9][0-3][0-9][A-Z])", line
            )[0]
            short_name = "IC" + long_name[8:]
            gcns.append((short_name, gcn_no))

    return gcns


def parse_gcn_circular(gcn_number):
    url = f"https://gcn.gsfc.nasa.gov/gcn3/{gcn_number}.gcn3"
    response = requests.get(url)
    splittext = response.text.splitlines()
    for i, line in enumerate(splittext):
        if ("RA" in line or "Ra" in line) and (
            "DEC" in splittext[i + 1] or "Dec" in splittext[i + 1]
        ):
            ra, ra_upper, ra_lower = parse_radec(line)
            dec, dec_upper, dec_lower = parse_radec(splittext[i + 1])
            ra_err = [ra_upper, ra_lower]
            dec_err = [dec_upper, dec_lower]
            return ra, ra_err, dec, dec_err


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
    date = latest["EVENT"]["Date"][0].replace("/", "-")
    obstime = latest["EVENT"]["Time UT"][0]
    ra = latest["OBSERVATION"]["RA [deg]"][0]
    dec = latest["OBSERVATION"]["Dec [deg]"][0]
    arrivaltime = Time(f"20{date} {obstime}")
    return ra, dec, arrivaltime
