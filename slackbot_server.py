#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os, datetime, logging
from astropy.time import Time
from astropy import units as u

from datetime import datetime
from flask import Flask
from slack import WebClient
from slackeventsapi import SlackEventAdapter
from slackbot import Slackbot
from astropy.coordinates import EarthLocation

app = Flask(__name__)

slack_events_adapter = SlackEventAdapter(
    os.environ.get("SLACK_EVENTS_TOKEN"), "/slack/events", app
)
slack_web_client = WebClient(token=os.environ.get("SLACK_TOKEN"))


def do_obs_plan(
    channel,
    name,
    ra=None,
    dec=None,
    date=None,
    multiday=False,
    alertsource=None,
    site=None,
):
    """ """
    slack_bot = Slackbot(
        channel=channel,
        name=name,
        ra=ra,
        dec=dec,
        date=date,
        multiday=multiday,
        alertsource=alertsource,
        site=site,
    )
    slack_bot.create_plot()

    # Post a text summary
    slack_web_client.chat_postMessage(
        channel=channel,
        text=slack_bot.summary,
    )
    if slack_bot.fields is not None:
        slack_web_client.chat_postMessage(
            channel=channel,
            text=f"Available fields: {slack_bot.fields}",
        )

    if slack_bot.summary != "Not observable due to airmass constraint" and not multiday:
        # Post the airmass plot
        imgpath_plot = f"{name}/{name}_airmass.png"
        imgdata_plot = open(imgpath_plot, "rb")
        slack_web_client.files_upload(
            file=imgdata_plot, filename=imgpath_plot, channels=channel
        )

    if slack_bot.summary != "Not observable due to airmass constraint":
        # Post the ZTF grid plots
        if site in [None, "Palomar"]:
            for field in slack_bot.fields:
                imgpath = f"{name}/{name}_grid_{field}.png"
                imgdata = open(imgpath, "rb")
                slack_web_client.files_upload(
                    file=imgdata, filename=imgpath, channels=channel
                )
        if multiday:
            imgpath_plot = f"{name}/{name}_multiday.pdf"
            imgdata_plot = open(imgpath_plot, "rb")
            slack_web_client.files_upload(
                file=imgdata_plot, filename=imgpath_plot, channels=channel
            )

            slack_web_client.chat_postMessage(
                channel=channel,
                text=slack_bot.multiday_summary,
            )


def fuzzy_parameters(param_list):
    """ """
    fuzzy_parameters = []
    for param in param_list:
        for character in ["", "-", "--", "–"]:
            fuzzy_parameters.append(f"{character}{param}")
    return fuzzy_parameters


ts_old = []


@slack_events_adapter.on("message")
def message(payload):
    """ """
    event = payload.get("event", {})
    text = event.get("text")
    user = event.get("user")
    ts = event.get("ts")
    if ts not in ts_old:
        ts_old.append(ts)

        text = text.replace("*", "")
        split_text = text.split()

        if len(split_text) == 0:
            return

        if split_text[0] == "Plan" or split_text[0] == "plan":
            do_plan = True
            ra = None
            dec = None
            date = None
            radec_given = False
            multiday = False
            alertsource = None
            site = "Palomar"
            channel_id = event.get("channel")
            name = split_text[1]

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(["ra", "RA", "Ra"]):
                    ra = float(split_text[i + 1])
                if parameter in fuzzy_parameters(["dec", "Dec", "DEC"]):
                    dec = float(split_text[i + 1])
                    radec_given = True

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(["date", "DATE", "Date"]):
                    date = split_text[i + 1]

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(["tomorrow", "TOMORROW", "Tomorrow"]):
                    tomorrow = Time(datetime.utcnow()) + 24 * u.h
                    date = str(tomorrow.datetime.date())

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(
                    ["multiday", "MULTIDAY", "Multiday", "multi", "MULTI", "Multi"]
                ):
                    multiday = True

            if not radec_given:
                if not multiday:
                    if date:
                        message = f"Hi there; creating your observability plot for *{name}*. Starting date is {date}. One moment please."
                    else:
                        message = f"Hi there; creating your observability plot for *{name}*. Starting date is today. One moment please."
                else:
                    if date:
                        f"Hi there; creating your multiday observability plot for *{name}*. Starting date is {date}. One moment please."
                    else:
                        message = f"Hi there; creating your multiday observability plot for *{name}*. Starting date is today. One moment please."
            else:
                if not multiday:
                    if date:
                        message = f"Hi there; creating your observability plot for *{name}*. You specified RA={ra} and Dec={dec}. Starting date is {date}. One moment please."
                    else:
                        message = f"Hi there; creating your observability plot for *{name}*. You specified RA={ra} and Dec={dec}. Starting date is today. One moment please."
                else:
                    if date:
                        message = f"Hi there; creating your multiday observability plot for *{name}*. You specified RA={ra} and Dec={dec}. Starting date is {date}. One moment please."
                    else:
                        message = f"Hi there; creating your multiday observability plot for *{name}*. You specified RA={ra} and Dec={dec}. Starting date is today. One moment please."

            available_sites = EarthLocation.get_site_names()
            available_sites_reformatted = [
                entry.replace(" ", "_") for entry in available_sites if entry is not ""
            ]

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(
                    ["site", "Site", "telescope", "Telescope"]
                ):
                    site = str(split_text[i + 1]).replace("_", " ")
                    if site not in available_sites:
                        message = f"Your site/telescope needs to be in the following list: {available_sites_reformatted}"
                        do_plan = False
                    else:
                        message += f" Chosen site: {site}"

            if ra is None:
                from ztf_plan_obs.plan import is_icecube_name, is_ztf_name

                if is_icecube_name(name):
                    alertsource = "icecube"

                elif is_ztf_name(name):
                    alertsource = "ZTF"

                else:
                    message = f"When not giving radec, you have to provide an IceCube name (ICYYMMDD[A-Z]) or a ZTF name (ZTFYY[7*a-z])"
                    do_plan = False

            slack_web_client.chat_postMessage(
                channel=channel_id,
                text=message,
            )

            if do_plan:
                do_obs_plan(
                    channel=channel_id,
                    name=name,
                    ra=ra,
                    dec=dec,
                    date=date,
                    multiday=multiday,
                    alertsource=alertsource,
                    site=site,
                )


if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())

    app.run(host="168.119.229.141", port=3000)
