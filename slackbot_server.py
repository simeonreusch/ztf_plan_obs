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
from ztf_plan_obs.api import Queue

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
    submit_trigger=False,
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
        submit_trigger=submit_trigger,
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
    if slack_bot.recommended_field is not None:
        slack_web_client.chat_postMessage(
            channel=channel,
            text=f"Recommended field: {slack_bot.recommended_field} ({slack_bot.coverage[slack_bot.recommended_field]:.2f}% coverage)",
        )

    if slack_bot.summary != "Not observable due to airmass constraint" and not multiday:
        # Post the airmass plot
        if site == "Palomar":
            imgpath_plot = f"{name}/{name}_airmass.png"
        else:
            imgpath_plot = f"{name}/{name}_airmass_{site}.png"
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


def get_submitted_too() -> list:
    q = Queue(user="DESY")
    existing_too_queue = q.get_too_queues(names_only=True)
    message = ""
    for entry in existing_too_queue:
        message += f"{entry}\n"
    message = message[:-1]
    return message


def get_submitted_full() -> list:
    q = Queue(user="DESY")
    existing_queue = q.get_all_queues(names_only=True)
    message = ""
    for entry in existing_queue:
        message += f"{entry}\n"
    message = message[:-1]
    return message


def delete_trigger(triggername) -> None:
    q = Queue(user="DESY")
    q.delete_trigger(triggername)


def fuzzy_parameters(param_list) -> list:
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
        logger.info(split_text)

        if len(split_text) == 0:
            return

        elif split_text[0] == "Plan" or split_text[0] == "plan":
            do_plan = True
            ra = None
            dec = None
            date = None
            radec_given = False
            multiday = False
            submit_trigger = False
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

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(
                    ["submit", "trigger", "Submit", "Trigger", "SUBMIT", "TRIGGER"]
                ):
                    submit_trigger = True

            if not radec_given:
                if not multiday:
                    if date:
                        message = f"Hi there; creating your observability plot for *{name}*. Starting date is {date}. One moment please."
                    else:
                        message = f"Hi there; creating your observability plot for *{name}*. Starting date is today. One moment please."
                else:
                    if date:
                        if not submit_trigger:
                            message = f"Hi there; creating your multiday observability plot for *{name}*. Starting date is {date}. One moment please."
                        else:
                            message = f"Hi there; creating your multiday observability plot for *{name}*. Starting date is {date}.\n\n!! The full multiday plan will be triggered !!\nPlease check the ZTF queue with 'Queue -get'."
                    else:
                        if not submit_trigger:
                            message = f"Hi there; creating your multiday observability plot for *{name}*. Starting date is today. One moment please."
                        else:
                            message = f"Hi there; creating your multiday observability plot for *{name}*. Starting date is today.\n\n!! The full multiday plan will be triggered !!\nPlease check the ZTF queue with 'Queue -get'."
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
                from ztf_plan_obs.utils import is_icecube_name, is_ztf_name

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
                    submit_trigger=submit_trigger,
                    alertsource=alertsource,
                    site=site,
                )

        elif split_text[0] in ["QUEUE", "Queue", "queue"]:

            channel_id = event.get("channel")

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(["get"]):
                    queue = get_submitted_too()
                    if len(queue) == 0:
                        message = "Currently, no ToO triggers are in the ZTF observation queue."
                    else:
                        message = f"The current ZTF ToO observation queue:\n{queue}"
                    slack_web_client.chat_postMessage(
                        channel=channel_id,
                        text=message,
                    )

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(["getfull"]):
                    queue = get_submitted_full()
                    message = f"The complete current ZTF observation queue:\n{queue}"
                    slack_web_client.chat_postMessage(
                        channel=channel_id,
                        text=message,
                    )

            for i, parameter in enumerate(split_text):
                if parameter in fuzzy_parameters(["delete"]):
                    triggername = split_text[i + 1]
                    delete_trigger(triggername)
                    message = (
                        f"Deleting the following trigger from the queue:\n{triggername}"
                    )
                    slack_web_client.chat_postMessage(
                        channel=channel_id,
                        text=message,
                    )


if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())

    app.run(host="168.119.229.141", port=3000)
