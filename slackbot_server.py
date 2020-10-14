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
from slackbot import ObsBot

app = Flask(__name__)

slack_events_adapter = SlackEventAdapter(
    os.environ.get("SLACK_EVENTS_TOKEN"), "/slack/events", app
)
slack_web_client = WebClient(token=os.environ.get("SLACK_TOKEN"))


def do_obs_plan(channel, name, ra=None, dec=None, date=None, alertsource=None):
    """ """
    obs_bot = ObsBot(
        channel=channel, name=name, ra=ra, dec=dec, date=date, alertsource=alertsource
    )
    obs_bot.create_plot()

    # Post a text summary
    slack_web_client.chat_postMessage(
        channel=channel,
        text=obs_bot.summary,
    )
    if obs_bot.fields is not None:
        slack_web_client.chat_postMessage(
            channel=channel,
            text=f"Available fields: {obs_bot.fields}",
        )

    if obs_bot.summary != "Not observable due to airmass constraint":
        # Post the airmass plot
        imgpath_plot = f"{name}/{name}_airmass.png"
        imgdata_plot = open(imgpath_plot, "rb")
        slack_web_client.files_upload(
            file=imgdata_plot, filename=imgpath_plot, channels=channel
        )

        # Post the ZTF grid plots
        for field in obs_bot.fields:
            imgpath = f"{name}/{name}_grid_{field}.png"
            imgdata = open(imgpath, "rb")
            slack_web_client.files_upload(
                file=imgdata, filename=imgpath, channels=channel
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
            ra = None
            dec = None
            date = None
            radec_given = False
            alertsource = None
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

            if not radec_given:
                message = f"Hi there; creating your observability plot for *{name}*. One moment please."
            else:
                message = f"Hi there; creating your observability plot for *{name}*. You specified RA={ra} and Dec={dec}. One moment please."

            if ra is None:
                alertsource = "icecube"

            slack_web_client.chat_postMessage(
                channel=channel_id,
                text=message,
            )

            do_obs_plan(
                channel=channel_id,
                name=name,
                ra=ra,
                dec=dec,
                date=date,
                alertsource=alertsource,
            )


if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())

    app.run(host="130.255.78.114", port=3000)
