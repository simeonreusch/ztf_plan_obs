#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import os
import logging
from flask import Flask
from slack import WebClient
from slackeventsapi import SlackEventAdapter
from obsbot import ObsBot

app = Flask(__name__)

slack_events_adapter = SlackEventAdapter(
    os.environ.get("SLACK_EVENTS_TOKEN"), "/slack/events", app
)
slack_web_client = WebClient(token=os.environ.get("SLACK_TOKEN"))


def do_obs_plan(channel, name):
    """ """
    obs_bot = ObsBot(channel, name)
    obs_bot.create_plot()

    imgpath = f"{name}/{name}_airmass.png"
    imgdata = open(imgpath, "rb")
    slack_web_client.files_upload(file=imgdata, filename=imgpath, channels=channel)


ts_old = []


@slack_events_adapter.on("message")
def message(payload):
    """ """
    event = payload.get("event", {})
    text = event.get("text")
    user = event.get("user")
    ts = event.get("ts")
    name = "IC200929A"
    if ts not in ts_old:
        ts_old.append(ts)
        print(text)
        print(user)
        print(ts)
        print("############################")

        if "plan" in text.lower():
            channel_id = event.get("channel")
            slack_web_client.chat_postMessage(
                channel=channel_id,
                text=f"Hi there; creating your observability plot. One moment please.",
            )
            do_obs_plan(channel_id, name)


if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())

    app.run(host="130.255.78.114", port=3000)
