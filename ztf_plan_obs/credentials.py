#!/usr/bin/env python
# coding: utf-8

import os
import warnings
import logging
from ztfquery import io

# Manage ztfquery logins from environment variables


def load_credentials(name: str, token_based: bool = False):
    """ZTFquery wrapper for loading credentials."""
    return io._load_id_(name, token_based=token_based)


try:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        io.set_account(
            "irsa",
            username=os.environ["IRSA_USER"],
            password=os.environ["IRSA_PASSWORD"],
        )
        logging.info('Set up "irsa" credentials')

except KeyError:
    logging.info(
        'No Credentials for "IRSA" found in environment' "Assuming they are set."
    )

try:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        io.set_account("fritz", token=os.environ["FRITZ_TOKEN"], token_based=True)
        logging.info('Set up "Fritz" token')

except KeyError:
    logging.info("No token for Fritz API found in environment" "Assume it is set.")
