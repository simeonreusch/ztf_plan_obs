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
    io.set_account(
        "irsa", username=os.environ["IRSA_USER"], password=os.environ["IRSA_PASSWORD"]
    )
    logging.info('Set up "irsa" credentials')

except KeyError:
    logging.info(
        'No Credentials for "irsa" found in environment' "Assuming they are set."
    )
