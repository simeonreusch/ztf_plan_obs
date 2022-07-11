#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import logging
import multiprocessing
import numpy as np
from ztfquery import io
from tqdm import tqdm

MARSHAL_BASEURL = "http://skipper.caltech.edu:8080/cgi-bin/growth/view_avro.cgi?name="

logger = logging.getLogger(__name__)


class FritzInfo:
    """Testing only"""

    def __init__(self, ztf_names):
        self.ztf_names = ztf_names

        self.queryresult = self.get_info()

    def get_info(self):
        from ztfquery import fritz

        returndict = {}

        object_count = len(self.ztf_names)

        queryresult = []

        for i, name in enumerate(tqdm(self.ztf_names)):
            query_res = fritz.download_alerts(name)
            queryresult.append(query_res)

        ras = []
        decs = []

        for entry in queryresult[0]:
            ras.append(entry["candidate"]["ra"])
            decs.append(entry["candidate"]["dec"])

        ra = np.median(ras)
        dec = np.median(decs)

        returndict.update({"ra": ra, "dec": dec})
        return returndict
