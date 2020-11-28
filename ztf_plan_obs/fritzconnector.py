#!/usr/bin/env python3
# Author: Simeon Reusch (simeon.reusch@desy.de)
# License: BSD-3-Clause

import multiprocessing
import numpy as np
from ztfquery import io
from astropy.utils.console import ProgressBar

MARSHAL_BASEURL = "http://skipper.caltech.edu:8080/cgi-bin/growth/view_avro.cgi?name="


class MarshalInfo:
    """ """

    def __init__(self, ztf_names, nprocess=16, logger=None):
        print("LALALALALALLALA")
        import requests
        import pandas as pd

        auth = io._load_id_("marshal")
        urls = []
        for ztf_name in ztf_names:
            url = MARSHAL_BASEURL + ztf_name
            urls.append(url)
        object_count = len(ztf_names)
        auth_ = [auth] * object_count
        from astropy.utils.console import ProgressBar

        bar = ProgressBar(object_count)
        results = []
        with multiprocessing.Pool(nprocess) as p:
            for index, result in enumerate(
                p.map(self.get_info_multiprocessor, zip(ztf_names, urls, auth_))
            ):
                bar.update(index)
                results.append(result)
            bar.update(object_count)
        self.queryresult = results

    @staticmethod
    def get_info_multiprocessor(args):
        """ """
        import requests
        import pandas as pd

        ztf_name, url, auth = args
        request = requests.get(url, auth=auth)
        tables = pd.read_html(request.content)
        mtb = tables[len(tables) - 1]
        ndet = len(mtb)

        if ndet == 0:
            ra = 999
            dec = 999
            jd = 999
        else:
            ra = np.zeros(ndet)
            dec = np.zeros(ndet)
            jd = np.zeros(ndet)
            mag = np.full(ndet, 99.0)
            magerr = np.zeros(ndet)
            maglim = np.zeros(ndet)
            jd = np.zeros(ndet)
            fid = np.full(ndet, 99)
            magzp = np.zeros(ndet)
            magzp_err = np.zeros(ndet)

            for i in range(ndet):
                isdiffpos = True
                try:
                    line = mtb.values[i][0].split(",")
                except:
                    print(mtb.values[i][0])
                for j in range(len(line)):
                    if line[j][:14] == '  "isdiffpos":':
                        isdiffpos = str(line[j].split(":")[1])
                        if isdiffpos[2:-1] == "f":
                            isdiffpos = False
                    if line[j][:7] == '  "ra":':
                        ra[i] = float(line[j].split(":")[1])
                    elif line[j][:8] == '  "dec":':
                        dec[i] = float(line[j].split(":")[1])

                # Throw away all alert datapoints
                # with negative diff images
                if isdiffpos == False:
                    ra[i] = 0

            ras = ra[ra != 0]
            decs = dec[ra != 0]
            jds = jd[ra != 0]
            ind = np.argsort(jds)
            ra_median = np.median(ras[ind])
            dec_median = np.median(decs[ind])

        return ra_median, dec_median


class FritzInfo:
    """ Testing only """

    def __init__(self, ztf_names):
        self.ztf_names = ztf_names

        self.queryresult = self.get_info()

    def get_info(self):
        from ztfquery import fritz

        returndict = {}

        object_count = len(self.ztf_names)
        bar = ProgressBar(object_count)

        queryresult = []

        for i, name in enumerate(self.ztf_names):
            query_res = fritz.download_alerts(name)
            queryresult.append(query_res)
            bar.update(i)

        bar.update(object_count)

        ras = []
        decs = []

        for entry in queryresult[0]:
            ras.append(entry["candidate"]["ra"])
            decs.append(entry["candidate"]["dec"])

        ra = np.median(ras)
        dec = np.median(decs)

        returndict.update({"ra": ra, "dec": dec})
        return returndict
