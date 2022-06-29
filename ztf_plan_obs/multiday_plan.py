#!/usr/bin/env python3
import matplotlib.pyplot as plt
import os
from datetime import datetime, date
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
from astropy.time import Time
from astropy import units as u
from ztf_plan_obs.plan import PlanObservation
from ztf_plan_obs.utils import (
    round_time,
    short_time,
    isotime_delta_to_seconds,
    isotime_to_mjd,
    mjd_to_isotime,
)

NIGHTS = [1, 2, 3, 5, 7, 9]
SHORT_NIGHTS = NIGHTS[1:]
ONE_FILTER_NIGHTS = NIGHTS[1:-1]


class MultiDayObservation:
    """ """

    def __init__(
        self,
        name: str,
        ra: float = None,
        dec: float = None,
        startdate=None,
        verbose: bool = True,
        **kwargs,
    ):

        self.name = name
        self.ra = ra
        self.dec = dec

        self.triggers: list = []

        today = date.today()
        now = datetime.now()

        if self.ra is None:
            plan_initial = PlanObservation(name=name, alertsource="icecube")
        else:
            plan_initial = PlanObservation(name=name, ra=self.ra, dec=self.dec)

        if startdate is None:
            first_obs = plan_initial.g_band_recommended_time_start
            first_obs_day = Time(
                first_obs, format="iso", scale="utc", out_subfmt="date"
            )
            next_days = [(first_obs_day + i - 1).value for i in NIGHTS]
        else:
            startdate_astropy = Time(
                str(startdate), format="iso", scale="utc", out_subfmt="date"
            )
            next_days = [(startdate_astropy + i - 1).value for i in NIGHTS]

        # if startdate is None:
        #     now_astropy = Time(str(now), format="iso", scale="utc", out_subfmt="date")
        #     next_days = [(now_astropy + i - 1).value for i in NIGHTS]
        # else:
        #     startdate_astropy = Time(
        #         str(startdate), format="iso", scale="utc", out_subfmt="date"
        #     )
        #     next_days = [(startdate_astropy + i - 1).value for i in NIGHTS]

        # if self.ra is None:
        #     plan_initial = PlanObservation(
        #         name=name, date=str(today), alertsource="icecube"
        #     )
        # else:
        #     plan_initial = PlanObservation(
        #         name=name, date=str(today), ra=self.ra, dec=self.dec
        #     )

        # print(plan_initial.g_band_recommended_time_start)
        # quit()

        ra = plan_initial.ra
        dec = plan_initial.dec

        observable = []
        g_band_start = []
        g_band_end = []
        r_band_start = []
        r_band_end = []

        plan_initial.request_ztf_fields()

        if plan_initial.ra_err:
            recommended_field = plan_initial.recommended_field

        pdf_outfile = os.path.join(name, f"{name}_multiday.pdf")

        with PdfPages(pdf_outfile) as pdf:
            for i, day in enumerate(tqdm(next_days)):
                if NIGHTS[i] not in SHORT_NIGHTS:
                    plan = PlanObservation(
                        name=name, date=day, ra=ra, dec=dec, verbose=False
                    )
                else:
                    if NIGHTS[i] in ONE_FILTER_NIGHTS:
                        bands = ["g"]
                    else:
                        bands = ["g", "r"]
                    plan = PlanObservation(
                        name=name,
                        date=day,
                        ra=ra,
                        dec=dec,
                        observationlength=30,
                        bands=bands,
                        verbose=False,
                    )

                observable.append(plan.observable)

                if observable:
                    g_band_start.append(plan.g_band_recommended_time_start)
                    g_band_end.append(plan.g_band_recommended_time_end)
                    r_band_start.append(plan.r_band_recommended_time_start)
                    r_band_end.append(plan.r_band_recommended_time_end)
                else:
                    g_band_start.append(None)
                    g_band_end.append(None)
                    r_band_start.append(None)
                    r_band_end.append(None)

                ax = plan.plot_target()
                plt.tight_layout()
                pdf.savefig()
                plt.close()

        self.summarytext = f"\nYour multi-day observation plan for {name}\n"

        self.summarytext += "-------------------------------------------------\n"
        self.summarytext += "g-band observations\n"
        for i, item in enumerate(g_band_start):
            if item is not None:
                if observable[i]:
                    self.summarytext += f"Night {NIGHTS[i]} {short_time(item.value)} - {short_time(g_band_end[i].value)}\n"
                    exposure_time = isotime_delta_to_seconds(
                        isotime_start=item.value, isotime_end=g_band_end[i].value
                    )
                    self.triggers.append(
                        {
                            "field_id": recommended_field,
                            "filter_id": 1,
                            "mjd_start": isotime_to_mjd(item.value),
                            "exposure_time": exposure_time,
                        }
                    )
            else:
                self.summarytext += f"Night {NIGHTS[i]} NOT OBSERVABLE\n"
        self.summarytext += "-------------------------------------------------\n"

        self.summarytext += "\n-------------------------------------------------\n"
        self.summarytext += "r-band observations\n"

        for i, item in enumerate(r_band_start):
            if NIGHTS[i] not in ONE_FILTER_NIGHTS:
                if item is not None:
                    if observable[i]:
                        self.summarytext += f"Night {NIGHTS[i]} {short_time(item.value)} - {short_time(r_band_end[i].value)}\n"
                        exposure_time = isotime_delta_to_seconds(
                            isotime_start=item.value, isotime_end=r_band_end[i].value
                        )
                        self.triggers.append(
                            {
                                "field_id": recommended_field,
                                "filter_id": 2,
                                "mjd_start": isotime_to_mjd(item.value),
                                "exposure_time": exposure_time,
                            }
                        )

                else:
                    self.summarytext += f"Night {NIGHTS[i]} NOT OBSERVABLE\n"
        self.summarytext += "-------------------------------------------------\n\n"

    def print_plan(self):
        print(self.summarytext)

    def print_triggers(self):
        bands = {1: "g", 2: "r", 3: "i"}
        message = ""
        for i, trigger in enumerate(self.triggers):
            t_start = short_time(mjd_to_isotime(trigger["mjd_start"]))
            message += f"{t_start} // {trigger['exposure_time']} s exposure // filter={bands[trigger['filter_id']]} // field={trigger['field_id']}\n"
        message = message[:-1]
        print(message)
        return message
