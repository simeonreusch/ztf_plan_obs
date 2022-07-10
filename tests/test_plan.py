import unittest
import logging

from ztf_plan_obs.plan import PlanObservation
from ztf_plan_obs.multiday_plan import MultiDayObservation
from ztf_plan_obs.api import Queue


class TestPlan(unittest.TestCase):
    def setUp(self):

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        self.max_distance_diff_arcsec = 2

    def test_plan(self):

        self.logger.info("\n\n Testing Plan \n\n")

        neutrino_name = "IC220624A"
        date = "2022-06-24"

        self.logger.info(f"creatin an observation plan for neutrino {neutrino_name}")
        plan = PlanObservation(name=neutrino_name, date=date, alertsource="icecube")
        plan.plot_target()  # Plots the observing conditions
        plan.request_ztf_fields()

        recommended = plan.recommended_field
        expected_field = 720

        self.logger.info(f"recommended field: {recommended}, expected {expected_field}")
        self.assertEqual(recommended, expected_field)
