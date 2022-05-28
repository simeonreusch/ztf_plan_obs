import os, time
from typing import Union

from penquins import Kowalski


class APIError(Exception):
    pass


class Queue:
    """
    Submit observation triggers to Kowalski, query the queue and delete observation triggers
    """

    def __init__(
        self,
        user: str,
    ) -> None:

        self.user = user
        self.protocol: str = "https"
        self.host: str = os.environ.get("KOWALSKI_HOST", default="localhost")
        self.port: int = 443
        self.api_token: str = os.environ.get("KOWALSKI_API_TOKEN")

        self.queue: dict = {}

        if self.api_token is None:
            err = (
                "No kowalski API token found. Set the environment variable with \n"
                "export KOWALSKI_API_TOKEN=api_token"
            )
            raise APIError(err)

        self.kowalski = Kowalski(
            token=self.api_token, protocol=self.protocol, host=self.host, port=self.port
        )
        if not self.kowalski.ping():
            err = f"Ping of Kowalski with specified token failed. Are you sure this token is correct? Provided token: {self.api_token}"
            raise APIError(err)

    def get_all_queues(self, names_only: bool = False) -> Union[list, dict]:
        """
        Get all the queues
        """
        res = self.kowalski.api("get", "/api/triggers/ztf")
        if res["status"] != "success":
            err = f"API call failed with status '{res['status']}'' and message '{res['message']}''"
            raise APIError(err)

        if names_only:
            res = [x["queue_name"] for x in res["data"]]

        return res

    def get_too_queues(self, names_only: bool = False) -> Union[list, dict]:
        """
        Get all the queues and return ToO triggers only
        """
        res = self.get_all_queues()
        res["data"] = [x for x in res["data"] if x["is_TOO"]]

        if names_only:
            res = [x["queue_name"] for x in res["data"]]

        return res

    def add_trigger_to_queue(
        self,
        trigger_name: str,
        validity_window_start_mjd: float,
        field_id: list,
        filter_id: list,
        request_id: int = 1,
        subprogram_name: str = "ToO_Neutrino",
        exposure_time: int = [30],
        validity_window_end_mjd: float = None,
        program_id: int = 2,
        program_pi: str = "Kulkarni",
    ) -> None:
        """
        Add one trigger (requesting a single observation)
        to the queue (containing all the triggers that will be
        subbmitted)
        """
        if trigger_name[:4] != "ToO_":
            raise ValueError(
                f"Trigger names must begin with 'ToO_', but you entered '{trigger_name}'"
            )

        if validity_window_end_mjd is None:
            validity_window_end_mjd = validity_window_start_mjd + exposure_time / 86400

        targets = [
            {
                "request_id": request_id,
                "field_id": field_id,
                "filter_id": filter_id,
                "subprogram_name": subprogram_name,
                "program_pi": program_pi,
                "program_id": program_id,
                "exposure_time": exposure_time,
            }
        ]

        trigger_id = len(self.queue)

        trigger = {
            trigger_id: {
                "user": self.user,
                "queue_name": f"{trigger_name}_{trigger_id}",
                "queue_type": "list",
                "validity_window_mjd": [
                    validity_window_start_mjd,
                    validity_window_end_mjd,
                ],
                "targets": targets,
            }
        }
        self.queue.update(trigger)

    def submit_queue(self) -> None:
        """
        Submit the queue of triggers via the Kowalski API
        """
        results = []
        for i, trigger in self.queue.items():
            res = self.kowalski.api(
                method="put", endpoint="/api/triggers/ztf", data=trigger
            )
            results.append(res)
            if res["status"] != "success":
                err = "something went wrong with submitting."
                raise APIError(err)

        print(f"Submitted {len(self.queue)} triggers to Kowalski.")

        return results

    def delete_trigger(self, trigger_name) -> None:
        """
        Delete a trigger that has been submitted
        """
        req = {"user": self.user, "queue_name": trigger_name}

        res = self.kowalski.api(method="delete", endpoint="/api/triggers/ztf", data=req)

        if res["status"] != "success":
            err = "something went wrong with deleting the trigger."
            raise APIError(err)

        return res
