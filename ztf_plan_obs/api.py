from penquins import Kowalski
import os


class APIError(Exception):
    pass


protocol = "https"
host = os.environ.get("KOWALSKI_HOST", default="localhost")
port = 443

api_token = os.environ.get("KOWALSKI_API_TOKEN")

if api_token is None:
    err = "No kowalski API token found. Set the envirnoment variable with \n" \
          " export KOWALSKI_API_TOKEN=api_token"
    raise APIError(err)

kowalski = Kowalski(
    token=api_token,
    protocol=protocol,
    host=host,
    port=port
)
if not kowalski.ping():
    raise APIError("Ping of kowalski with specified token failed. Are you sure this token is correct?")


def get_all_queues():
    res = kowalski.api("get", "/api/triggers/ztf")
    if res["status"] != "success":
        raise APIError(f"API call failed with status '{res['status']}'' and message '{res['message']}''")
    return res


def get_too_queues():
    res = get_all_queues()
    res["data"] = [x for x in res["data"] if x["is_TOO"]]
    return res


def build_queue_entry(
        request_id: int,
        field_id: int,
        filter_id: int,
        subprogram_name: str,
        program_pi: str = "Kulkarni",
        program_id: int = 2,
        exposure_time: float = None,
        ra: float = None,
        dec: float = None,
        n_repeats: int = None,
        max_airmass: float = None

) -> dict:
    args = locals()

    entry = {}

    for key, val in args.items():
        if val is not None:
            entry[key] = val

    if subprogram_name[:4] != "ToO_":
        raise ValueError(f"Queue subprogram names must begin with 'ToO_', but you entered '{subprogram_name}'")

    return entry


def build_request(
        user: str,
        queue_name: str,
        validity_window_start_mjd: float,
        validity_window_end_mjd: float,
        subprogram_name: str,
        field_ids: list,
        filter_ids: list,
        exposure_times: list = None,
        program_id: int = 2,
        program_pi: str = "Kulkarni"
):
    if queue_name[:4] != "ToO_":
        raise ValueError(f"Queue names must begin with 'ToO_', but you entered '{queue_name}'")

    targets = []

    for i, field_id in enumerate(field_ids):

        if exposure_times is None:
            exp_time = None
        else:
            exp_time = exposure_times[i]

        target = build_queue_entry(
            request_id=i,
            field_id=field_id,
            filter_id=filter_ids[i],
            subprogram_name=subprogram_name,
            program_pi=program_pi,
            program_id=program_id,
            exposure_time=exp_time
        )
        targets.append(target)

    payload = {
        "user": user,
        "queue_name": queue_name,
        "queue_type": "list",
        "validity_window_mjd": [
            validity_window_start_mjd,
            validity_window_end_mjd
        ],
        "targets": targets
    }
    return payload


def submit_request(payload):
    return kowalski.api(method="put", endpoint="/api/triggers/ztf", data=payload)


def delete_request(
        user,
        queue_name
):
    req = {
        "user": user,
        "queue_name": queue_name
    }

    return kowalski.api(method="delete", endpoint="/api/triggers/ztf", data=req)




