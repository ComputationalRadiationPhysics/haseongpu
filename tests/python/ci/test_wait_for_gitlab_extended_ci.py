import importlib.util
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[3]
    / "script"
    / "ci"
    / "wait-for-gitlab-extended-ci.py"
)
SPEC = importlib.util.spec_from_file_location("wait_for_gitlab_extended_ci", SCRIPT_PATH)
assert SPEC is not None
assert SPEC.loader is not None
WAITER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(WAITER)


def test_select_status_matches_source_attempt_and_uses_latest_status():
    marker = WAITER.status_marker(123, 2)
    selected = WAITER.select_status(
        [
            {
                "id": 9,
                "context": WAITER.GITHUB_STATUS_CONTEXT,
                "description": "Succeeded [CI run 123/1]",
            },
            {
                "id": 10,
                "context": WAITER.GITHUB_STATUS_CONTEXT,
                "description": f"Pending {marker}",
            },
            {
                "id": 11,
                "context": WAITER.GITHUB_STATUS_CONTEXT,
                "description": f"Succeeded {marker}",
            },
            {
                "id": 12,
                "context": "unrelated",
                "description": f"Succeeded {marker}",
            },
        ],
        marker,
    )

    assert selected["id"] == 11


def test_wait_for_status_reports_success_after_pending():
    marker = WAITER.status_marker(123, 1)
    responses = iter(
        [
            [],
            [
                {
                    "id": 20,
                    "context": WAITER.GITHUB_STATUS_CONTEXT,
                    "description": f"In progress {marker}",
                    "state": "pending",
                }
            ],
            [
                {
                    "id": 21,
                    "context": WAITER.GITHUB_STATUS_CONTEXT,
                    "description": f"Succeeded {marker}",
                    "state": "success",
                }
            ],
        ]
    )

    result = WAITER.wait_for_status(
        lambda: next(responses),
        marker,
        poll_interval_seconds=0,
        sleep=lambda _: None,
    )

    assert result == 0


def test_wait_for_status_propagates_failure():
    marker = WAITER.status_marker(123, 1)

    result = WAITER.wait_for_status(
        lambda: [
            {
                "id": 20,
                "context": WAITER.GITHUB_STATUS_CONTEXT,
                "description": f"Failed {marker}",
                "state": "failure",
            }
        ],
        marker,
        poll_interval_seconds=0,
        sleep=lambda _: None,
    )

    assert result == 1


def test_wait_for_status_times_out_when_status_is_not_reported():
    times = iter([0.0, 1.0])

    result = WAITER.wait_for_status(
        lambda: [],
        WAITER.status_marker(123, 1),
        timeout_seconds=1,
        poll_interval_seconds=0,
        monotonic=lambda: next(times),
        sleep=lambda _: None,
    )

    assert result == 1
