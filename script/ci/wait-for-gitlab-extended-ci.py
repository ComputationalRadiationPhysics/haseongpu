#!/usr/bin/env python3

"""Mirror the GitLab bridge commit status into the normal CI workflow."""

from __future__ import annotations

import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable


GITHUB_API_VERSION = "2026-03-10"
GITHUB_STATUS_CONTEXT = "ci/gitlab/codebase.helmholtz.cloud"
GITHUB_REPOSITORY_PATTERN = re.compile(
    r"[A-Za-z0-9](?:[A-Za-z0-9-]{0,37}[A-Za-z0-9])?/"
    r"[A-Za-z0-9._-]{1,100}"
)
RETRY_DELAY_SECONDS = 5
POLL_INTERVAL_SECONDS = 15
# Leave time to report an error before the 360-minute workflow timeout.
TIMEOUT_SECONDS = 355 * 60


def required_environment(name: str) -> str:
    value = os.environ.get(name, "").strip()
    if not value:
        raise RuntimeError(f"required environment variable is empty: {name}")
    return value


def positive_integer_environment(name: str) -> int:
    value = required_environment(name)
    try:
        parsed = int(value)
    except ValueError as error:
        raise RuntimeError(f"{name} must be a positive integer") from error
    if parsed < 1:
        raise RuntimeError(f"{name} must be a positive integer")
    return parsed


def github_repository_environment(name: str) -> str:
    value = required_environment(name)
    if not GITHUB_REPOSITORY_PATTERN.fullmatch(value):
        raise RuntimeError(f"{name} must have the form owner/repository")
    _, repository = value.split("/", maxsplit=1)
    if repository in {".", ".."}:
        raise RuntimeError(f"{name} must have the form owner/repository")
    return value


def status_marker(source_run_id: int, source_run_attempt: int) -> str:
    return f"[CI run {source_run_id}/{source_run_attempt}]"


def request_statuses(
    api_url: str,
    repository: str,
    source_sha: str,
    token: str,
    *,
    attempts: int = 3,
) -> list[dict[str, object]]:
    encoded_sha = urllib.parse.quote(source_sha, safe="")
    request = urllib.request.Request(
        f"{api_url}/repos/{repository}/commits/{encoded_sha}/statuses?per_page=100",
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "X-GitHub-Api-Version": GITHUB_API_VERSION,
        },
    )
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = json.load(response)
            if not isinstance(payload, list):
                raise RuntimeError("GitHub statuses API returned a non-list response")
            return [status for status in payload if isinstance(status, dict)]
        except urllib.error.HTTPError as error:
            raise RuntimeError(
                f"GitHub API request failed: HTTP {error.code} {error.reason}"
            ) from error
        except (urllib.error.URLError, TimeoutError) as error:
            if attempt == attempts:
                raise RuntimeError(f"GitHub API request failed: {error}") from error
            time.sleep(RETRY_DELAY_SECONDS * attempt)
        except (json.JSONDecodeError, UnicodeDecodeError) as error:
            raise RuntimeError("GitHub API returned invalid JSON") from error
    raise AssertionError("unreachable")


def select_status(
    statuses: list[dict[str, object]], marker: str
) -> dict[str, object] | None:
    matches = [
        status
        for status in statuses
        if status.get("context") == GITHUB_STATUS_CONTEXT
        and str(status.get("description", "")).endswith(marker)
    ]
    if not matches:
        return None
    return max(matches, key=lambda status: int(status.get("id", 0)))


def wait_for_status(
    fetch_statuses: Callable[[], list[dict[str, object]]],
    marker: str,
    *,
    timeout_seconds: float = TIMEOUT_SECONDS,
    poll_interval_seconds: float = POLL_INTERVAL_SECONDS,
    monotonic: Callable[[], float] = time.monotonic,
    sleep: Callable[[float], None] = time.sleep,
) -> int:
    deadline = monotonic() + timeout_seconds
    previous_state = ""

    while True:
        status = select_status(fetch_statuses(), marker)
        state = "not reported" if status is None else str(status.get("state", ""))
        if state != previous_state:
            print(f"GitLab extended CI status: {state}", flush=True)
            if status is not None and status.get("target_url"):
                print(f"GitLab pipeline: {status['target_url']}", flush=True)
            previous_state = state

        if state == "success":
            return 0
        if state in {"error", "failure"}:
            description = str(status.get("description", ""))
            print(
                f"GitLab extended CI did not succeed: {description or state}",
                file=sys.stderr,
            )
            return 1
        if state not in {"not reported", "pending"}:
            print(f"Unexpected GitHub commit status state: {state}", file=sys.stderr)
            return 1
        if monotonic() >= deadline:
            print("Timed out waiting for GitLab extended CI status", file=sys.stderr)
            return 1
        sleep(poll_interval_seconds)


def main() -> int:
    api_url = required_environment("GITHUB_API_URL").rstrip("/")
    repository = github_repository_environment("GITHUB_REPOSITORY")
    token = required_environment("HASE_GITHUB_TOKEN")
    source_sha = required_environment("HASE_SOURCE_SHA").lower()
    if not re.fullmatch(r"[0-9a-f]{40}", source_sha):
        raise RuntimeError("HASE_SOURCE_SHA must be a full hexadecimal Git commit")
    source_run_id = positive_integer_environment("HASE_SOURCE_RUN_ID")
    source_run_attempt = positive_integer_environment("HASE_SOURCE_RUN_ATTEMPT")
    marker = status_marker(source_run_id, source_run_attempt)

    print(f"Waiting for GitLab extended CI status {marker}", flush=True)
    return wait_for_status(
        lambda: request_statuses(api_url, repository, source_sha, token),
        marker,
    )


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(error, file=sys.stderr)
        raise SystemExit(1) from error
