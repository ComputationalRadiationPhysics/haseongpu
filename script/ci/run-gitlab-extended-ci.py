#!/usr/bin/env python3

"""Trigger HASEonGPU extended CI on GitLab and wait for its result."""

from __future__ import annotations

import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


GITLAB_API_URL = "https://codebase.helmholtz.cloud/api/v4"
GITLAB_PIPELINE_REF = "master"
POLL_INTERVAL_SECONDS = 15
TIMEOUT_SECONDS = 6 * 60 * 60

TRANSIENT_STATUSES = {
    "created",
    "waiting_for_resource",
    "preparing",
    "waiting_for_callback",
    "pending",
    "running",
    "scheduled",
}


def required_environment(name: str) -> str:
    value = os.environ.get(name, "").strip()
    if not value:
        raise RuntimeError(f"required environment variable is empty: {name}")
    return value


def request_json(
    request: urllib.request.Request,
    *,
    attempts: int = 1,
) -> dict[str, object]:
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.load(response)
        except urllib.error.HTTPError as error:
            raise RuntimeError(
                f"GitLab API request failed: HTTP {error.code} {error.reason}"
            ) from error
        except (urllib.error.URLError, TimeoutError) as error:
            if attempt == attempts:
                raise RuntimeError(f"GitLab API request failed: {error}") from error
            time.sleep(min(5 * attempt, POLL_INTERVAL_SECONDS))
    raise AssertionError("unreachable")


def trigger_pipeline(
    project_id: str,
    pipeline_ref: str,
    trigger_token: str,
    variables: dict[str, str],
) -> dict[str, object]:
    fields = {
        "token": trigger_token,
        "ref": pipeline_ref,
        **{f"variables[{name}]": value for name, value in variables.items()},
    }
    request = urllib.request.Request(
        f"{GITLAB_API_URL}/projects/{project_id}/trigger/pipeline",
        data=urllib.parse.urlencode(fields).encode("utf-8"),
        method="POST",
    )
    return request_json(request)


def get_pipeline(
    project_id: str, pipeline_id: int, read_token: str
) -> dict[str, object]:
    request = urllib.request.Request(
        f"{GITLAB_API_URL}/projects/{project_id}/pipelines/{pipeline_id}",
        headers={"PRIVATE-TOKEN": read_token},
    )
    return request_json(request, attempts=3)


def append_step_summary(pipeline_id: int, source_sha: str) -> None:
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary_path:
        return
    with Path(summary_path).open("a", encoding="utf-8") as summary:
        summary.write("## GitLab extended CI\n\n")
        summary.write(f"- Source commit: `{source_sha}`\n")
        summary.write(f"- Pipeline ID: `{pipeline_id}`\n")


def main() -> int:
    project_id = required_environment("HASE_GITLAB_PROJECT_ID")
    if not re.fullmatch(r"[1-9][0-9]*", project_id):
        raise RuntimeError("HASE_GITLAB_PROJECT_ID must be a positive integer")
    pipeline_ref = os.environ.get("HASE_GITLAB_PIPELINE_REF", GITLAB_PIPELINE_REF)
    trigger_token = required_environment("HASE_GITLAB_TRIGGER_TOKEN")
    read_token = required_environment("HASE_GITLAB_READ_API_TOKEN")
    source_sha = required_environment("HASE_SOURCE_SHA").lower()
    if not re.fullmatch(r"[0-9a-f]{40}", source_sha):
        raise RuntimeError("HASE_SOURCE_SHA must be a full hexadecimal Git commit")

    variables = {
        "HASE_SOURCE_SHA": source_sha,
        "HASE_SOURCE_EVENT": os.environ.get("HASE_SOURCE_EVENT", "unknown"),
        "HASE_SOURCE_REF": os.environ.get("HASE_SOURCE_REF", "unknown"),
        "HASE_PR_NUMBER": os.environ.get("HASE_PR_NUMBER", ""),
        "HASE_GITHUB_RUN_URL": os.environ.get("HASE_GITHUB_RUN_URL", ""),
    }
    pipeline = trigger_pipeline(project_id, pipeline_ref, trigger_token, variables)
    try:
        pipeline_id = int(pipeline["id"])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError("GitLab trigger response did not identify a pipeline") from error

    print(f"Triggered GitLab pipeline {pipeline_id}", flush=True)
    append_step_summary(pipeline_id, source_sha)

    deadline = time.monotonic() + TIMEOUT_SECONDS
    previous_status = ""
    while True:
        pipeline = get_pipeline(project_id, pipeline_id, read_token)
        status = str(pipeline.get("status", ""))
        if status != previous_status:
            print(f"GitLab pipeline {pipeline_id}: {status}", flush=True)
            previous_status = status

        if status == "success":
            return 0
        if status not in TRANSIENT_STATUSES:
            print(f"GitLab pipeline {pipeline_id} did not succeed", file=sys.stderr)
            return 1
        if time.monotonic() >= deadline:
            print(f"Timed out waiting for GitLab pipeline {pipeline_id}", file=sys.stderr)
            return 1
        time.sleep(POLL_INTERVAL_SECONDS)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(error, file=sys.stderr)
        raise SystemExit(1) from error
