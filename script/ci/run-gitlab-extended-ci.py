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
GITLAB_WEB_HOST = "codebase.helmholtz.cloud"
GITHUB_API_URL = "https://api.github.com"
GITHUB_API_VERSION = "2026-03-10"
GITHUB_STATUS_CONTEXT = "ci/gitlab/codebase.helmholtz.cloud"
GITHUB_STATUS_DESCRIPTION_LIMIT = 140
GITHUB_REPOSITORY_PATTERN = re.compile(
    r"[A-Za-z0-9](?:[A-Za-z0-9-]{0,37}[A-Za-z0-9])?/"
    r"[A-Za-z0-9._-]{1,100}"
)
RETRY_DELAY_SECONDS = 5
POLL_INTERVAL_SECONDS = 15
# Leave time to publish an error before the 360-minute workflow timeout.
TIMEOUT_SECONDS = 355 * 60

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


def positive_integer_environment(name: str) -> int:
    value = required_environment(name)
    try:
        parsed = int(value)
    except ValueError as error:
        raise RuntimeError(f"{name} must be a positive integer") from error
    if parsed < 1:
        raise RuntimeError(f"{name} must be a positive integer")
    return parsed


def github_status_description(
    message: str, source_run_id: int, source_run_attempt: int
) -> str:
    suffix = f" [CI run {source_run_id}/{source_run_attempt}]"
    return f"{message[: GITHUB_STATUS_DESCRIPTION_LIMIT - len(suffix)]}{suffix}"


def github_repository_environment(name: str) -> str:
    value = required_environment(name)
    if not GITHUB_REPOSITORY_PATTERN.fullmatch(value):
        raise RuntimeError(f"{name} must have the form owner/repository")
    _, repository = value.split("/", maxsplit=1)
    if repository in {".", ".."}:
        raise RuntimeError(f"{name} must have the form owner/repository")
    return value


def request_json(
    request: urllib.request.Request,
    *,
    attempts: int = 1,
    service: str = "GitLab",
) -> dict[str, object]:
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise RuntimeError(f"{service} API returned a non-object response")
            return payload
        except urllib.error.HTTPError as error:
            raise RuntimeError(
                f"{service} API request failed: HTTP {error.code} {error.reason}"
            ) from error
        except (urllib.error.URLError, TimeoutError) as error:
            if attempt == attempts:
                raise RuntimeError(f"{service} API request failed: {error}") from error
            time.sleep(RETRY_DELAY_SECONDS * attempt)
        except (json.JSONDecodeError, UnicodeDecodeError) as error:
            raise RuntimeError(f"{service} API returned invalid JSON") from error
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


def post_commit_status(
    repository: str,
    source_sha: str,
    github_token: str,
    state: str,
    description: str,
    target_url: str,
) -> None:
    payload = {
        "state": state,
        "context": GITHUB_STATUS_CONTEXT,
        "description": description,
        "target_url": target_url,
    }
    request = urllib.request.Request(
        f"{GITHUB_API_URL}/repos/{repository}/statuses/{source_sha}",
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {github_token}",
            "Content-Type": "application/json",
            "X-GitHub-Api-Version": GITHUB_API_VERSION,
        },
        method="POST",
    )
    request_json(request, attempts=3, service="GitHub")


def append_step_summary(
    pipeline_id: int, source_sha: str, pipeline_url: str
) -> None:
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary_path:
        return
    with Path(summary_path).open("a", encoding="utf-8") as summary:
        summary.write("## GitLab extended CI\n\n")
        summary.write(f"- Source commit: `{source_sha}`\n")
        summary.write(f"- Pipeline: [#{pipeline_id}]({pipeline_url})\n")


def pipeline_web_url(
    pipeline: dict[str, object], pipeline_id: int, project_id: str
) -> str:
    try:
        response_project_id = int(pipeline["project_id"])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError(
            "GitLab trigger response did not identify the expected project"
        ) from error
    if response_project_id != int(project_id):
        raise RuntimeError(
            "GitLab trigger response did not identify the expected project"
        )

    value = pipeline.get("web_url")
    if isinstance(value, str):
        parsed = urllib.parse.urlparse(value)
        if (
            parsed.scheme == "https"
            and parsed.netloc == GITLAB_WEB_HOST
            and parsed.path.endswith(f"/-/pipelines/{pipeline_id}")
            and not parsed.params
            and not parsed.query
            and not parsed.fragment
        ):
            return value
    raise RuntimeError(
        "GitLab trigger response did not include the expected pipeline URL"
    )


def main() -> int:
    source_sha = required_environment("HASE_SOURCE_SHA").lower()
    if not re.fullmatch(r"[0-9a-f]{40}", source_sha):
        raise RuntimeError("HASE_SOURCE_SHA must be a full hexadecimal Git commit")
    source_run_id = positive_integer_environment("HASE_SOURCE_RUN_ID")
    source_run_attempt = positive_integer_environment("HASE_SOURCE_RUN_ATTEMPT")
    github_repository = github_repository_environment("HASE_GITHUB_REPOSITORY")
    github_token = required_environment("HASE_GITHUB_TOKEN")
    github_run_url = required_environment("HASE_GITHUB_RUN_URL")
    github_bridge_run_url = required_environment("HASE_GITHUB_BRIDGE_RUN_URL")
    pipeline_url = github_bridge_run_url

    def publish_status(state: str, description: str, target_url: str) -> None:
        post_commit_status(
            github_repository,
            source_sha,
            github_token,
            state,
            github_status_description(
                description, source_run_id, source_run_attempt
            ),
            target_url,
        )

    try:
        project_id = required_environment("HASE_GITLAB_PROJECT_ID")
        if not re.fullmatch(r"[1-9][0-9]*", project_id):
            raise RuntimeError("HASE_GITLAB_PROJECT_ID must be a positive integer")
        pipeline_ref = os.environ.get(
            "HASE_GITLAB_PIPELINE_REF", GITLAB_PIPELINE_REF
        )
        trigger_token = required_environment("HASE_GITLAB_TRIGGER_TOKEN")
        read_token = required_environment("HASE_GITLAB_READ_API_TOKEN")
        source_repository = github_repository_environment("HASE_SOURCE_REPOSITORY")

        variables = {
            "HASE_SOURCE_SHA": source_sha,
            "HASE_SOURCE_REPOSITORY": source_repository,
            "HASE_SOURCE_EVENT": os.environ.get("HASE_SOURCE_EVENT", "unknown"),
            "HASE_SOURCE_REF": os.environ.get("HASE_SOURCE_REF", "unknown"),
            "HASE_PR_NUMBER": os.environ.get("HASE_PR_NUMBER", ""),
            "HASE_GITHUB_RUN_URL": github_run_url,
        }
        pipeline = trigger_pipeline(project_id, pipeline_ref, trigger_token, variables)
        try:
            pipeline_id = int(pipeline["id"])
        except (KeyError, TypeError, ValueError) as error:
            raise RuntimeError(
                "GitLab trigger response did not identify a pipeline"
            ) from error
        if pipeline_id <= 0:
            raise RuntimeError("GitLab trigger response did not identify a pipeline")

        pipeline_url = pipeline_web_url(pipeline, pipeline_id, project_id)
        print(f"Triggered GitLab pipeline {pipeline_id}", flush=True)
        append_step_summary(pipeline_id, source_sha, pipeline_url)
        publish_status(
            "pending",
            f"GitLab pipeline {pipeline_id} is in progress",
            pipeline_url,
        )

        deadline = time.monotonic() + TIMEOUT_SECONDS
        previous_status = ""
        while True:
            pipeline = get_pipeline(project_id, pipeline_id, read_token)
            status = str(pipeline.get("status", ""))
            if status != previous_status:
                print(f"GitLab pipeline {pipeline_id}: {status}", flush=True)
                previous_status = status

            if status == "success":
                publish_status(
                    "success",
                    f"GitLab pipeline {pipeline_id} succeeded",
                    pipeline_url,
                )
                return 0
            if status not in TRANSIENT_STATUSES:
                publish_status(
                    "failure",
                    f"GitLab pipeline {pipeline_id} finished with status {status}",
                    pipeline_url,
                )
                print(
                    f"GitLab pipeline {pipeline_id} did not succeed", file=sys.stderr
                )
                return 1
            if time.monotonic() >= deadline:
                publish_status(
                    "error",
                    f"Timed out waiting for GitLab pipeline {pipeline_id}",
                    pipeline_url,
                )
                print(
                    f"Timed out waiting for GitLab pipeline {pipeline_id}",
                    file=sys.stderr,
                )
                return 1
            time.sleep(POLL_INTERVAL_SECONDS)
    except RuntimeError:
        try:
            publish_status(
                "error",
                "GitLab bridge failed before the pipeline completed",
                pipeline_url,
            )
        except RuntimeError as status_error:
            print(
                f"Additionally failed to update GitHub status: {status_error}",
                file=sys.stderr,
            )
        raise


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(error, file=sys.stderr)
        raise SystemExit(1) from error
