import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT_PATH = (
    Path(__file__).resolve().parents[3] / "script" / "ci" / "run-gitlab-extended-ci.py"
)
SPEC = importlib.util.spec_from_file_location("run_gitlab_extended_ci", SCRIPT_PATH)
assert SPEC is not None
assert SPEC.loader is not None
BRIDGE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BRIDGE)

SOURCE_SHA = "a" * 40
GITHUB_RUN_URL = "https://github.com/example/haseongpu/actions/runs/123"
GITLAB_PIPELINE_URL = (
    "https://codebase.helmholtz.cloud/th168408/"
    "haseongpu-extended-ci/-/pipelines/42"
)


def configure_environment(monkeypatch):
    environment = {
        "HASE_SOURCE_SHA": SOURCE_SHA,
        "HASE_SOURCE_REPOSITORY": "contributor/haseongpu",
        "HASE_GITHUB_REPOSITORY": "example/haseongpu",
        "HASE_GITHUB_TOKEN": "github-token",
        "HASE_GITHUB_RUN_URL": GITHUB_RUN_URL,
        "HASE_GITLAB_PROJECT_ID": "1234",
        "HASE_GITLAB_TRIGGER_TOKEN": "trigger-token",
        "HASE_GITLAB_READ_API_TOKEN": "read-token",
        "HASE_SOURCE_EVENT": "pull_request",
        "HASE_SOURCE_REF": "feature",
        "HASE_PR_NUMBER": "206",
    }
    for name, value in environment.items():
        monkeypatch.setenv(name, value)


def record_statuses(monkeypatch):
    statuses = []

    def record(repository, source_sha, token, state, description, target_url):
        statuses.append(
            {
                "repository": repository,
                "source_sha": source_sha,
                "token": token,
                "state": state,
                "description": description,
                "target_url": target_url,
            }
        )

    monkeypatch.setattr(BRIDGE, "post_commit_status", record)
    return statuses


def configure_pipeline(monkeypatch):
    configure_environment(monkeypatch)
    triggered = {}

    def trigger(project_id, pipeline_ref, trigger_token, variables):
        triggered.update(
            {
                "project_id": project_id,
                "pipeline_ref": pipeline_ref,
                "trigger_token": trigger_token,
                "variables": variables,
            }
        )
        return {
            "id": 42,
            "project_id": 1234,
            "web_url": GITLAB_PIPELINE_URL,
        }

    monkeypatch.setattr(
        BRIDGE,
        "trigger_pipeline",
        trigger,
    )
    monkeypatch.setattr(BRIDGE.time, "sleep", lambda seconds: None)
    monkeypatch.setattr(BRIDGE, "append_step_summary", lambda *args: None)
    return triggered


def test_post_commit_status_uses_github_status_api(monkeypatch):
    captured = {}

    def capture(request, *, attempts=1, service="GitLab"):
        captured["request"] = request
        captured["attempts"] = attempts
        captured["service"] = service
        return {}

    monkeypatch.setattr(BRIDGE, "request_json", capture)

    BRIDGE.post_commit_status(
        "example/haseongpu",
        SOURCE_SHA,
        "github-token",
        "pending",
        "GitLab pipeline is pending",
        GITLAB_PIPELINE_URL,
    )

    request = captured["request"]
    assert request.full_url == (
        f"https://api.github.com/repos/example/haseongpu/statuses/{SOURCE_SHA}"
    )
    assert request.method == "POST"
    assert request.get_header("Authorization") == "Bearer github-token"
    assert request.get_header("X-github-api-version") == "2026-03-10"
    assert captured["attempts"] == 3
    assert captured["service"] == "GitHub"
    assert json.loads(request.data) == {
        "state": "pending",
        "context": "ci/gitlab/codebase.helmholtz.cloud",
        "description": "GitLab pipeline is pending",
        "target_url": GITLAB_PIPELINE_URL,
    }


def test_success_updates_pending_status_and_pipeline_result(monkeypatch):
    statuses = record_statuses(monkeypatch)
    triggered = configure_pipeline(monkeypatch)
    pipelines = iter([{"status": "running"}, {"status": "success"}])
    monkeypatch.setattr(
        BRIDGE,
        "get_pipeline",
        lambda project_id, pipeline_id, read_token: next(pipelines),
    )

    assert BRIDGE.main() == 0

    assert triggered["project_id"] == "1234"
    assert triggered["pipeline_ref"] == "master"
    assert triggered["trigger_token"] == "trigger-token"
    assert triggered["variables"] == {
        "HASE_SOURCE_SHA": SOURCE_SHA,
        "HASE_SOURCE_REPOSITORY": "contributor/haseongpu",
        "HASE_SOURCE_EVENT": "pull_request",
        "HASE_SOURCE_REF": "feature",
        "HASE_PR_NUMBER": "206",
        "HASE_GITHUB_RUN_URL": GITHUB_RUN_URL,
    }
    assert [status["state"] for status in statuses] == ["pending", "success"]
    assert statuses[0]["target_url"] == GITLAB_PIPELINE_URL
    assert statuses[1]["target_url"] == GITLAB_PIPELINE_URL


def test_failed_pipeline_updates_failure_status(monkeypatch):
    statuses = record_statuses(monkeypatch)
    configure_pipeline(monkeypatch)
    monkeypatch.setattr(
        BRIDGE,
        "get_pipeline",
        lambda project_id, pipeline_id, read_token: {"status": "failed"},
    )

    assert BRIDGE.main() == 1

    assert [status["state"] for status in statuses] == ["pending", "failure"]
    assert statuses[-1]["target_url"] == GITLAB_PIPELINE_URL


def test_timeout_updates_error_status(monkeypatch):
    statuses = record_statuses(monkeypatch)
    configure_pipeline(monkeypatch)
    monkeypatch.setattr(BRIDGE, "TIMEOUT_SECONDS", 0)
    monkeypatch.setattr(
        BRIDGE,
        "get_pipeline",
        lambda project_id, pipeline_id, read_token: {"status": "running"},
    )

    assert BRIDGE.main() == 1

    assert [status["state"] for status in statuses] == ["pending", "error"]
    assert "Timed out" in statuses[-1]["description"]
    assert statuses[-1]["target_url"] == GITLAB_PIPELINE_URL


def test_trigger_error_updates_error_status(monkeypatch):
    configure_environment(monkeypatch)
    statuses = record_statuses(monkeypatch)

    def fail_trigger(project_id, pipeline_ref, trigger_token, variables):
        raise RuntimeError("trigger failed")

    monkeypatch.setattr(BRIDGE, "trigger_pipeline", fail_trigger)

    with pytest.raises(RuntimeError, match="trigger failed"):
        BRIDGE.main()

    assert [status["state"] for status in statuses] == ["error"]
    assert statuses[-1]["target_url"] == GITHUB_RUN_URL


@pytest.mark.parametrize(
    "repository",
    [
        "https://github.com/contributor/haseongpu.git",
        "contributor",
        "contributor/../haseongpu",
        "contributor/..",
        "contributor/haseongpu extra",
        "-contributor/haseongpu",
        "contributor-/haseongpu",
    ],
)
def test_invalid_source_repository_is_rejected(monkeypatch, repository):
    configure_environment(monkeypatch)
    monkeypatch.setenv("HASE_SOURCE_REPOSITORY", repository)
    statuses = record_statuses(monkeypatch)

    def unexpected_trigger(*args):
        raise AssertionError("invalid source repository reached the GitLab trigger")

    monkeypatch.setattr(BRIDGE, "trigger_pipeline", unexpected_trigger)

    with pytest.raises(RuntimeError, match="HASE_SOURCE_REPOSITORY"):
        BRIDGE.main()

    assert [status["state"] for status in statuses] == ["error"]
    assert statuses[0]["target_url"] == GITHUB_RUN_URL


@pytest.mark.parametrize(
    "pipeline",
    [
        {"id": 42, "project_id": 1234},
        {
            "id": 42,
            "project_id": 1234,
            "web_url": "https://codebase.example/pipelines/42",
        },
        {
            "id": 42,
            "project_id": 4321,
            "web_url": (
                "https://codebase.helmholtz.cloud/th168408/"
                "haseongpu-extended-ci/-/pipelines/42"
            ),
        },
    ],
)
def test_invalid_pipeline_url_updates_error_status(monkeypatch, pipeline):
    configure_environment(monkeypatch)
    statuses = record_statuses(monkeypatch)
    monkeypatch.setattr(
        BRIDGE,
        "trigger_pipeline",
        lambda project_id, pipeline_ref, trigger_token, variables: pipeline,
    )

    with pytest.raises(RuntimeError, match="expected"):
        BRIDGE.main()

    assert [status["state"] for status in statuses] == ["error"]
    assert statuses[0]["target_url"] == GITHUB_RUN_URL
