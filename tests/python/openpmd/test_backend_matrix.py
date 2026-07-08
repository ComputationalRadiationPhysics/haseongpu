import pytest

import openpmd_backend_matrix


def test_openpmd_backend_matrix_uses_probe_artifact_by_default(monkeypatch, tmp_path):
    monkeypatch.delenv("HASE_OPENPMD_TEST_BACKENDS", raising=False)
    (tmp_path / "openpmd_backends.txt").write_text("hdf5\n", encoding="utf-8")
    monkeypatch.setattr(openpmd_backend_matrix, "_binding_package_dirs", lambda: (tmp_path,))
    monkeypatch.setattr(openpmd_backend_matrix, "_build_artifact_candidates", lambda: ())

    assert openpmd_backend_matrix.openpmd_test_backends() == ["hdf5"]


def test_openpmd_backend_matrix_accepts_explicit_plural_selector(monkeypatch):
    monkeypatch.setenv("HASE_OPENPMD_TEST_BACKENDS", "hdf5, adios-sst")

    assert openpmd_backend_matrix.openpmd_test_backends() == ["adios-sst", "hdf5"]


def test_openpmd_backend_matrix_fails_without_probe_artifact(monkeypatch):
    monkeypatch.delenv("HASE_OPENPMD_TEST_BACKENDS", raising=False)
    monkeypatch.setattr(openpmd_backend_matrix, "_binding_package_dirs", lambda: ())
    monkeypatch.setattr(openpmd_backend_matrix, "_build_artifact_candidates", lambda: ())

    with pytest.raises(RuntimeError, match="No openPMD backend probe artifact found"):
        openpmd_backend_matrix.openpmd_test_backends()


def test_openpmd_backend_matrix_fails_on_empty_probe_artifact(monkeypatch, tmp_path):
    monkeypatch.delenv("HASE_OPENPMD_TEST_BACKENDS", raising=False)
    (tmp_path / "openpmd_backends.txt").write_text("", encoding="utf-8")
    monkeypatch.setattr(openpmd_backend_matrix, "_binding_package_dirs", lambda: (tmp_path,))
    monkeypatch.setattr(openpmd_backend_matrix, "_build_artifact_candidates", lambda: ())

    with pytest.raises(RuntimeError, match="probe artifact is empty"):
        openpmd_backend_matrix.openpmd_test_backends()
