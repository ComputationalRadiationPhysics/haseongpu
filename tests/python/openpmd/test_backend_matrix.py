import pytest

import openpmd_backend_matrix


def test_backendMatrixUsesFrontEndBackends(monkeypatch):
    monkeypatch.delenv("HASE_OPENPMD_TEST_BACKENDS", raising=False)
    monkeypatch.setattr(openpmd_backend_matrix, "_front_end_backends", lambda: ["hdf5"])

    assert openpmd_backend_matrix.openpmd_test_backends() == ["hdf5"]


def test_backendMatrixAcceptsSelector(monkeypatch):
    monkeypatch.setenv("HASE_OPENPMD_TEST_BACKENDS", "hdf5, adios-sst")

    assert openpmd_backend_matrix.openpmd_test_backends() == ["adios-sst", "hdf5"]


def test_backendMatrixFailsWithoutProbe(monkeypatch):
    monkeypatch.delenv("HASE_OPENPMD_TEST_BACKENDS", raising=False)
    monkeypatch.setattr(openpmd_backend_matrix, "_front_end_backends", lambda: [])

    with pytest.raises(RuntimeError, match="frontend did not report any available openPMD backends"):
        openpmd_backend_matrix.openpmd_test_backends()


def test_backendMatrixPropagatesFrontEndProbeFailure(monkeypatch):
    monkeypatch.delenv("HASE_OPENPMD_TEST_BACKENDS", raising=False)

    def fail():
        raise RuntimeError("probe library missing")

    monkeypatch.setattr(openpmd_backend_matrix, "_front_end_backends", fail)

    with pytest.raises(RuntimeError, match="probe library missing"):
        openpmd_backend_matrix.openpmd_test_backends()


def test_backendMatrixDropsUnknownFrontEndBackends(monkeypatch):
    monkeypatch.delenv("HASE_OPENPMD_TEST_BACKENDS", raising=False)
    monkeypatch.setattr(openpmd_backend_matrix, "_front_end_backends", lambda: ["hdf5", "unsupported", "adios"])

    assert openpmd_backend_matrix.openpmd_test_backends() == ["adios", "hdf5"]
