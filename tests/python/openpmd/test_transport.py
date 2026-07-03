import os
import shlex
import subprocess
import sys
import textwrap
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import pyInclude.openpmd.transport as transport
from pyInclude import AlpakaBackends
from pyInclude.geometry import GainMedium, MeshTopology
from pyInclude.laser import CrossSectionData
from pyInclude.openpmd import HASE_TRANSPORT_VERSION, PrimitiveFieldSpec, PrismSchema, backendFlat, backendFlatArray, fieldSpec, haseTransportAttributes, primitiveView, spectralContext, unitDimension
from pyInclude.simulation import PhiASE




MESH_FIELD_VALUES = {
    "points": np.array([0.0, 1.5, 0.25, 2.25, -0.75, 0.0, -0.5, 1.25, 2.5, 3.75], dtype=np.float64),
    "connectivity": np.array([0, 2, 0, 1, 3, 2, 2, 4, 4], dtype=np.uint32),
    "neighbors": np.array([-1, 0, 1, 1, -1, 0, 2, 2, -1], dtype=np.int32),
    "forbiddenEdges": np.array([-1, -1, 2, 0, 2, -1, 1, -1, 0], dtype=np.int32),
    "normalPoints": np.array([4, 2, 0, 3, 1, 2, 0, 4, 3], dtype=np.uint32),
    "cellCenterX": np.array([0.5, 1.25, -0.25], dtype=np.float64),
    "cellCenterY": np.array([0.75, 1.5, 2.25], dtype=np.float64),
    "cellNormalX": np.array([0.10, 0.40, 0.70, 0.20, 0.50, 0.80, 0.30, 0.60, 0.90], dtype=np.float64),
    "cellNormalY": np.array([-0.10, -0.40, -0.70, -0.20, -0.50, -0.80, -0.30, -0.60, -0.90], dtype=np.float64),
    "surface": np.array([1.25, 2.50, 3.75], dtype=np.float32),
    "betaVolume": np.array([0.11, 0.21, 0.31, 0.12, 0.22, 0.32, 0.13, 0.23, 0.33, 0.14, 0.24, 0.34, 0.15, 0.25, 0.35], dtype=np.float64),
    "pointBeta": np.array([100.0 + 10.0 * point + level for level in range(6) for point in range(5)], dtype=np.float64),
    "claddingCellType": np.array([0, 2, 1], dtype=np.uint32),
    "refractiveIndex": np.array([1.80, 1.20, 1.65, 1.05], dtype=np.float32),
    "reflectivity": np.array([0.01, 0.03, 0.05, 0.02, 0.04, 0.06], dtype=np.float32),
}


def _launch_backend():
    backends = AlpakaBackends.all()
    for backend in backends:
        if "CpuOmpBlocks" in backend:
            return backend
    if not backends:
        pytest.skip("no Alpaka backend is available in this build")
    return backends[0]

SPECTRAL_FIELD_VALUES = {
    "lambdaAbsorption": np.array([900e-9, 910e-9, 930e-9], dtype=np.float64),
    "lambdaEmission": np.array([1000e-9, 1015e-9, 1040e-9], dtype=np.float64),
    "sigmaAbsorption": np.array([0.010, 0.025, 0.040], dtype=np.float64),
    "sigmaEmission": np.array([0.050, 0.035, 0.020], dtype=np.float64),
}

SCALAR_RECORD_SPECS = {
    "connectivity": "connectivity",
    "neighbors": "neighbors",
    "forbidden_edges": "forbiddenEdges",
    "normal_points": "normalPoints",
    "cell_normal_x": "cellNormalX",
    "cell_normal_y": "cellNormalY",
    "surface": "surface",
    "beta_volume": "betaVolume",
    "point_beta": "pointBeta",
    "cladding_cell_type": "claddingCellType",
    "refractive_index": "refractiveIndex",
    "reflectivity": "reflectivity",
    "lambda_absorption": "lambdaAbsorption",
    "lambda_emission": "lambdaEmission",
    "sigma_absorption": "sigmaAbsorption",
    "sigma_emission": "sigmaEmission",
}


def asymmetric_topology():
    points = np.column_stack((MESH_FIELD_VALUES["points"][:5], MESH_FIELD_VALUES["points"][5:]))
    triangles = MESH_FIELD_VALUES["connectivity"].reshape((3, 3), order="F")
    return MeshTopology(points, triangles, levels=6, thickness=0.375)


def asymmetric_medium():
    return GainMedium(asymmetric_topology()).withPhysicalProperties(
        betaVolume=backendFlat(MESH_FIELD_VALUES["betaVolume"]),
        betaCells=backendFlat(MESH_FIELD_VALUES["pointBeta"]),
        claddingCellTypes=MESH_FIELD_VALUES["claddingCellType"],
        refractiveIndices=MESH_FIELD_VALUES["refractiveIndex"],
        reflectivities=backendFlat(MESH_FIELD_VALUES["reflectivity"]),
        nTot=7.5,
        crystalTFluo=1.75,
        claddingNumber=3,
        claddingAbsorption=0.075,
    )


def asymmetric_cross_sections():
    return CrossSectionData(
        wavelengthsAbsorption=SPECTRAL_FIELD_VALUES["lambdaAbsorption"],
        crossSectionAbsorption=SPECTRAL_FIELD_VALUES["sigmaAbsorption"],
        wavelengthsEmission=SPECTRAL_FIELD_VALUES["lambdaEmission"],
        crossSectionEmission=SPECTRAL_FIELD_VALUES["sigmaEmission"],
        resolution=3,
    )


def asymmetric_phi_ase():
    return PhiASE(
        crossSections=asymmetric_cross_sections(),
        minRaysPerSample=1,
        maxRaysPerSample=1,
        mseThreshold=0.25,
        repetitions=1,
        adaptiveSteps=1,
        useReflections=True,
        backend=_launch_backend(),
        parallelMode="single",
        numDevices=1,
        minSampleRange=0,
        maxSampleRange=0,
        rngSeed=1234,
    )


def launch_smoke_topology():
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    triangles = np.array([[0, 1, 2]], dtype=np.uint32)
    return MeshTopology(points, triangles, levels=2, thickness=1.0)


def launch_smoke_medium():
    topology = launch_smoke_topology()
    return GainMedium(topology).withPhysicalProperties(
        betaVolume=backendFlat(np.array([0.0], dtype=np.float64)),
        betaCells=backendFlat(np.zeros(topology.numberOfPoints * topology.levels, dtype=np.float64)),
        claddingCellTypes=np.array([0], dtype=np.uint32),
        refractiveIndices=np.array([1.5, 1.0, 1.5, 1.0], dtype=np.float32),
        reflectivities=backendFlat(np.array([0.0, 0.0], dtype=np.float32)),
        nTot=1.0,
        crystalTFluo=1.0,
        claddingNumber=99,
        claddingAbsorption=0.0,
    )


def launch_smoke_cross_sections():
    return CrossSectionData(
        wavelengthsAbsorption=np.array([900e-9], dtype=np.float64),
        crossSectionAbsorption=np.array([0.0], dtype=np.float64),
        wavelengthsEmission=np.array([1000e-9], dtype=np.float64),
        crossSectionEmission=np.array([0.0], dtype=np.float64),
        resolution=1,
    )


def launch_smoke_phi_ase():
    cross_sections = launch_smoke_cross_sections()
    return PhiASE(
        crossSections=cross_sections,
        minRaysPerSample=1,
        maxRaysPerSample=1,
        mseThreshold=0.25,
        repetitions=1,
        adaptiveSteps=1,
        useReflections=False,
        backend=_launch_backend(),
        parallelMode="single",
        numDevices=1,
        minSampleRange=0,
        maxSampleRange=0,
        rngSeed=1234,
    )


def _file_backend_for_tests():
    backend = os.environ.get("HASE_OPENPMD_TEST_BACKEND", "adios").strip().lower()
    if backend in {"hdf5", "adios"}:
        return backend
    return "adios"


def _file_suffix_for_tests():
    return transport._backend_spec(_file_backend_for_tests()).suffix


def asymmetric_mesh():
    medium = asymmetric_medium()
    topology = medium.topology
    derived = topology._topology()
    return SimpleNamespace(
        numberOfPoints=topology.numberOfPoints,
        numberOfTriangles=topology.numberOfTriangles,
        numberOfLevels=int(topology.levels),
        points=np.asarray(topology.points).reshape(-1, order="F"),
        trianglePointIndices=np.asarray(topology.trianglePointIndices).reshape(-1, order="F"),
        triangleNeighbors=derived["triangleNeighbors"],
        forbiddenEdge=derived["forbiddenEdge"],
        triangleNormalPoint=derived["triangleNormalPoint"],
        triangleCenterX=derived["triangleCenterX"],
        triangleCenterY=derived["triangleCenterY"],
        triangleNormalsX=derived["triangleNormalsX"],
        triangleNormalsY=derived["triangleNormalsY"],
        triangleSurfaces=derived["triangleSurfaces"],
        betaVolume=medium.get("betaVolume").value,
        betaCells=medium.get("betaCells").value,
        claddingCellTypes=medium.get("claddingCellTypes").value,
        refractiveIndices=medium.get("refractiveIndices").value,
        reflectivities=medium.get("reflectivities").value,
    )


def _field_context():
    topology = asymmetric_topology()
    return SimpleNamespace(
        numberOfPoints=topology.numberOfPoints,
        numberOfTriangles=topology.numberOfTriangles,
        numberOfLevels=int(topology.levels),
    )


def _transport_scalar_record_values():
    return {
        "beta_volume": MESH_FIELD_VALUES["betaVolume"],
        "point_beta": MESH_FIELD_VALUES["pointBeta"],
        "cladding_cell_type": MESH_FIELD_VALUES["claddingCellType"],
        "refractive_index": MESH_FIELD_VALUES["refractiveIndex"],
        "reflectivity": MESH_FIELD_VALUES["reflectivity"],
        "lambda_absorption": SPECTRAL_FIELD_VALUES["lambdaAbsorption"],
        "lambda_emission": SPECTRAL_FIELD_VALUES["lambdaEmission"],
        "sigma_absorption": SPECTRAL_FIELD_VALUES["sigmaAbsorption"],
        "sigma_emission": SPECTRAL_FIELD_VALUES["sigmaEmission"],
    }


def _mesh_field_values(mesh):
    return {
        "points": np.asarray(mesh.points),
        "connectivity": np.asarray(mesh.trianglePointIndices),
        "neighbors": np.asarray(mesh.triangleNeighbors),
        "forbiddenEdges": np.asarray(mesh.forbiddenEdge),
        "normalPoints": np.asarray(mesh.triangleNormalPoint),
        "cellCenterX": np.asarray(mesh.triangleCenterX),
        "cellCenterY": np.asarray(mesh.triangleCenterY),
        "cellNormalX": np.asarray(mesh.triangleNormalsX),
        "cellNormalY": np.asarray(mesh.triangleNormalsY),
        "surface": np.asarray(mesh.triangleSurfaces),
        "betaVolume": np.asarray(mesh.betaVolume),
        "pointBeta": np.asarray(mesh.betaCells),
        "claddingCellType": np.asarray(mesh.claddingCellTypes),
        "refractiveIndex": np.asarray(mesh.refractiveIndices),
        "reflectivity": np.asarray(mesh.reflectivities),
    }


def _scalar_record_values(mesh):
    values = _mesh_field_values(mesh)
    return {
        "connectivity": values["connectivity"],
        "neighbors": values["neighbors"],
        "forbidden_edges": values["forbiddenEdges"],
        "normal_points": values["normalPoints"],
        "cell_normal_x": values["cellNormalX"],
        "cell_normal_y": values["cellNormalY"],
        "surface": values["surface"],
        "beta_volume": values["betaVolume"],
        "point_beta": values["pointBeta"],
        "cladding_cell_type": values["claddingCellType"],
        "refractive_index": values["refractiveIndex"],
        "reflectivity": values["reflectivity"],
        "lambda_absorption": SPECTRAL_FIELD_VALUES["lambdaAbsorption"],
        "lambda_emission": SPECTRAL_FIELD_VALUES["lambdaEmission"],
        "sigma_absorption": SPECTRAL_FIELD_VALUES["sigmaAbsorption"],
        "sigma_emission": SPECTRAL_FIELD_VALUES["sigmaEmission"],
    }


@pytest.fixture
def contract_input(tmp_path):
    _require_openpmd_transport_io()
    output = tmp_path / ("contract" + _file_suffix_for_tests())
    with transport.OpenPmdInputSeries(output, backend=_file_backend_for_tests()) as writer:
        writer.write(asymmetric_phi_ase(), asymmetric_medium(), asymmetric_cross_sections())
    return output


def _io():
    try:
        return transport._io()
    except RuntimeError as exc:
        if "CMake-built openpmd_api" in str(exc):
            pytest.skip(str(exc))
        raise


def _require_openpmd_transport_io():
    _io()


def _openpmd_subprocess_env():
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    openpmd_path = str(Path(_io().__file__).resolve().parents[1])
    pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = openpmd_path if not pythonpath else os.pathsep.join((openpmd_path, pythonpath))
    return env


def _read_scalar(series, iteration, name):
    io = _io()
    component = iteration.meshes[name][io.Mesh_Record_Component.SCALAR]
    chunk = component.load_chunk()
    series.flush()
    return np.array(chunk, copy=True)


def _read_component(series, component):
    chunk = component.load_chunk()
    series.flush()
    return np.array(chunk, copy=True)


def _attribute_list(value):
    if isinstance(value, str):
        return [value]
    try:
        return list(value)
    except TypeError:
        return [value]


def _record_metadata(record, spec):
    return {
        "transport_version": record.get_attribute("haseTransportVersion"),
        "entity": record.get_attribute("haseEntity"),
        "axes": _attribute_list(record.get_attribute("haseAxes")),
        "axes_string": record.get_attribute("haseAxesString"),
        "layout": record.get_attribute("haseLayoutOrder"),
        "primitive_shape": _attribute_list(record.get_attribute("hasePrimitiveShape")),
        "static": record.get_attribute("haseStatic"),
        "dynamic": record.get_attribute("haseDynamic"),
        "backend_required": record.get_attribute("haseBackendRequired"),
        "unit": record.get_attribute("haseUnit"),
        "axis_labels": list(record.axis_labels),
        "axis_labels_string": record.get_attribute("haseAxisLabelsString"),
        "unit_si": record[_io().Mesh_Record_Component.SCALAR].unit_SI,
    }


def _unit_dimension_values(record, io):
    value = record.unit_dimension
    if isinstance(value, (list, tuple)):
        return tuple(float(item) for item in value)
    labels = (
        io.Unit_Dimension.L,
        io.Unit_Dimension.M,
        io.Unit_Dimension.T,
        io.Unit_Dimension.I,
        io.Unit_Dimension.theta,
        io.Unit_Dimension.N,
        io.Unit_Dimension.J,
    )
    return tuple(float(value.get(label, 0.0)) for label in labels)


def _assert_base_openpmd_scalar_metadata(record, *, axis_labels, unit_si=1.0, unit_dimension=None):
    io = _io()
    expected_dimension = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) if unit_dimension is None else unit_dimension

    assert record.get_attribute("geometry") == "other"
    assert record.get_attribute("dataOrder") == "C"
    assert list(record.axis_labels) == list(axis_labels)
    assert list(record.grid_spacing) == [1.0] * len(axis_labels)
    assert list(record.grid_global_offset) == [0.0] * len(axis_labels)
    assert record.grid_unit_SI == 1.0
    assert _unit_dimension_values(record, io) == expected_dimension

    component = record[io.Mesh_Record_Component.SCALAR]
    assert component.unit_SI == unit_si
    assert list(component.position) == [0.0] * len(axis_labels)


def _assert_hase_metadata(record, spec, context):
    assert _record_metadata(record, spec) == {
        "transport_version": HASE_TRANSPORT_VERSION,
        "entity": spec.entity,
        "axes": list(spec.axes),
        "axes_string": ",".join(spec.axes),
        "layout": "backendFlat",
        "primitive_shape": list(spec.expectedShape(context)),
        "static": not spec.dynamic,
        "dynamic": spec.dynamic,
        "backend_required": spec.backendRequired,
        "unit": spec.unit,
        "axis_labels": ["flatIndex"],
        "axis_labels_string": "flatIndex",
        "unit_si": spec.unitSI,
    }
    assert _unit_dimension_values(record, _io()) == spec.unitDimension


def _context_for_spec(spec_name):
    if spec_name in SPECTRAL_FIELD_VALUES:
        return spectralContext(SPECTRAL_FIELD_VALUES[spec_name])
    return _field_context()


def test_layout_helpers_define_exact_backend_flat_contract():
    mesh = asymmetric_mesh()
    for name, values in _mesh_field_values(mesh).items():
        spec = fieldSpec(name)
        flat = backendFlatArray(values, spec, mesh, layoutOrder="backendFlat")
        np.testing.assert_array_equal(flat, values.astype(spec.dtypeObject, copy=False))
        view = primitiveView(backendFlat(values), spec, mesh)
        assert view.shape == spec.expectedShape(mesh)
        assert view.dtype == spec.dtypeObject
        np.testing.assert_array_equal(view.reshape(-1, order="F"), flat)

    for name, values in SPECTRAL_FIELD_VALUES.items():
        spec = fieldSpec(name)
        context = spectralContext(values)
        flat = backendFlatArray(values, spec, context, layoutOrder="backendFlat")
        np.testing.assert_array_equal(flat, values.astype(spec.dtypeObject, copy=False))
        assert primitiveView(backendFlat(values), spec, context).shape == spec.expectedShape(context)


def test_openpmd_backend_names_map_to_expected_suffixes_and_configs(monkeypatch):
    assert transport._backend_spec().name == "adios-sst"
    assert transport._backend_spec("adios").suffix == ".bp"
    assert transport._backend_spec("adios").config == {"backend": "adios2"}
    assert transport._backend_spec("adios-sst").suffix == ".sst"
    sst_engine = transport._backend_spec("adios-sst").config["adios2"]["engine"]
    assert sst_engine["type"] == "sst"
    assert "QueueFullPolicy" not in sst_engine["parameters"]
    assert transport._backend_spec("adios-sst").streaming is True
    assert transport._backend_spec("hdf5").config == {"backend": "hdf5"}

    monkeypatch.setenv("HASE_OPENPMD_TEST_BACKEND", "hdf5")
    assert _file_backend_for_tests() == "hdf5"
    monkeypatch.setenv("HASE_OPENPMD_TEST_BACKEND", "adios")
    assert _file_backend_for_tests() == "adios"
    monkeypatch.setenv("HASE_OPENPMD_TEST_BACKEND", "adios-sst")
    assert _file_backend_for_tests() == "adios"

    for backend in ("bp", "unsupported"):
        with pytest.raises(ValueError, match="unsupported openPMD backend"):
            transport._backend_spec(backend)


def test_openpmd_watchdog_interval_defaults_to_thirty_seconds(monkeypatch):
    monkeypatch.delenv("HASE_OPENPMD_WATCHDOG_INTERVAL", raising=False)

    assert transport._watchdog_interval() == 30.0


def test_openpmd_watchdog_interval_accepts_explicit_and_environment_values(monkeypatch):
    assert transport._watchdog_interval(12) == 12.0
    assert transport._watchdog_interval("2.5") == 2.5
    assert transport._watchdog_interval("none") is None

    monkeypatch.setenv("HASE_OPENPMD_WATCHDOG_INTERVAL", "45")
    assert transport._watchdog_interval() == 45.0
    monkeypatch.setenv("HASE_OPENPMD_WATCHDOG_INTERVAL", "0")
    assert transport._watchdog_interval() is None


def test_openpmd_watchdog_interval_rejects_invalid_values(monkeypatch):
    monkeypatch.setenv("HASE_OPENPMD_WATCHDOG_INTERVAL", "not-a-number")
    with pytest.raises(ValueError, match="HASE_OPENPMD_WATCHDOG_INTERVAL"):
        transport._watchdog_interval()

    with pytest.raises(ValueError, match="HASE_OPENPMD_WATCHDOG_INTERVAL"):
        transport._watchdog_interval(-1)


def test_streaming_wait_ignores_watchdog_alive_events():
    session = transport.OpenPmdPhiAseSession(transport="adios-sst", watchdog_interval="none")
    result = object()
    session._sender_errors = transport.queue.Queue()
    session._watchdog_events = transport.queue.Queue()
    session._result_queue = transport.queue.Queue()
    session._proc = SimpleNamespace(poll=lambda: None)

    session._watchdog_events.put((True, None))
    session._result_queue.put((True, (0, result)))

    assert session._wait_for_result(0) is result


def test_streaming_wait_raises_on_watchdog_failure():
    session = transport.OpenPmdPhiAseSession(transport="adios-sst", watchdog_interval="none")
    session._sender_errors = transport.queue.Queue()
    session._watchdog_events = transport.queue.Queue()
    session._result_queue = transport.queue.Queue()
    session._proc = SimpleNamespace(poll=lambda: None, stderr=None)

    session._watchdog_events.put((False, RuntimeError("backend liveness probe failed")))

    with pytest.raises(RuntimeError, match="openPMD backend watchdog failed"):
        session._wait_for_result(0)


def test_streaming_reader_queues_eof_sentinel_when_stream_ends(monkeypatch, tmp_path):
    session = transport.OpenPmdPhiAseSession(transport="adios-sst", watchdog_interval="none")
    session._result_queue = transport.queue.Queue()
    session._output_path = tmp_path / "results.sst"

    result = object()
    closed = []

    class FakeSeries:
        def __init__(self, path, access, config):
            self.path = path
            self.access = access
            self.config = config

        def read_iterations(self):
            return iter([object()])

        def close(self):
            closed.append(True)

    fake_io = SimpleNamespace(Series=FakeSeries)
    monkeypatch.setattr(transport, "_io", lambda: fake_io)
    monkeypatch.setattr(transport, "_access", lambda name: name)
    monkeypatch.setattr(transport, "_series_config", lambda path: {"config": str(path)})
    monkeypatch.setattr(
        transport,
        "_read_result_iteration",
        lambda series, iteration, fallback_index=None: (fallback_index, result),
    )

    session._read_streaming_results()

    assert session._result_queue.get_nowait() == (True, (0, result))
    ok, payload = session._result_queue.get_nowait()
    assert ok is True
    assert payload is transport._STREAMING_RESULT_EOF
    assert closed == [True]


def test_streaming_wait_fails_when_result_stream_ends_without_expected_iteration():
    session = transport.OpenPmdPhiAseSession(transport="adios-sst", watchdog_interval="none")
    session._sender_errors = transport.queue.Queue()
    session._watchdog_events = transport.queue.Queue()
    session._result_queue = transport.queue.Queue()
    session._proc = SimpleNamespace(poll=lambda: None)

    session._result_queue.put((True, (1, object())))
    session._result_queue.put((True, transport._STREAMING_RESULT_EOF))

    message = r"Expected result iteration 0.*result stream ended.*buffered iterations: \[1\]"
    with pytest.raises(RuntimeError, match=message):
        session._wait_for_result(0)


def test_streaming_wait_fails_when_backend_finishes_before_expected(monkeypatch):
    monkeypatch.setattr(transport, "_STREAMING_RESULT_POLL_SECONDS", 0.001)
    session = transport.OpenPmdPhiAseSession(transport="adios-sst", watchdog_interval="none")
    session._sender_errors = transport.queue.Queue()
    session._watchdog_events = transport.queue.Queue()
    session._result_queue = transport.queue.Queue()
    session._proc = SimpleNamespace(poll=lambda: 0)

    with pytest.raises(RuntimeError, match="calcPhiASE completed before result iteration 2"):
        session._wait_for_result(2)


def test_streaming_wait_fails_when_reader_thread_stops(monkeypatch):
    monkeypatch.setattr(transport, "_STREAMING_RESULT_POLL_SECONDS", 0.001)
    session = transport.OpenPmdPhiAseSession(transport="adios-sst", watchdog_interval="none")
    session._sender_errors = transport.queue.Queue()
    session._watchdog_events = transport.queue.Queue()
    session._result_queue = transport.queue.Queue()
    session._proc = SimpleNamespace(poll=lambda: None)
    session._reader = SimpleNamespace(is_alive=lambda: False)

    with pytest.raises(RuntimeError, match="result receiver thread stopped"):
        session._wait_for_result(3)


def test_streaming_thread_join_uses_finite_timeout(monkeypatch):
    monkeypatch.setenv("HASE_OPENPMD_THREAD_JOIN_TIMEOUT", "0.25")
    session = transport.OpenPmdPhiAseSession(transport="adios-sst", watchdog_interval="none")

    class HangingThread:
        def __init__(self):
            self.timeout = None

        def join(self, timeout=None):
            self.timeout = timeout

        def is_alive(self):
            return True

    thread = HangingThread()
    error = session._join_streaming_thread(thread, "result receiver")

    assert thread.timeout == 0.25
    assert isinstance(error, RuntimeError)
    assert "did not stop within 0.25 seconds" in str(error)


@pytest.mark.integration
def test_adios_sst_python_pair_publishes_one_iteration(tmp_path):
    io = _io()
    if "sst" not in getattr(io, "file_extensions", []):
        pytest.skip("openPMD-api was not built with ADIOS2 SST support")

    stream = tmp_path / "sst_pair.sst"
    config = repr(transport.SST_CONFIG)
    values = [1.25, 2.5, 3.75]
    reader_code = textwrap.dedent(
        f"""
        import numpy as np
        import openpmd_api as io
        import sys

        def access(name):
            return getattr(io.Access_Type, name) if hasattr(io, "Access_Type") else getattr(io.Access, name)

        series = io.Series(sys.argv[1], access("read_linear"), {config})
        try:
            iterations = series.read_iterations() if hasattr(series, "read_iterations") else series.snapshots().items()
            item = next(iter(iterations))
            iteration = item[1] if isinstance(item, tuple) else item
            component = iteration.meshes["probe"][io.Mesh_Record_Component.SCALAR]
            chunk = component.load_chunk()
            series.flush()
            data = np.asarray(chunk, dtype=np.float64).reshape(-1)
            iteration.close()
        finally:
            series.close()

        expected = np.asarray({values!r}, dtype=np.float64)
        if not np.allclose(data, expected):
            raise SystemExit(f"unexpected SST payload: {{data!r}}")
        print("reader done", flush=True)
        """
    )
    writer_code = textwrap.dedent(
        f"""
        import numpy as np
        import openpmd_api as io
        import sys

        def access(name):
            return getattr(io.Access_Type, name) if hasattr(io, "Access_Type") else getattr(io.Access, name)

        series = io.Series(sys.argv[1], access("create_linear"), {config})
        try:
            iteration = series.snapshots()[0]
            data = np.asarray({values!r}, dtype=np.float64)
            component = iteration.meshes["probe"][io.Mesh_Record_Component.SCALAR]
            component.reset_dataset(io.Dataset(data.dtype, data.shape))
            component.store_chunk(data)
            iteration.close()
        finally:
            series.close()
        print("writer done", flush=True)
        """
    )

    def finish(process, name, timeout=20):
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            pytest.fail(f"{name} timed out\nstdout:\n{stdout}\nstderr:\n{stderr}")
        assert process.returncode == 0, f"{name} failed\nstdout:\n{stdout}\nstderr:\n{stderr}"
        return stdout, stderr

    env = _openpmd_subprocess_env()
    reader = subprocess.Popen(
        [sys.executable, "-u", "-c", reader_code, str(stream)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    time.sleep(0.5)
    writer = subprocess.Popen(
        [sys.executable, "-u", "-c", writer_code, str(stream)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )

    try:
        writer_stdout, _ = finish(writer, "SST writer")
        reader_stdout, _ = finish(reader, "SST reader")
    finally:
        for process in (writer, reader):
            if process.poll() is None:
                process.kill()

    assert "writer done" in writer_stdout
    assert "reader done" in reader_stdout



@pytest.mark.integration
def test_adios_sst_threads_can_send_and_receive_independent_streams(tmp_path):
    io = _io()
    if "sst" not in getattr(io, "file_extensions", []):
        pytest.skip("openPMD-api was not built with ADIOS2 SST support")

    input_stream = tmp_path / "threaded_input.sst"
    output_stream = tmp_path / "threaded_output.sst"
    config = repr(transport.SST_CONFIG)
    input_values = [1.0, 2.0, 3.0]
    output_values = [10.0, 20.0, 30.0]

    common_code = f"""
import numpy as np
import openpmd_api as io

config = {config}

def access(name):
    return getattr(io.Access_Type, name) if hasattr(io, "Access_Type") else getattr(io.Access, name)

def write_probe(series, values):
    iteration = series.snapshots()[0]
    data = np.asarray(values, dtype=np.float64)
    component = iteration.meshes["probe"][io.Mesh_Record_Component.SCALAR]
    component.reset_dataset(io.Dataset(data.dtype, data.shape))
    component.store_chunk(data)
    iteration.close()
    series.flush()

def read_probe(series):
    iterations = series.read_iterations() if hasattr(series, "read_iterations") else series.snapshots().items()
    item = next(iter(iterations))
    iteration = item[1] if isinstance(item, tuple) else item
    component = iteration.meshes["probe"][io.Mesh_Record_Component.SCALAR]
    chunk = component.load_chunk()
    series.flush()
    data = np.asarray(chunk, dtype=np.float64).reshape(-1).copy()
    iteration.close()
    return data
"""
    sender_code = textwrap.dedent(
        common_code
        + f"""
import sys

series = io.Series(sys.argv[1], access("create_linear"), config)
try:
    write_probe(series, {input_values!r})
finally:
    series.close()
print("sender done", flush=True)
"""
    )
    backend_code = textwrap.dedent(
        common_code
        + f"""
import sys

input_series = io.Series(sys.argv[1], access("read_linear"), config)
try:
    data = read_probe(input_series)
finally:
    input_series.close()
np.testing.assert_allclose(data, np.asarray({input_values!r}, dtype=np.float64))
print("backend input read", flush=True)

output_series = io.Series(sys.argv[2], access("create_linear"), config)
try:
    write_probe(output_series, {output_values!r})
    print("backend output written", flush=True)
    sys.stdin.read(1)
finally:
    output_series.close()
"""
    )
    receiver_code = textwrap.dedent(
        common_code
        + f"""
import sys

series = io.Series(sys.argv[1], access("read_linear"), config)
try:
    data = read_probe(series)
finally:
    series.close()
np.testing.assert_allclose(data, np.asarray({output_values!r}, dtype=np.float64))
print("receiver done", flush=True)
"""
    )

    def finish(process, name, timeout=20):
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            pytest.fail(f"{name} timed out\nstdout:\n{stdout}\nstderr:\n{stderr}")
        assert process.returncode == 0, f"{name} failed\nstdout:\n{stdout}\nstderr:\n{stderr}"
        return stdout, stderr

    env = _openpmd_subprocess_env()
    backend = subprocess.Popen(
        [sys.executable, "-u", "-c", backend_code, str(input_stream), str(output_stream)],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    time.sleep(0.5)
    sender = subprocess.Popen(
        [sys.executable, "-u", "-c", sender_code, str(input_stream)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )

    try:
        sender_stdout, _ = finish(sender, "SST sender")
        receiver = subprocess.Popen(
            [sys.executable, "-u", "-c", receiver_code, str(output_stream)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
        receiver_stdout, _ = finish(receiver, "SST receiver")
        try:
            backend_stdout, backend_stderr = backend.communicate(input="\n", timeout=20)
        except subprocess.TimeoutExpired:
            backend.kill()
            backend_stdout, backend_stderr = backend.communicate()
            pytest.fail(f"SST backend timed out\nstdout:\n{backend_stdout}\nstderr:\n{backend_stderr}")
        assert backend.returncode == 0, f"SST backend failed\nstdout:\n{backend_stdout}\nstderr:\n{backend_stderr}"
    finally:
        for process in (sender, backend):
            if process.poll() is None:
                process.kill()

    assert "sender done" in sender_stdout
    assert "backend input read" in backend_stdout
    assert "backend output written" in backend_stdout
    assert "receiver done" in receiver_stdout


def test_run_phi_ase_uses_openpmd_session_manager(monkeypatch):
    result = object()
    events = []

    class FakeSession:
        def __init__(self, *, transport=None):
            events.append(("init", transport))

        def __enter__(self):
            events.append(("enter",))
            return self

        def __exit__(self, exc_type, exc, traceback):
            events.append(("exit", exc_type))
            return False

        def run(self, phi_ase, gain_medium, cross_sections):
            events.append(("run", phi_ase, gain_medium, cross_sections))
            return result

    monkeypatch.setattr(transport, "OpenPmdPhiAseSession", FakeSession)
    phi_ase = object()
    medium = object()
    cross_sections = object()

    assert transport.runPhiASE(phi_ase, medium, cross_sections, transport="adios") is result
    assert events == [
        ("init", "adios"),
        ("enter",),
        ("run", phi_ase, medium, cross_sections),
        ("exit", None),
    ]


def testRunPhiAseUsesProvidedOpenPmdSessionWithoutClosingIt(monkeypatch):
    result = object()
    events = []

    class FakeSession:
        def run(self, phiAse, gainMedium, crossSections):
            events.append(("run", phiAse, gainMedium, crossSections))
            return result

    monkeypatch.setattr(
        transport,
        "OpenPmdPhiAseSession",
        lambda **kwargs: pytest.fail("provided openpmdSession should bypass session construction"),
    )

    phiAse = object()
    gainMedium = object()
    crossSections = object()
    openpmdSession = FakeSession()

    assert transport.runPhiASE(
        phiAse,
        gainMedium,
        crossSections,
        openpmdSession=openpmdSession,
    ) is result
    assert events == [("run", phiAse, gainMedium, crossSections)]


def test_openpmd_session_assigns_monotonic_request_iterations(monkeypatch):
    calls = []

    monkeypatch.setattr(transport, "findCalcPhiAse", lambda: Path("calcPhiASE"))
    monkeypatch.setattr(transport, "_ensure_backend_available", lambda backend: None)

    def fake_run_file_iteration(self, iteration_index, phi_ase, gain_medium, cross_sections):
        calls.append((iteration_index, phi_ase, gain_medium, cross_sections))
        return SimpleNamespace(iteration=iteration_index)

    monkeypatch.setattr(transport.OpenPmdPhiAseSession, "_run_file_iteration", fake_run_file_iteration)

    with transport.OpenPmdPhiAseSession(transport="adios") as session:
        first = session.run("phi0", "medium0", "cross0")
        second = session.run("phi1", "medium1", "cross1")

    assert first.iteration == 0
    assert second.iteration == 1
    assert calls == [(0, "phi0", "medium0", "cross0"), (1, "phi1", "medium1", "cross1")]


def test_forward_backend_logging_replays_streams_when_enabled(monkeypatch, capsys):
    monkeypatch.setenv("HASE_FORWARD_LOGGING", "ON")

    transport._forward_backend_logging(stdout="backend stdout\n", stderr="backend stderr\n")

    captured = capsys.readouterr()
    assert captured.out == "backend stdout\n"
    assert captured.err == "backend stderr\n"


def test_forward_backend_logging_is_quiet_when_disabled(monkeypatch, capsys):
    monkeypatch.setenv("HASE_FORWARD_LOGGING", "OFF")

    transport._forward_backend_logging(stdout="backend stdout\n", stderr="backend stderr\n")

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_backend_failure_detail_includes_stdout_and_stderr():
    detail = transport._backend_failure_detail(
        stdout=" parser progress \n",
        stderr=" backend error \n",
    )

    assert "calcPhiASE stdout:\nparser progress" in detail
    assert "calcPhiASE stderr:\nbackend error" in detail


def test_openpmd_api_preference_rejects_missing_cmake_selected_module(monkeypatch, tmp_path):
    active = tmp_path / "site-packages" / "openpmd_api" / "__init__.py"
    monkeypatch.setitem(transport.sys.modules, "openpmd_api", SimpleNamespace(__file__=str(active)))
    monkeypatch.setattr(transport, "_candidate_python_paths", lambda executable: iter([tmp_path / "missing"]))
    monkeypatch.setattr(
        transport,
        "_binding_config",
        lambda: SimpleNamespace(HASE_USE_SYSTEM_OPENPMD=False),
    )

    with pytest.raises(RuntimeError, match=r"compatible with the openPMD C\+\+ provider"):
        transport._prefer_matching_openpmd_api(Path("calcPhiASE"))


def test_openpmd_api_preference_allows_normal_import_for_external_openpmd(monkeypatch):
    monkeypatch.setattr(transport, "_candidate_python_paths", lambda executable: iter([]))
    monkeypatch.setattr(
        transport,
        "_binding_config",
        lambda: SimpleNamespace(HASE_USE_SYSTEM_OPENPMD=True),
    )

    transport._prefer_matching_openpmd_api(Path("calcPhiASE"))


def test_openpmd_api_candidate_paths_include_configured_provider(monkeypatch, tmp_path):
    configured = tmp_path / "site-packages"
    (configured / "openpmd_api").mkdir(parents=True)
    monkeypatch.delenv("HASE_OPENPMD_PYTHONPATH", raising=False)
    monkeypatch.delenv("HASE_OPENPMD_PYTHON_PACKAGE_DIR", raising=False)
    monkeypatch.setattr(
        transport,
        "_binding_config",
        lambda: SimpleNamespace(HASE_OPENPMD_PYTHON_PACKAGE_DIR=str(configured)),
    )

    candidates = list(transport._candidate_python_paths(tmp_path / "build" / "calcPhiASE"))

    assert configured in candidates


def test_openpmd_api_candidate_paths_include_explicit_runtime_override(monkeypatch, tmp_path):
    configured = tmp_path / "configured"
    runtime_package = tmp_path / "runtime" / "openpmd_api"
    runtime_package.mkdir(parents=True)
    monkeypatch.setenv("HASE_OPENPMD_PYTHONPATH", str(runtime_package))
    monkeypatch.setattr(
        transport,
        "_binding_config",
        lambda: SimpleNamespace(HASE_OPENPMD_PYTHON_PACKAGE_DIR=str(configured)),
    )

    candidates = list(transport._candidate_python_paths(tmp_path / "build" / "calcPhiASE"))

    assert candidates[0] == runtime_package.parent


def test_openpmd_api_preference_rejects_mismatched_bundled_build(monkeypatch, tmp_path):
    candidate = tmp_path / "build" / "site-packages"
    candidate.mkdir(parents=True)
    active = tmp_path / "other" / "site-packages" / "openpmd_api" / "__init__.py"
    monkeypatch.setitem(transport.sys.modules, "openpmd_api", SimpleNamespace(__file__=str(active)))
    monkeypatch.setattr(transport, "_candidate_python_paths", lambda executable: iter([candidate]))

    with pytest.raises(RuntimeError, match="same openPMD-api build"):
        transport._prefer_matching_openpmd_api(Path("calcPhiASE"))


def test_layout_helpers_reject_accidental_transpose_views():
    mesh = asymmetric_mesh()
    spec = fieldSpec("betaVolume")
    primitive = primitiveView(backendFlat(MESH_FIELD_VALUES["betaVolume"]), spec, mesh)
    transposed_same_size = np.asfortranarray(primitive.T)
    with pytest.raises(ValueError, match="expects primitive shape"):
        backendFlatArray(transposed_same_size, spec, mesh)


def test_openpmd_input_series_contains_exact_values_order_shape_dtype_units_and_metadata(contract_input):
    io = _io()
    series = io.Series(str(contract_input), io.Access.read_only)
    for name, value in haseTransportAttributes.items():
        assert series.get_attribute(name) == value
    assert series.get_attribute("haseTransportVersion") == HASE_TRANSPORT_VERSION
    iteration = series.iterations[0]

    assert iteration.get_attribute("number_of_points") == 5
    assert iteration.get_attribute("number_of_cells") == 3
    assert iteration.get_attribute("number_of_levels") == 6
    assert iteration.get_attribute("thickness") == pytest.approx(0.375)
    assert iteration.get_attribute("n_tot") == pytest.approx(7.5)
    assert iteration.get_attribute("crystal_t_fluo") == pytest.approx(1.75)
    assert iteration.get_attribute("spectral_resolution") == 3
    assert iteration.get_attribute("rng_seed") == 1234

    canonical_context = transport._canonical_topology_context(asymmetric_topology())
    points = iteration.meshes["core_points"]
    assert points.get_attribute("haseTransportVersion") == HASE_TRANSPORT_VERSION
    assert points.get_attribute("haseEntity") == "coordinate_mesh_point"
    assert _attribute_list(points.get_attribute("haseAxes")) == ["coordinate", "mesh_point"]
    assert _attribute_list(points.get_attribute("hasePrimitiveShape")) == [3, 30]
    assert points.get_attribute("haseUnit") == "m"
    assert points.unit_dimension[io.Unit_Dimension.L] == 1.0
    components = transport._canonical_point_components(asymmetric_topology())
    np.testing.assert_array_equal(_read_component(series, points["x"]), components["x"])
    np.testing.assert_array_equal(_read_component(series, points["y"]), components["y"])
    np.testing.assert_array_equal(_read_component(series, points["z"]), components["z"])
    series.flush()

    canonical_scalars = {
        "cells_connectivity": (transport.CANONICAL_CONNECTIVITY_SPEC, transport._canonical_cell_connectivity(asymmetric_topology())),
        "cells_offsets": (
            transport.CANONICAL_OFFSETS_SPEC,
            np.arange(canonical_context.numberOfPrisms + 1, dtype=np.uint32) * np.uint32(6),
        ),
        "cells_types": (
            transport.CANONICAL_CELL_TYPES_SPEC,
            np.full(canonical_context.numberOfPrisms, 13, dtype=np.uint32),
        ),
    }
    for record_name, (spec, expected) in canonical_scalars.items():
        record = iteration.meshes["core_" + record_name]
        values = _read_scalar(series, iteration, "core_" + record_name)
        np.testing.assert_array_equal(values, expected.astype(spec.dtypeObject, copy=False))
        _assert_hase_metadata(record, spec, canonical_context)

    for removed_record in [
        "core_vertices",
        "core_connectivity",
        "core_neighbors",
        "core_forbidden_edges",
        "core_normal_points",
        "core_cell_center",
        "core_cell_normal_x",
        "core_cell_normal_y",
        "core_surface",
    ]:
        assert removed_record not in iteration.meshes

    present_specs = {
        name: spec_name
        for name, spec_name in SCALAR_RECORD_SPECS.items()
        if name not in {
            "connectivity",
            "neighbors",
            "forbidden_edges",
            "normal_points",
            "cell_normal_x",
            "cell_normal_y",
            "surface",
        }
    }
    for record_name, spec_name in present_specs.items():
        spec = fieldSpec(spec_name)
        context = _context_for_spec(spec_name)
        record = iteration.meshes["core_" + record_name]
        values = _read_scalar(series, iteration, "core_" + record_name)
        expected = _transport_scalar_record_values()[record_name].astype(spec.dtypeObject, copy=False)
        assert values.shape == (expected.size,)
        assert values.dtype == spec.dtypeObject
        np.testing.assert_array_equal(values, expected)
        _assert_hase_metadata(record, spec, context)

    series.close()


def _build_dir_candidates():
    root = Path(__file__).resolve().parents[3]
    env = os.environ.get("BUILD_DIR")
    if env:
        path = Path(env)
        yield path if path.is_absolute() else root / path

    yield from sorted(root.glob("build/cp*"))
    yield root / "build"
    yield root / "build" / "ci"


def _parser_validation_candidates():
    env = os.environ.get("HASE_OPENPMD_PARSER_VALIDATION")
    if env:
        yield Path(env)
    for build_dir in _build_dir_candidates():
        yield build_dir / "tests" / "tests_openpmdParserValidation"


def _parser_validation_binary():
    configured = os.environ.get("HASE_OPENPMD_PARSER_VALIDATION")
    for helper in _parser_validation_candidates():
        if helper.is_file() and os.access(helper, os.X_OK):
            return helper
    if configured:
        pytest.fail(f"configured HASE_OPENPMD_PARSER_VALIDATION is not executable: {configured}")
    pytest.skip("no openPMD parser validation binary found; build tests_openpmdParserValidation or set HASE_OPENPMD_PARSER_VALIDATION")


def test_openpmd_input_series_omits_static_topology_after_first_iteration(tmp_path):
    _require_openpmd_transport_io()
    output = tmp_path / ("dynamic_split" + _file_suffix_for_tests())
    with transport.OpenPmdInputSeries(output, backend=_file_backend_for_tests()) as writer:
        writer.write(asymmetric_phi_ase(), asymmetric_medium(), asymmetric_cross_sections())
        writer.write(asymmetric_phi_ase(), asymmetric_medium(), asymmetric_cross_sections())

    series = _io().Series(str(output), _io().Access.read_only)
    first = series.iterations[0]
    second = series.iterations[1]
    assert first.get_attribute("haseStaticUpdate") is True
    assert second.get_attribute("haseStaticUpdate") is False
    assert "core_points" in first.meshes
    assert "core_cells_connectivity" in first.meshes
    assert "core_points" not in second.meshes
    assert "core_cells_connectivity" not in second.meshes
    assert sorted(second.meshes) == ["core_beta_volume", "core_point_beta"]
    series.close()


def test_python_writer_openpmd_cpp_parser_result_round_trip(contract_input, tmp_path):
    output = tmp_path / ("round_trip_result" + _file_suffix_for_tests())
    env = os.environ.copy()
    env["HASE_OPENPMD_PYTHON_CONTRACT_INPUT"] = str(contract_input)
    env["HASE_OPENPMD_PYTHON_CONTRACT_OUTPUT"] = str(output)
    helper = _parser_validation_binary()
    completed = subprocess.run(
        [str(helper), "openPMD parser round-trips a Python writer contract input"],
        check=False,
        text=True,
        capture_output=True,
        env=env,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr

    result = transport.read_result(output)
    expected_phi = np.array([0.5 + i for i in range(30)], dtype=np.float32)
    expected_mse = np.array([1000.0 + i for i in range(30)], dtype=np.float64)
    expected_total_rays = np.array([200 + i for i in range(30)], dtype=np.uint32)
    expected_dndt_ase = np.array([-10.0 - i for i in range(30)], dtype=np.float64)
    np.testing.assert_array_equal(result.phiAse, expected_phi)
    np.testing.assert_array_equal(result.mse, expected_mse)
    np.testing.assert_array_equal(result.totalRays, expected_total_rays)
    np.testing.assert_array_equal(result.dndtAse, expected_dndt_ase)

    io = _io()
    series = io.Series(str(output), io.Access.read_only)
    iteration = series.iterations[0]
    expected_units = {
        "phi_ase": "cm^-2 s^-1",
        "mse": "1",
        "total_rays": "count",
        "dndt_ase": "s^-1",
    }
    for name in ["phi_ase", "mse", "total_rays", "dndt_ase"]:
        record = iteration.meshes["core_result_" + name]
        assert record.get_attribute("haseTransportVersion") == HASE_TRANSPORT_VERSION
        assert record.get_attribute("haseEntity") == "point_level"
        assert _attribute_list(record.get_attribute("haseAxes")) == ["point", "level"]
        assert record.get_attribute("haseLayoutOrder") == "recordC"
        assert _attribute_list(record.get_attribute("hasePrimitiveShape")) == [5, 6]
        assert record.get_attribute("haseUnit") == expected_units[name]
        assert list(record.axis_labels) == ["point", "level"]
    series.close()


def _openpmd_backend_values():
    raw = os.environ.get("HASE_OPENPMD_TEST_BACKENDS")
    if raw:
        return [backend.strip() for backend in raw.split(",") if backend.strip()]
    return ["adios"]


def _cache_value(cache_path, name):
    if cache_path is None or not cache_path.is_file():
        return None
    for line in cache_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(name + ":"):
            return line.split("=", 1)[1].strip()
    return None


def _build_dir_for_executable(executable):
    path = Path(executable).resolve()
    for parent in [path.parent, *path.parents]:
        if (parent / "CMakeCache.txt").is_file():
            return parent
        if parent == Path.cwd().resolve():
            break

    for build_dir in _build_dir_candidates():
        if (build_dir / "CMakeCache.txt").is_file():
            return build_dir
    return None


def _target_uses_openpmd_main(build_dir):
    if build_dir is None:
        return True
    manifests = [build_dir / name for name in ("build.ninja", "Makefile", "compile_commands.json")]
    existing = [path for path in manifests if path.is_file()]
    if not existing:
        return True
    return any("src/openpmd_main.cpp" in path.read_text(encoding="utf-8", errors="ignore") for path in existing)


def _installed_calc_phi_ase_candidates():
    try:
        import HASEonGPU_Bindings
    except ImportError:
        return []

    return [Path(path) / "calcPhiASE" for path in HASEonGPU_Bindings.__path__]


def _calc_phi_ase_candidates():
    env = os.environ.get("HASE_OPENPMD_CALCPHIASE") or os.environ.get("HASE_CALCPHIASE")
    if env:
        yield Path(env)
    yield from _installed_calc_phi_ase_candidates()
    for build_dir in _build_dir_candidates():
        yield build_dir / "python" / "HASEonGPU_Bindings" / "calcPhiASE"
        yield build_dir / "calcPhiASE"


def _openpmd_transport_sources():
    root = Path(__file__).resolve().parents[3]
    yield root / "src" / "openpmd_main.cpp"
    yield root / "src" / "openpmd" / "OpenPmdParser.cpp"
    yield root / "include" / "openpmd" / "OpenPmdParser.hpp"


def _is_current_openpmd_calc_phi_ase(executable):
    build_dir = _build_dir_for_executable(executable)
    if not _target_uses_openpmd_main(build_dir):
        return False
    executable_mtime = executable.stat().st_mtime
    return all(
        source.stat().st_mtime <= executable_mtime
        for source in _openpmd_transport_sources()
        if source.is_file()
    )


def _openpmd_calc_phi_ase():
    configured = os.environ.get("HASE_OPENPMD_CALCPHIASE") or os.environ.get("HASE_CALCPHIASE")
    for executable in _calc_phi_ase_candidates():
        if (
            executable.is_file()
            and os.access(executable, os.X_OK)
            and _is_current_openpmd_calc_phi_ase(executable)
        ):
            return executable.resolve(), _build_dir_for_executable(executable)
    if configured:
        pytest.fail(f"configured openPMD calcPhiASE is missing, not executable, or stale: {configured}")
    pytest.skip("no current openPMD calcPhiASE binary found; build calcPhiASE from src/openpmd_main.cpp or set HASE_CALCPHIASE")


@pytest.mark.integration
def testAdiosSstWatchdogEmitsFiveAliveBeatsFromNormalSession(monkeypatch):
    configuredBackends = _openpmd_backend_values()
    if "adios-sst" not in configuredBackends:
        pytest.skip("watchdog beat test only runs when HASE_OPENPMD_TEST_BACKENDS includes adios-sst")

    io = _io()
    if "sst" not in getattr(io, "file_extensions", []):
        pytest.skip("openPMD-api was not built with ADIOS2 SST support")

    executable, _ = _openpmd_calc_phi_ase()
    monkeypatch.setenv("HASE_CALCPHIASE", str(executable))

    expectedBeats = 5
    aliveBeats = 0
    with transport.OpenPmdPhiAseSession(transport="adios-sst", watchdog_interval=0.05) as session:
        while aliveBeats < expectedBeats:
            try:
                isAlive, payload = session._watchdog_events.get(timeout=1.0)
            except transport.queue.Empty:
                pytest.fail(f"watchdog emitted only {aliveBeats} alive beats")
            if not isAlive:
                raise AssertionError("watchdog reported backend failure") from payload
            aliveBeats += 1

    assert aliveBeats == expectedBeats


def _mpi_enabled(build_dir):
    cache = None if build_dir is None else build_dir / "CMakeCache.txt"
    if _cache_value(cache, "DISABLE_MPI") == "ON":
        return False
    compile_commands = None if build_dir is None else build_dir / "compile_commands.json"
    if compile_commands is not None and compile_commands.is_file():
        commands = compile_commands.read_text(encoding="utf-8", errors="ignore")
        if "src/openpmd_main.cpp" in commands and "-DDISABLE_MPI" in commands:
            return False
    mpi_found = _cache_value(cache, "MPI_CXX_FOUND") or _cache_value(cache, "MPI_FOUND")
    return mpi_found is None or mpi_found.upper() in {"TRUE", "ON", "1", "YES"}


def _assert_mpi_mode_rejected_by_non_mpi_build(tmp_path, executable):
    _require_openpmd_transport_io()
    phi_ase = launch_smoke_phi_ase()
    phi_ase.parallelMode = "mpi"
    input_path = tmp_path / f"mpi_input{_file_suffix_for_tests()}"
    output_path = tmp_path / f"mpi_output{_file_suffix_for_tests()}"

    with transport.OpenPmdInputSeries(input_path, backend=_file_backend_for_tests()) as writer:
        writer.write(phi_ase, launch_smoke_medium(), launch_smoke_cross_sections())
    completed = subprocess.run(
        [str(executable), f"--input-path={input_path}", f"--output-path={output_path}"],
        check=False,
        text=True,
        capture_output=True,
        timeout=90,
    )
    assert completed.returncode != 0


def _matrix_mpi_rank_count():
    value = os.environ.get("HASE_MPI_TEST_RANKS", "").strip()
    if not value:
        pytest.skip("HASE_MPI_TEST_RANKS is not configured for this CI row")
    if "," in value:
        pytest.fail("HASE_MPI_TEST_RANKS must contain exactly one rank count per CI row")
    try:
        ranks = int(value)
    except ValueError:
        pytest.fail(f"HASE_MPI_TEST_RANKS must be an integer, got {value!r}")
    if ranks < 1 or ranks > 4:
        pytest.fail(f"HASE_MPI_TEST_RANKS must be between 1 and 4, got {ranks}")
    return ranks


def _mpiexec_command_prefix(ranks):
    return ["mpiexec", *shlex.split(os.environ.get("HASE_MPIEXEC_EXTRA_ARGS", "")), "-n", str(ranks)]


def _assert_launch_smoke_result(result):
    expected_size = launch_smoke_topology().numberOfPoints * launch_smoke_topology().levels
    assert result.phiAse.shape == (expected_size,)
    assert result.mse.shape == (expected_size,)
    assert result.totalRays.shape == (expected_size,)
    assert result.dndtAse.shape == (expected_size,)
    assert np.all(np.isfinite(result.phiAse))
    assert np.all(np.isfinite(result.mse))
    assert np.all(result.totalRays >= 0)


def _round_trip_calc_phi_ase(tmp_path, parallel_mode):
    _require_openpmd_transport_io()
    executable, build_dir = _openpmd_calc_phi_ase()
    if parallel_mode == "mpi" and not _mpi_enabled(build_dir):
        return _assert_mpi_mode_rejected_by_non_mpi_build(tmp_path, executable)

    phi_ase = launch_smoke_phi_ase()
    phi_ase.parallelMode = parallel_mode
    input_path = tmp_path / f"{parallel_mode}_input{_file_suffix_for_tests()}"
    output_path = tmp_path / f"{parallel_mode}_output{_file_suffix_for_tests()}"

    with transport.OpenPmdInputSeries(input_path, backend=_file_backend_for_tests()) as writer:
        writer.write(phi_ase, launch_smoke_medium(), launch_smoke_cross_sections())
    completed = subprocess.run(
        [str(executable), f"--input-path={input_path}", f"--output-path={output_path}"],
        check=False,
        text=True,
        capture_output=True,
        timeout=90,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr

    result = transport.read_result(output_path)
    _assert_launch_smoke_result(result)

    io = _io()
    series = io.Series(str(input_path), io.Access.read_only)
    try:
        assert series.iterations[0].get_attribute("parallel_mode") == parallel_mode
    finally:
        series.close()


@pytest.mark.integration
@pytest.mark.parametrize("openpmd_backend", _openpmd_backend_values())
def test_python_api_launches_configured_openpmd_backend_once(monkeypatch, openpmd_backend):
    _require_openpmd_transport_io()
    executable, _ = _openpmd_calc_phi_ase()
    monkeypatch.setenv("HASE_CALCPHIASE", str(executable))
    phi_ase = launch_smoke_phi_ase()
    phi_ase.openpmdBackend = openpmd_backend
    phi_ase.run(gainMedium=launch_smoke_medium(), crossSections=launch_smoke_cross_sections())
    result = phi_ase.getResults()

    _assert_launch_smoke_result(result)


@pytest.mark.integration
def test_calc_phi_ase_single_openpmd_round_trip(tmp_path):
    _round_trip_calc_phi_ase(tmp_path, "single")


@pytest.mark.integration
def test_calc_phi_ase_mpi_openpmd_round_trip_when_mpi_enabled(tmp_path):
    _round_trip_calc_phi_ase(tmp_path, "mpi")


@pytest.mark.integration
def test_calc_phi_ase_mpi_openpmd_round_trip_with_matrix_rank(monkeypatch, tmp_path):
    _require_openpmd_transport_io()
    ranks = _matrix_mpi_rank_count()
    executable, build_dir = _openpmd_calc_phi_ase()
    if not _mpi_enabled(build_dir):
        pytest.fail("HASE_MPI_TEST_RANKS is configured, but calcPhiASE was not built with MPI")

    monkeypatch.setenv("HASE_CALCPHIASE", str(executable))
    phi_ase = launch_smoke_phi_ase()
    phi_ase.parallelMode = "mpi"

    result = transport.runPhiASE(
        phi_ase,
        launch_smoke_medium(),
        launch_smoke_cross_sections(),
        transport=_file_backend_for_tests(),
        command_prefix=_mpiexec_command_prefix(ranks),
        workspace_dir=tmp_path,
    )

    _assert_launch_smoke_result(result)


def test_openpmd_input_series_preserves_custom_fields(tmp_path):
    _require_openpmd_transport_io()

    temperature_dimension = unitDimension.tFluo

    class ThermalPrism(PrismSchema):
        temperature = PrimitiveFieldSpec(
            "temperature",
            "custom_temperature",
            np.float64,
            unit="K",
            unitDimension=temperature_dimension,
            backendRequired=False,
        )

    values = np.array(
        [
            [300.0, 301.0, 302.0, 303.0, 304.0],
            [310.0, 311.0, 312.0, 313.0, 314.0],
            [320.0, 321.0, 322.0, 323.0, 324.0],
        ],
        dtype=np.float64,
    )
    medium = asymmetric_medium().withPrimitiveSchema(ThermalPrism, temperature=values)
    assert next(iter(medium.getPrisms())).temperature == pytest.approx(300.0)
    output = tmp_path / ("custom" + _file_suffix_for_tests())

    with transport.OpenPmdInputSeries(output, backend=_file_backend_for_tests()) as writer:
        writer.write(asymmetric_phi_ase(), medium, asymmetric_cross_sections())

    io = _io()
    series = io.Series(str(output), io.Access.read_only)
    iteration = series.iterations[0]
    record = iteration.meshes["custom_temperature"]
    np.testing.assert_array_equal(_read_scalar(series, iteration, "custom_temperature"), values.reshape(-1, order="F"))
    assert record.get_attribute("haseTransportVersion") == HASE_TRANSPORT_VERSION
    assert record.get_attribute("haseEntity") == "cell_layer"
    assert _attribute_list(record.get_attribute("haseAxes")) == ["cell", "layer"]
    assert record.get_attribute("haseAxesString") == "cell,layer"
    assert record.get_attribute("haseLayoutOrder") == "backendFlat"
    assert _attribute_list(record.get_attribute("hasePrimitiveShape")) == [3, 5]
    assert record.get_attribute("haseStatic") is True
    assert record.get_attribute("haseDynamic") is False
    assert record.get_attribute("haseBackendRequired") is False
    assert record.get_attribute("haseUserDefined") is True
    assert record.get_attribute("haseUserFieldName") == "temperature"
    assert record.get_attribute("haseUnit") == "K"
    _assert_base_openpmd_scalar_metadata(record, axis_labels=["flatIndex"], unit_dimension=temperature_dimension)
    series.close()
