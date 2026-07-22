"""Discover the durable native runtime used by the thin Python frontend."""

from __future__ import annotations

import os
from pathlib import Path
from types import ModuleType


def _path_entries(value):
    for entry in str(value or "").split(os.pathsep):
        if entry:
            yield Path(entry).expanduser()


def _packaged_runtime_dir():
    try:
        from . import _config
    except ImportError:
        return ""
    return str(getattr(_config, "HASE_RUNTIME_DIR", "") or "")


def _source_root():
    candidate = Path(__file__).resolve().parents[2]
    if (candidate / "CMakeLists.txt").is_file() and (candidate / "pyproject.toml").is_file():
        return candidate
    return None


def _unique(paths):
    seen = set()
    for path in paths:
        path = Path(path).expanduser()
        if not path.is_absolute():
            source_root = _source_root()
            path = (source_root if source_root is not None else Path.cwd()) / path
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        yield resolved


def runtime_roots():
    """Yield runtime roots in explicit-to-default precedence order."""
    paths = [
        *_path_entries(os.environ.get("HASE_RUNTIME_DIR")),
        *_path_entries(_packaged_runtime_dir()),
    ]
    source_root = _source_root()
    if source_root is not None:
        paths.append(source_root / "build")
    yield from _unique(paths)


def runtime_root():
    """Return the selected runtime root, even before it has been built."""
    try:
        return next(runtime_roots())
    except StopIteration as exc:
        raise RuntimeError(
            "The HASE Python frontend has no configured native runtime. Reinstall with "
            "`python3 -m pip install .` or set HASE_RUNTIME_DIR."
        ) from exc


def runtime_metadata_path(root=None):
    root = runtime_root() if root is None else Path(root).expanduser().resolve()
    return root / "python" / "pyInclude" / "_runtime" / "_config.py"


def runtime_config():
    """Load live provider metadata from the selected native runtime."""
    root = runtime_root()
    path = runtime_metadata_path(root)
    if not path.is_file():
        raise RuntimeError(
            f"HASE runtime metadata is missing: {path}. Configure/build the selected "
            f"runtime '{root}' before launching HASE."
        )
    module = ModuleType(f"pyInclude._runtime._live_config_{abs(hash(str(path)))}")
    module.__file__ = str(path)
    exec(compile(path.read_bytes(), str(path), "exec"), module.__dict__)
    recorded_root_text = str(getattr(module, "HASE_RUNTIME_DIR", "") or "")
    if not recorded_root_text or Path(recorded_root_text).expanduser().resolve() != root:
        raise RuntimeError(
            f"HASE runtime metadata {path} does not describe the selected runtime '{root}'. "
            "Reconfigure that runtime before launching HASE."
        )
    return module


def native_dirs_from_root(root):
    root = Path(root).expanduser().resolve()
    yield root
    yield root / "bin"
    yield root / "lib"
    yield root / "python" / "pyInclude" / "_runtime"


def runtime_native_dirs():
    yield from _unique(native_dirs_from_root(runtime_root()))


def runtime_executable_candidates(names):
    for directory in runtime_native_dirs():
        for name in names:
            yield directory / name


def runtime_library_candidates(names):
    root = runtime_root()
    directories = (
        root / "lib",
        root / "bin",
        root,
        root / "python" / "pyInclude" / "_runtime",
    )
    for directory in _unique(directories):
        for name in names:
            yield directory / name
