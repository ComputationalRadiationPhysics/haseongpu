#!/usr/bin/env python3
"""Print openPMD backends supported by the installed HASE frontend."""

from __future__ import annotations

import argparse

from HASEonGPU import OpenPmdBackends


def available_backends() -> list[str]:
    return OpenPmdBackends.all()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", action="store_true", help="print a comma-separated backend list")
    args = parser.parse_args(argv)

    backends = available_backends()
    if args.csv:
        print(",".join(backends))
    else:
        for backend in backends:
            print(backend)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
