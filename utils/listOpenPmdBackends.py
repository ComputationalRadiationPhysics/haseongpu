#!/usr/bin/env python3

# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from HASEonGPU import OpenPmdBackends


def available_backends():
    return OpenPmdBackends.all()


def main():
    for backend in available_backends():
        print(backend)


if __name__ == "__main__":
    main()
