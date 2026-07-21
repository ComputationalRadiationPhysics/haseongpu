/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#define HASEONGPU_VERSION_MAJOR 2
#define HASEONGPU_VERSION_MINOR 1
#define HASEONGPU_VERSION_PATCH 1

#define HASEONGPU_PP_STRINGIFY_DETAIL(value) #value
#define HASEONGPU_PP_STRINGIFY(value) HASEONGPU_PP_STRINGIFY_DETAIL(value)

#define HASEONGPU_VERSION_STRING                                                                                      \
    HASEONGPU_PP_STRINGIFY(HASEONGPU_VERSION_MAJOR)                                                                   \
    "." HASEONGPU_PP_STRINGIFY(HASEONGPU_VERSION_MINOR) "." HASEONGPU_PP_STRINGIFY(HASEONGPU_VERSION_PATCH)

/* version number encoding
 * 4 digits for major version (max 9999)
 * 3 digits for minor version (max 999)
 * 5 digits for patch version (max 99999)
 * example: version 1.2.3 -> 0001 002 00003
 */
#define HASEONGPU_VERSION_NUMBER(major, minor, patch)                                                                 \
    ((((major) % 10000llu) * 100'000'000llu) + (((minor) % 1000llu) * 100000llu) + ((patch) % 100000llu))

//! The HASEonGPU version number
#define HASEONGPU_VERSION                                                                                             \
    HASEONGPU_VERSION_NUMBER(HASEONGPU_VERSION_MAJOR, HASEONGPU_VERSION_MINOR, HASEONGPU_VERSION_PATCH)
