include_guard(GLOBAL)

set(HASE_BUILD_OPENPMD_FROM_SOURCE_DESCRIPTION
    "Fetch and build the pinned openPMD-api C++ library and Python bindings instead of using an externally installed openPMD-api package."
)
option(
    HASE_BUILD_OPENPMD_FROM_SOURCE
    "${HASE_BUILD_OPENPMD_FROM_SOURCE_DESCRIPTION}"
    OFF
)
if(DEFINED CACHE{HASE_USE_SYSTEM_OPENPMD})
    set(HASE_LEGACY_USE_SYSTEM_OPENPMD "${HASE_USE_SYSTEM_OPENPMD}")
    message(
        DEPRECATION
        "HASE_USE_SYSTEM_OPENPMD is deprecated. Use "
        "HASE_BUILD_OPENPMD_FROM_SOURCE instead. "
        "HASE_USE_SYSTEM_OPENPMD=ON maps to "
        "HASE_BUILD_OPENPMD_FROM_SOURCE=OFF; "
        "HASE_USE_SYSTEM_OPENPMD=OFF maps to "
        "HASE_BUILD_OPENPMD_FROM_SOURCE=ON."
    )
    unset(HASE_USE_SYSTEM_OPENPMD CACHE)
    if(HASE_LEGACY_USE_SYSTEM_OPENPMD)
        set(HASE_BUILD_OPENPMD_FROM_SOURCE
            OFF
            CACHE BOOL
            "${HASE_BUILD_OPENPMD_FROM_SOURCE_DESCRIPTION}"
            FORCE
        )
    else()
        set(HASE_BUILD_OPENPMD_FROM_SOURCE
            ON
            CACHE BOOL
            "${HASE_BUILD_OPENPMD_FROM_SOURCE_DESCRIPTION}"
            FORCE
        )
    endif()
endif()
if(DEFINED CACHE{HASE_OPENPMD_BACKEND})
    message(
        DEPRECATION
        "HASE_OPENPMD_BACKEND is no longer a CMake build option. "
        "The openPMD storage backend is selected at runtime through "
        "PhiASE.openpmdBackend, YAML openpmd_backend, or direct transport= "
        "arguments. The stale CMake cache entry will be ignored."
    )
    unset(HASE_OPENPMD_BACKEND CACHE)
endif()
if(HASE_BUILD_OPENPMD_FROM_SOURCE)
    set(HASE_USE_SYSTEM_OPENPMD OFF)
else()
    set(HASE_USE_SYSTEM_OPENPMD ON)
endif()

set(HASE_OPENPMD_FILE_EXTENSION "sst")
set(HASE_OPENPMD_TEST_FILE_EXTENSION "bp")

set(HASE_OPENPMD_GIT_REPOSITORY "https://github.com/openPMD/openPMD-api.git")
set(HASE_OPENPMD_GIT_TAG "0.17.0")
set(HASE_ADIOS2_GIT_REPOSITORY "https://github.com/ornladios/ADIOS2.git")
set(HASE_ADIOS2_GIT_TAG "v2.12.1")
set(HASE_HDF5_GIT_REPOSITORY "https://github.com/HDFGroup/hdf5.git")
set(HASE_HDF5_GIT_TAG "hdf5_1.14.6")

if(DEFINED openPMD_SUPERBUILD)
    set(HASE_OPENPMD_SUPERBUILD_DEFAULT ${openPMD_SUPERBUILD})
else()
    set(HASE_OPENPMD_SUPERBUILD_DEFAULT ON)
endif()
option(
    HASE_OPENPMD_SUPERBUILD
    "Allow openPMD-api to fetch/build its bundled helper dependencies"
    ${HASE_OPENPMD_SUPERBUILD_DEFAULT}
)
option(
    HASE_OPENPMD_BUILD_PYTHON_BINDINGS
    "Build openPMD-api Python bindings as part of the HASE CMake build tree"
    OFF
)
set(HASE_OPENPMD_PYTHON_PACKAGE_DIR
    ""
    CACHE PATH
    "Directory containing the openpmd_api Python package matching openPMD::openPMD"
)
if(NOT HASE_BUILD_OPENPMD_FROM_SOURCE)
    message(STATUS "Using external openPMD-api for the HASE openPMD transport")
    find_package(openPMD CONFIG QUIET)

    if(NOT openPMD_FOUND)
        message(
            FATAL_ERROR
            "HASEonGPU now expects an external openPMD-api C++ package by default, "
            "but CMake could not find openPMD. Install or build an external "
            "openPMD-api C++ package, then point CMake to it with openPMD_DIR "
            "or CMAKE_PREFIX_PATH. To use the bundled source-build path instead, "
            "configure with "
            "-DHASE_BUILD_OPENPMD_FROM_SOURCE=ON."
        )
    endif()

    if(NOT TARGET openPMD::openPMD)
        message(
            FATAL_ERROR
            "The external openPMD package did not provide openPMD::openPMD. "
            "Use an openPMD-api CMake installation that exports openPMD::openPMD, "
            "or configure with -DHASE_BUILD_OPENPMD_FROM_SOURCE=ON to use the "
            "bundled source-build path."
        )
    endif()

    if(HASE_OPENPMD_BUILD_PYTHON_BINDINGS)
        message(
            STATUS
            "HASE_OPENPMD_BUILD_PYTHON_BINDINGS is ignored with the default external openPMD contract; "
            "the Python openpmd_api module must come from the same external installation."
        )
    endif()
    set(HASE_OPENPMD_BUILD_PYTHON_BINDINGS
        OFF
        CACHE BOOL
        "Build openPMD-api Python bindings as part of the HASE CMake build"
        FORCE
    )

    if(HASE_ENABLE_PYTHON)
        if(HASE_OPENPMD_PYTHON_PACKAGE_DIR)
            if(NOT EXISTS "${HASE_OPENPMD_PYTHON_PACKAGE_DIR}/openpmd_api")
                message(
                    FATAL_ERROR
                    "HASE_OPENPMD_PYTHON_PACKAGE_DIR='${HASE_OPENPMD_PYTHON_PACKAGE_DIR}' "
                    "does not contain an openpmd_api Python package. Set this cache variable "
                    "to the site-packages directory from the same external openPMD-api installation "
                    "as openPMD::openPMD."
                )
            endif()
            message(
                STATUS
                "HASE external openPMD Python package directory: ${HASE_OPENPMD_PYTHON_PACKAGE_DIR}"
            )
        else()
            message(
                STATUS
                "HASE_ENABLE_PYTHON=ON with the external openPMD contract will use "
                "the runtime Python environment's openpmd_api package. It must match "
                "the external openPMD::openPMD C++ provider and support the "
                "runtime openPMD backend selected by Python or YAML. Set "
                "-DHASE_OPENPMD_PYTHON_PACKAGE_DIR=<site-packages directory containing openpmd_api> "
                "only when you need HASEonGPU to prefer a specific provider path."
            )
        endif()
    endif()
    return()
endif()

message(STATUS "Fetching pinned openPMD-api for the HASE openPMD transport")

include(FetchContent)

set(HASE_OPENPMD_USE_ADIOS2 ON)
set(HASE_OPENPMD_USE_HDF5 ON)
set(HASE_OPENPMD_USE_SST ON)

if(HASE_OPENPMD_USE_HDF5 AND HASE_OPENPMD_SUPERBUILD)
    message(STATUS "Fetching pinned HDF5 for the HASE openPMD transport")
    set(BUILD_SHARED_LIBS
        ON
        CACHE BOOL
        "Build shared third-party libraries for the HASE openPMD transport"
        FORCE
    )
    set(HDF5_BUILD_CPP_LIB
        OFF
        CACHE BOOL
        "Disable HDF5 C++ bindings in the HASE superbuild"
        FORCE
    )
    set(HDF5_BUILD_FORTRAN
        OFF
        CACHE BOOL
        "Disable HDF5 Fortran bindings in the HASE superbuild"
        FORCE
    )
    set(HDF5_BUILD_HL_LIB
        OFF
        CACHE BOOL
        "Disable the HDF5 high-level library in the HASE superbuild"
        FORCE
    )
    set(HDF5_BUILD_JAVA
        OFF
        CACHE BOOL
        "Disable HDF5 Java bindings in the HASE superbuild"
        FORCE
    )
    set(HDF5_BUILD_TOOLS
        OFF
        CACHE BOOL
        "Disable HDF5 command-line tools in the HASE superbuild"
        FORCE
    )
    set(HDF5_BUILD_EXAMPLES
        OFF
        CACHE BOOL
        "Disable HDF5 examples in the HASE superbuild"
        FORCE
    )
    set(HDF5_ENABLE_SZIP_SUPPORT
        OFF
        CACHE BOOL
        "Disable optional SZIP support in the HASE HDF5 superbuild"
        FORCE
    )
    set(HDF5_ENABLE_Z_LIB_SUPPORT
        OFF
        CACHE BOOL
        "Disable optional zlib support in the HASE HDF5 superbuild"
        FORCE
    )
    set(BUILD_TESTING
        OFF
        CACHE BOOL
        "Disable third-party tests while configuring the HASE superbuild dependencies"
        FORCE
    )
    if(MPI_FOUND)
        find_package(MPI COMPONENTS C REQUIRED)
        set(HDF5_ENABLE_PARALLEL
            ON
            CACHE BOOL
            "Enable parallel HDF5 when HASE MPI is available"
            FORCE
        )
        set(HASE_INTERNAL_HDF5_IS_PARALLEL TRUE)
    else()
        set(HDF5_ENABLE_PARALLEL
            OFF
            CACHE BOOL
            "Disable parallel HDF5 when HASE MPI is unavailable"
            FORCE
        )
        set(HASE_INTERNAL_HDF5_IS_PARALLEL FALSE)
    endif()

    FetchContent_Declare(
        HDF5
        GIT_REPOSITORY "${HASE_HDF5_GIT_REPOSITORY}"
        GIT_TAG "${HASE_HDF5_GIT_TAG}"
    )
    FetchContent_MakeAvailable(HDF5)

    set(HASE_INTERNAL_HDF5_INCLUDE_DIRS
        "${hdf5_SOURCE_DIR}/src"
        "${hdf5_SOURCE_DIR}/src/H5FDsubfiling"
        "${hdf5_BINARY_DIR}/src"
    )
    set(HASE_INTERNAL_FIND_MODULE_DIR
        "${CMAKE_BINARY_DIR}/hase-cmake-overrides"
    )
    file(MAKE_DIRECTORY "${HASE_INTERNAL_FIND_MODULE_DIR}")
    file(
        WRITE
        "${HASE_INTERNAL_FIND_MODULE_DIR}/FindHDF5.cmake"
        "if(NOT TARGET hdf5-shared)\n"
        "    message(FATAL_ERROR \"HASE internal HDF5 target hdf5-shared is not available\")\n"
        "endif()\n"
        "set(HDF5_FOUND TRUE)\n"
        "set(HDF5_C_FOUND TRUE)\n"
        "set(HDF5_VERSION \"1.14.6\")\n"
        "set(HDF5_LIBRARIES hdf5-shared)\n"
        "set(HDF5_C_LIBRARIES hdf5-shared)\n"
        "set(HDF5_INCLUDE_DIRS \"${HASE_INTERNAL_HDF5_INCLUDE_DIRS}\")\n"
        "set(HDF5_C_INCLUDE_DIRS \"${HASE_INTERNAL_HDF5_INCLUDE_DIRS}\")\n"
        "set(HDF5_DEFINITIONS \"\")\n"
        "set(HDF5_IS_PARALLEL ${HASE_INTERNAL_HDF5_IS_PARALLEL})\n"
        "set(HDF5_ENABLE_PARALLEL ${HASE_INTERNAL_HDF5_IS_PARALLEL})\n"
        "set(HDF5_PROVIDES_PARALLEL ${HASE_INTERNAL_HDF5_IS_PARALLEL})\n"
        "foreach(component IN LISTS HDF5_FIND_COMPONENTS)\n"
        "    if(component STREQUAL \"C\")\n"
        "        set(HDF5_C_FOUND TRUE)\n"
        "    else()\n"
        "        set(HDF5_${component}_FOUND FALSE)\n"
        "    endif()\n"
        "endforeach()\n"
    )
    list(PREPEND CMAKE_MODULE_PATH "${HASE_INTERNAL_FIND_MODULE_DIR}")
endif()

if(HASE_OPENPMD_USE_ADIOS2)
    message(STATUS "Fetching pinned ADIOS2 for the HASE openPMD transport")
    set(ADIOS2_USE_Fortran
        OFF
        CACHE BOOL
        "Disable ADIOS2 Fortran bindings in the HASE superbuild"
        FORCE
    )
    set(ADIOS2_USE_Python
        OFF
        CACHE BOOL
        "Disable ADIOS2 Python bindings in the HASE superbuild"
        FORCE
    )
    set(ADIOS2_BUILD_EXAMPLES
        OFF
        CACHE BOOL
        "Disable ADIOS2 examples in the HASE superbuild"
        FORCE
    )
    set(ADIOS2_BUILD_TESTING
        OFF
        CACHE BOOL
        "Disable ADIOS2 tests in the HASE superbuild"
        FORCE
    )
    set(BUILD_TESTING
        OFF
        CACHE BOOL
        "Disable third-party tests while configuring the HASE superbuild dependencies"
        FORCE
    )
    set(ADIOS2_INSTALL_GENERATE_CONFIG
        OFF
        CACHE BOOL
        "Disable ADIOS2's install-time adios2-config helper generation in the HASE superbuild"
        FORCE
    )
    set(ADIOS2_USE_SST
        ${HASE_OPENPMD_USE_SST}
        CACHE STRING
        "Enable ADIOS2 SST for the HASE openPMD transport"
        FORCE
    )

    # Keep the ADIOS2 superbuild narrow. HASE's openPMD transport uses ADIOS2
    # ADIOS/SST-style openPMD series and does not need HDF5, compression plugins,
    # remote/cloud transports, visualization hooks, or profiling infrastructure.
    foreach(
        HASE_ADIOS2_DISABLED_OPTION
        IN
        ITEMS
            BZip2
            Blosc2
            Campaign
            Catalyst
            DAOS
            DataMan
            DataSpaces
            Endian_Reverse
            HDF5
            HDF5_VOL
            IME
            LIBPRESSIO
            MGARD
            MHS
            PNG
            Profiling
            SZ
            Sodium
            SysVShMem
            UCX
            ZFP
            ZeroMQ
    )
        set(ADIOS2_USE_${HASE_ADIOS2_DISABLED_OPTION}
            OFF
            CACHE STRING
            "Disable unused ADIOS2 ${HASE_ADIOS2_DISABLED_OPTION} support in the HASE superbuild"
            FORCE
        )
    endforeach()
    if(MPI_FOUND)
        set(ADIOS2_USE_MPI
            ON
            CACHE BOOL
            "Enable MPI in ADIOS2 when HASE MPI is available"
            FORCE
        )
    else()
        set(ADIOS2_USE_MPI
            OFF
            CACHE BOOL
            "Disable MPI in ADIOS2 when HASE MPI is unavailable"
            FORCE
        )
    endif()

    FetchContent_Declare(
        ADIOS2
        GIT_REPOSITORY "${HASE_ADIOS2_GIT_REPOSITORY}"
        GIT_TAG "${HASE_ADIOS2_GIT_TAG}"
        OVERRIDE_FIND_PACKAGE
    )
    FetchContent_MakeAvailable(ADIOS2)
    if(EXISTS "${ADIOS2_BINARY_DIR}/adios2-config.cmake")
        set(ADIOS2_DIR
            "${ADIOS2_BINARY_DIR}"
            CACHE PATH
            "ADIOS2 CMake config directory produced by the HASE FetchContent build"
            FORCE
        )
    endif()
    if(NOT DEFINED ADIOS2_VERSION OR "${ADIOS2_VERSION}" STREQUAL "")
        string(REGEX REPLACE "^v" "" ADIOS2_VERSION "${HASE_ADIOS2_GIT_TAG}")
        set(ADIOS2_VERSION
            "${ADIOS2_VERSION}"
            CACHE STRING
            "ADIOS2 version provided by the HASE FetchContent build"
            FORCE
        )
    endif()
endif()

set(openPMD_USE_ADIOS2
    ${HASE_OPENPMD_USE_ADIOS2}
    CACHE STRING
    "Enable ADIOS2 backend for the HASE openPMD transport"
    FORCE
)
set(openPMD_USE_HDF5
    ${HASE_OPENPMD_USE_HDF5}
    CACHE STRING
    "Enable HDF5 backend for the HASE openPMD transport"
    FORCE
)
set(openPMD_HAVE_PKGCONFIG
    OFF
    CACHE BOOL
    "Do not generate pkg-config metadata in the HASE superbuild"
    FORCE
)
set(openPMD_USE_VERIFY
    OFF
    CACHE BOOL
    "Disable openPMD internal VERIFY checks in the HASE superbuild"
    FORCE
)
set(openPMD_SUPERBUILD
    ${HASE_OPENPMD_SUPERBUILD}
    CACHE BOOL
    "Allow openPMD-api to fetch/build its bundled helper dependencies"
    FORCE
)
foreach(HASE_OPENPMD_INTERNAL_DEP IN ITEMS CATCH JSON TOML11 PYBIND11)
    set(openPMD_USE_INTERNAL_${HASE_OPENPMD_INTERNAL_DEP}
        ${HASE_OPENPMD_SUPERBUILD}
        CACHE BOOL
        "Use openPMD-api bundled ${HASE_OPENPMD_INTERNAL_DEP} dependency"
        FORCE
    )
endforeach()
set(openPMD_USE_PYTHON
    ${HASE_OPENPMD_BUILD_PYTHON_BINDINGS}
    CACHE BOOL
    "Build openPMD-api Python bindings from the HASE CMake build"
    FORCE
)
set(openPMD_BUILD_TESTING
    OFF
    CACHE BOOL
    "Disable openPMD-api tests in the HASE superbuild"
    FORCE
)
set(openPMD_BUILD_EXAMPLES
    OFF
    CACHE BOOL
    "Disable openPMD-api examples in the HASE superbuild"
    FORCE
)
set(openPMD_BUILD_CLI_TOOLS
    OFF
    CACHE BOOL
    "Disable openPMD-api CLI tools in the HASE superbuild"
    FORCE
)
set(openPMD_INSTALL
    OFF
    CACHE BOOL
    "Do not install openPMD-api from the HASE superbuild"
    FORCE
)
if(MPI_FOUND)
    set(openPMD_USE_MPI
        ON
        CACHE STRING
        "Enable MPI in openPMD-api when HASE MPI is available"
        FORCE
    )
else()
    set(openPMD_USE_MPI
        OFF
        CACHE STRING
        "Disable MPI in openPMD-api when HASE MPI is unavailable"
        FORCE
    )
endif()

FetchContent_Declare(
    openPMD
    GIT_REPOSITORY "${HASE_OPENPMD_GIT_REPOSITORY}"
    GIT_TAG "${HASE_OPENPMD_GIT_TAG}"
)
FetchContent_MakeAvailable(openPMD)

if(NOT TARGET openPMD::openPMD)
    message(FATAL_ERROR "openPMD::openPMD target was not created")
endif()

if(TARGET openPMD)
    set_target_properties(
        openPMD
        PROPERTIES INSTALL_RPATH "${HASE_INSTALL_LIB_RPATH}"
    )
endif()
if(TARGET hdf5-shared)
    set_target_properties(
        hdf5-shared
        PROPERTIES INSTALL_RPATH "${HASE_INSTALL_LIB_RPATH}"
    )
endif()

if(TARGET openPMD.py)
    set_target_properties(
        openPMD.py
        PROPERTIES
            INSTALL_RPATH "${HASE_INSTALL_RPATH}"
            INSTALL_RPATH_USE_LINK_PATH ON
    )
    add_custom_target(hase_openpmd_python DEPENDS openPMD.py)
endif()
