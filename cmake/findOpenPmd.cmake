include_guard(GLOBAL)

set(HASE_OPENPMD_PROVIDER_DESCRIPTION
    "openPMD-api provider selection: auto, bundled, or system"
)
set(HASE_OPENPMD_PROVIDER_WAS_SET FALSE)
if(DEFINED CACHE{HASE_OPENPMD_PROVIDER})
    set(HASE_OPENPMD_PROVIDER_WAS_SET TRUE)
else()
    set(HASE_OPENPMD_PROVIDER
        "auto"
        CACHE STRING
        "${HASE_OPENPMD_PROVIDER_DESCRIPTION}"
    )
endif()
set_property(CACHE HASE_OPENPMD_PROVIDER PROPERTY STRINGS auto bundled system)

if(NOT HASE_OPENPMD_PROVIDER_WAS_SET)
    if(DEFINED CACHE{HASE_BUILD_OPENPMD_FROM_SOURCE})
        message(
            DEPRECATION
            "HASE_BUILD_OPENPMD_FROM_SOURCE is deprecated. Use "
            "HASE_OPENPMD_PROVIDER=bundled instead of ON, or "
            "HASE_OPENPMD_PROVIDER=system instead of OFF."
        )
        if(HASE_BUILD_OPENPMD_FROM_SOURCE)
            set(HASE_OPENPMD_PROVIDER "bundled")
        else()
            set(HASE_OPENPMD_PROVIDER "system")
        endif()
        set(HASE_OPENPMD_PROVIDER
            "${HASE_OPENPMD_PROVIDER}"
            CACHE STRING
            "${HASE_OPENPMD_PROVIDER_DESCRIPTION}"
            FORCE
        )
    elseif(DEFINED CACHE{HASE_USE_SYSTEM_OPENPMD})
        message(
            DEPRECATION
            "HASE_USE_SYSTEM_OPENPMD is deprecated. Use "
            "HASE_OPENPMD_PROVIDER=system instead of ON, or "
            "HASE_OPENPMD_PROVIDER=bundled instead of OFF."
        )
        if(HASE_USE_SYSTEM_OPENPMD)
            set(HASE_OPENPMD_PROVIDER "system")
        else()
            set(HASE_OPENPMD_PROVIDER "bundled")
        endif()
        set(HASE_OPENPMD_PROVIDER
            "${HASE_OPENPMD_PROVIDER}"
            CACHE STRING
            "${HASE_OPENPMD_PROVIDER_DESCRIPTION}"
            FORCE
        )
    endif()
endif()

string(TOLOWER "${HASE_OPENPMD_PROVIDER}" HASE_OPENPMD_PROVIDER)
if(NOT HASE_OPENPMD_PROVIDER MATCHES "^(auto|bundled|system)$")
    message(
        FATAL_ERROR
        "Invalid HASE_OPENPMD_PROVIDER='${HASE_OPENPMD_PROVIDER}'. "
        "Allowed values are: auto, bundled, system."
    )
endif()
set(HASE_OPENPMD_PROVIDER
    "${HASE_OPENPMD_PROVIDER}"
    CACHE STRING
    "${HASE_OPENPMD_PROVIDER_DESCRIPTION}"
    FORCE
)

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

option(
    HASE_OPENPMD_USE_ADIOS2
    "Enable ADIOS2 support in the HASE-managed openPMD provider"
    ON
)
option(
    HASE_OPENPMD_USE_HDF5
    "Enable HDF5 support in the HASE-managed openPMD provider"
    OFF
)
option(
    HASE_OPENPMD_USE_SST
    "Enable ADIOS2 SST support in the HASE-managed openPMD provider"
    ON
)
option(
    HASE_OPENPMD_FETCH_ADIOS2
    "Fetch and build pinned ADIOS2 for the HASE-managed openPMD provider"
    ON
)
option(
    HASE_OPENPMD_FETCH_HDF5
    "Fetch and build pinned HDF5 for the HASE-managed openPMD provider"
    ON
)

if(HASE_OPENPMD_USE_SST AND NOT HASE_OPENPMD_USE_ADIOS2)
    message(
        FATAL_ERROR
        "HASE_OPENPMD_USE_SST requires HASE_OPENPMD_USE_ADIOS2=ON"
    )
endif()
if(NOT HASE_OPENPMD_USE_ADIOS2 AND NOT HASE_OPENPMD_USE_HDF5)
    message(
        FATAL_ERROR
        "Enable at least one openPMD backend: HASE_OPENPMD_USE_ADIOS2 or HASE_OPENPMD_USE_HDF5"
    )
endif()
if(HASE_OPENPMD_USE_ADIOS2 AND HASE_OPENPMD_USE_SST)
    set(HASE_OPENPMD_FILE_EXTENSION "sst")
    set(HASE_OPENPMD_TEST_FILE_EXTENSION "bp")
elseif(HASE_OPENPMD_USE_ADIOS2)
    set(HASE_OPENPMD_FILE_EXTENSION "bp")
    set(HASE_OPENPMD_TEST_FILE_EXTENSION "bp")
else()
    set(HASE_OPENPMD_FILE_EXTENSION "h5")
    set(HASE_OPENPMD_TEST_FILE_EXTENSION "h5")
endif()

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
    "Build and install openPMD-api Python bindings with the HASE-managed provider"
    OFF
)
set(HASE_OPENPMD_PYTHON_PACKAGE_DIR
    ""
    CACHE PATH
    "Directory containing the openpmd_api Python package matching openPMD::openPMD"
)
set(HASE_OPENPMD_BUNDLED_PREFIX
    "${CMAKE_BINARY_DIR}/hase-openpmd-provider/install"
    CACHE PATH
    "Install prefix for the HASE-managed bundled openPMD-api provider"
)
set(HASE_OPENPMD_BUNDLED_BUILD_DIR
    "${CMAKE_BINARY_DIR}/hase-openpmd-provider/build"
    CACHE PATH
    "Build directory for the HASE-managed bundled openPMD-api provider"
)
option(
    HASE_OPENPMD_BUNDLED_REBUILD
    "Force rebuilding the HASE-managed bundled openPMD-api provider"
    OFF
)

function(hase_openpmd_validate_found provider_kind)
    if(NOT openPMD_FOUND)
        message(
            FATAL_ERROR
            "Internal error: openPMD was not found for provider '${provider_kind}'."
        )
    endif()
    if(NOT TARGET openPMD::openPMD)
        message(
            FATAL_ERROR
            "The ${provider_kind} openPMD package did not provide openPMD::openPMD. "
            "Use an openPMD-api CMake installation that exports openPMD::openPMD."
        )
    endif()
endfunction()

function(hase_openpmd_configure_python provider_kind)
    if(
        HASE_OPENPMD_BUILD_PYTHON_BINDINGS
        AND "${provider_kind}" STREQUAL "system"
    )
        message(
            STATUS
            "HASE_OPENPMD_BUILD_PYTHON_BINDINGS is ignored with HASE_OPENPMD_PROVIDER=system; "
            "the Python openpmd_api module must come from the same system installation."
        )
        set(HASE_OPENPMD_BUILD_PYTHON_BINDINGS
            OFF
            CACHE BOOL
            "Build and install openPMD-api Python bindings with the HASE-managed provider"
            FORCE
        )
    endif()

    if(NOT HASE_ENABLE_PYTHON)
        return()
    endif()

    if(HASE_OPENPMD_PYTHON_PACKAGE_DIR)
        if(NOT EXISTS "${HASE_OPENPMD_PYTHON_PACKAGE_DIR}/openpmd_api")
            message(
                FATAL_ERROR
                "HASE_OPENPMD_PYTHON_PACKAGE_DIR='${HASE_OPENPMD_PYTHON_PACKAGE_DIR}' "
                "does not contain an openpmd_api Python package. Set this cache variable "
                "to the site-packages directory from the same openPMD-api installation "
                "as openPMD::openPMD."
            )
        endif()
        message(
            STATUS
            "HASE ${provider_kind} openPMD Python package directory: ${HASE_OPENPMD_PYTHON_PACKAGE_DIR}"
        )
    elseif("${provider_kind}" STREQUAL "system")
        message(
            STATUS
            "HASE_ENABLE_PYTHON=ON with HASE_OPENPMD_PROVIDER=system will use "
            "the runtime Python environment's openpmd_api package. It must match "
            "the system openPMD::openPMD C++ provider and support the runtime "
            "openPMD backend selected by Python or YAML. Set "
            "-DHASE_OPENPMD_PYTHON_PACKAGE_DIR=<site-packages directory containing openpmd_api> "
            "only when you need HASEonGPU to prefer a specific provider path."
        )
    endif()
endfunction()

function(hase_openpmd_find_bundled_config out_var)
    file(
        GLOB_RECURSE HASE_OPENPMD_BUNDLED_CONFIGS
        LIST_DIRECTORIES FALSE
        "${HASE_OPENPMD_BUNDLED_PREFIX}/openPMDConfig.cmake"
        "${HASE_OPENPMD_BUNDLED_PREFIX}/*/openPMDConfig.cmake"
    )
    list(SORT HASE_OPENPMD_BUNDLED_CONFIGS)
    if(HASE_OPENPMD_BUNDLED_CONFIGS)
        list(GET HASE_OPENPMD_BUNDLED_CONFIGS 0 HASE_OPENPMD_BUNDLED_CONFIG)
        get_filename_component(
            HASE_OPENPMD_BUNDLED_CONFIG_DIR
            "${HASE_OPENPMD_BUNDLED_CONFIG}"
            DIRECTORY
        )
        set(${out_var} "${HASE_OPENPMD_BUNDLED_CONFIG_DIR}" PARENT_SCOPE)
    else()
        set(${out_var} "" PARENT_SCOPE)
    endif()
endfunction()

function(hase_openpmd_find_bundled_python out_var)
    file(
        GLOB_RECURSE HASE_OPENPMD_BUNDLED_PYTHON_MARKERS
        LIST_DIRECTORIES FALSE
        "${HASE_OPENPMD_BUNDLED_PREFIX}/openpmd_api/__init__.py"
        "${HASE_OPENPMD_BUNDLED_PREFIX}/*/openpmd_api/__init__.py"
    )
    list(SORT HASE_OPENPMD_BUNDLED_PYTHON_MARKERS)
    if(HASE_OPENPMD_BUNDLED_PYTHON_MARKERS)
        list(
            GET
            HASE_OPENPMD_BUNDLED_PYTHON_MARKERS
            0
            HASE_OPENPMD_BUNDLED_PYTHON_MARKER
        )
        get_filename_component(
            HASE_OPENPMD_BUNDLED_PACKAGE_DIR
            "${HASE_OPENPMD_BUNDLED_PYTHON_MARKER}"
            DIRECTORY
        )
        get_filename_component(
            HASE_OPENPMD_BUNDLED_SITE_PACKAGES
            "${HASE_OPENPMD_BUNDLED_PACKAGE_DIR}"
            DIRECTORY
        )
        set(${out_var} "${HASE_OPENPMD_BUNDLED_SITE_PACKAGES}" PARENT_SCOPE)
    else()
        set(${out_var} "" PARENT_SCOPE)
    endif()
endfunction()

function(hase_openpmd_provider_stamp out_var)
    if(HASE_MPI_ENABLED)
        set(HASE_OPENPMD_PROVIDER_MPI ON)
    else()
        set(HASE_OPENPMD_PROVIDER_MPI OFF)
    endif()
    set(HASE_OPENPMD_PROVIDER_STAMP_CONTENT
        "openPMD=${HASE_OPENPMD_GIT_REPOSITORY}@${HASE_OPENPMD_GIT_TAG}\n"
        "ADIOS2=${HASE_ADIOS2_GIT_REPOSITORY}@${HASE_ADIOS2_GIT_TAG}\n"
        "HDF5=${HASE_HDF5_GIT_REPOSITORY}@${HASE_HDF5_GIT_TAG}\n"
        "USE_ADIOS2=${HASE_OPENPMD_USE_ADIOS2}\n"
        "USE_HDF5=${HASE_OPENPMD_USE_HDF5}\n"
        "USE_SST=${HASE_OPENPMD_USE_SST}\n"
        "FETCH_ADIOS2=${HASE_OPENPMD_FETCH_ADIOS2}\n"
        "FETCH_HDF5=${HASE_OPENPMD_FETCH_HDF5}\n"
        "SUPERBUILD=${HASE_OPENPMD_SUPERBUILD}\n"
        "PYTHON=${HASE_OPENPMD_BUILD_PYTHON_BINDINGS}\n"
        "MPI=${HASE_OPENPMD_PROVIDER_MPI}\n"
        "STAGING=separate-dependency-installs-v1\n"
        "CMAKE=${CMAKE_VERSION}\n"
        "C_COMPILER=${CMAKE_C_COMPILER}|${CMAKE_C_COMPILER_ID}|${CMAKE_C_COMPILER_VERSION}\n"
        "CXX_COMPILER=${CMAKE_CXX_COMPILER}|${CMAKE_CXX_COMPILER_ID}|${CMAKE_CXX_COMPILER_VERSION}\n"
    )
    string(
        SHA256
        HASE_OPENPMD_PROVIDER_STAMP_HASH
        "${HASE_OPENPMD_PROVIDER_STAMP_CONTENT}"
    )
    set(${out_var}
        "${HASE_OPENPMD_PROVIDER_STAMP_HASH}\n${HASE_OPENPMD_PROVIDER_STAMP_CONTENT}"
        PARENT_SCOPE
    )
endfunction()

function(hase_openpmd_run_provider_stage stage_name template_file)
    set(HASE_OPENPMD_STAGE_SOURCE_DIR
        "${CMAKE_BINARY_DIR}/hase-openpmd-provider/src/${stage_name}"
    )
    set(HASE_OPENPMD_STAGE_BUILD_DIR
        "${HASE_OPENPMD_BUNDLED_BUILD_DIR}/${stage_name}"
    )
    file(MAKE_DIRECTORY "${HASE_OPENPMD_STAGE_SOURCE_DIR}")
    configure_file(
        "${CMAKE_CURRENT_LIST_DIR}/${template_file}"
        "${HASE_OPENPMD_STAGE_SOURCE_DIR}/CMakeLists.txt"
        @ONLY
    )

    set(HASE_OPENPMD_STAGE_PREFIX_PATH "${HASE_OPENPMD_BUNDLED_PREFIX}")
    if(CMAKE_PREFIX_PATH)
        string(APPEND HASE_OPENPMD_STAGE_PREFIX_PATH ";${CMAKE_PREFIX_PATH}")
    endif()

    set(HASE_OPENPMD_STAGE_CONFIGURE_COMMAND
        "${CMAKE_COMMAND}"
        "-S"
        "${HASE_OPENPMD_STAGE_SOURCE_DIR}"
        "-B"
        "${HASE_OPENPMD_STAGE_BUILD_DIR}"
        "-DCMAKE_INSTALL_PREFIX=${HASE_OPENPMD_BUNDLED_PREFIX}"
        "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
        "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
        "-DCMAKE_PREFIX_PATH=${HASE_OPENPMD_STAGE_PREFIX_PATH}"
    )
    if(CMAKE_GENERATOR)
        list(
            APPEND
            HASE_OPENPMD_STAGE_CONFIGURE_COMMAND
            "-G"
            "${CMAKE_GENERATOR}"
        )
    endif()
    if(CMAKE_GENERATOR_PLATFORM)
        list(
            APPEND
            HASE_OPENPMD_STAGE_CONFIGURE_COMMAND
            "-A"
            "${CMAKE_GENERATOR_PLATFORM}"
        )
    endif()
    if(CMAKE_GENERATOR_TOOLSET)
        list(
            APPEND
            HASE_OPENPMD_STAGE_CONFIGURE_COMMAND
            "-T"
            "${CMAKE_GENERATOR_TOOLSET}"
        )
    endif()
    if(CMAKE_MAKE_PROGRAM)
        list(
            APPEND
            HASE_OPENPMD_STAGE_CONFIGURE_COMMAND
            "-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}"
        )
    endif()
    if(ADIOS2_DIR)
        list(
            APPEND
            HASE_OPENPMD_STAGE_CONFIGURE_COMMAND
            "-DADIOS2_DIR=${ADIOS2_DIR}"
        )
    endif()
    if(
        stage_name STREQUAL "openpmd"
        AND HASE_OPENPMD_USE_HDF5
        AND HASE_OPENPMD_FETCH_HDF5
    )
        list(
            APPEND
            HASE_OPENPMD_STAGE_CONFIGURE_COMMAND
            "-DHDF5_DIR=${HASE_OPENPMD_BUNDLED_PREFIX}/lib/cmake/hdf5"
        )
    elseif(HDF5_DIR)
        list(
            APPEND
            HASE_OPENPMD_STAGE_CONFIGURE_COMMAND
            "-DHDF5_DIR=${HDF5_DIR}"
        )
    endif()

    message(
        STATUS
        "Configuring HASE-managed bundled openPMD provider stage: ${stage_name}"
    )
    execute_process(
        COMMAND ${HASE_OPENPMD_STAGE_CONFIGURE_COMMAND}
        RESULT_VARIABLE HASE_OPENPMD_STAGE_CONFIGURE_RESULT
    )
    if(NOT HASE_OPENPMD_STAGE_CONFIGURE_RESULT EQUAL 0)
        message(
            FATAL_ERROR
            "Configuring the HASE-managed openPMD provider stage '${stage_name}' failed"
        )
    endif()

    set(HASE_OPENPMD_STAGE_BUILD_COMMAND
        "${CMAKE_COMMAND}"
        "--build"
        "${HASE_OPENPMD_STAGE_BUILD_DIR}"
        "--target"
        "install"
    )
    if(CMAKE_BUILD_TYPE)
        list(
            APPEND
            HASE_OPENPMD_STAGE_BUILD_COMMAND
            "--config"
            "${CMAKE_BUILD_TYPE}"
        )
    endif()
    message(
        STATUS
        "Building/installing HASE-managed bundled openPMD provider stage: ${stage_name}"
    )
    execute_process(
        COMMAND ${HASE_OPENPMD_STAGE_BUILD_COMMAND}
        RESULT_VARIABLE HASE_OPENPMD_STAGE_BUILD_RESULT
    )
    if(NOT HASE_OPENPMD_STAGE_BUILD_RESULT EQUAL 0)
        message(
            FATAL_ERROR
            "Building/installing the HASE-managed openPMD provider stage '${stage_name}' failed"
        )
    endif()
endfunction()

function(hase_openpmd_bootstrap_bundled_provider)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/hase-openpmd-provider")
    set(HASE_OPENPMD_BUNDLED_STAMP
        "${HASE_OPENPMD_BUNDLED_PREFIX}/hase-openpmd-provider.stamp"
    )
    hase_openpmd_provider_stamp(HASE_OPENPMD_EXPECTED_STAMP)
    hase_openpmd_find_bundled_config(HASE_OPENPMD_EXISTING_CONFIG_DIR)

    set(HASE_OPENPMD_REUSE_PROVIDER FALSE)
    if(
        NOT HASE_OPENPMD_BUNDLED_REBUILD
        AND HASE_OPENPMD_EXISTING_CONFIG_DIR
        AND EXISTS "${HASE_OPENPMD_BUNDLED_STAMP}"
    )
        file(READ "${HASE_OPENPMD_BUNDLED_STAMP}" HASE_OPENPMD_EXISTING_STAMP)
        if(
            "${HASE_OPENPMD_EXISTING_STAMP}"
                STREQUAL
                "${HASE_OPENPMD_EXPECTED_STAMP}"
        )
            set(HASE_OPENPMD_REUSE_PROVIDER TRUE)
        endif()
    endif()

    if(HASE_OPENPMD_REUSE_PROVIDER)
        message(
            STATUS
            "Reusing HASE-managed bundled openPMD provider at ${HASE_OPENPMD_BUNDLED_PREFIX}"
        )
    else()
        message(
            STATUS
            "Building HASE-managed bundled openPMD provider into ${HASE_OPENPMD_BUNDLED_PREFIX}"
        )
        file(REMOVE_RECURSE "${HASE_OPENPMD_BUNDLED_BUILD_DIR}")
        file(REMOVE_RECURSE "${HASE_OPENPMD_BUNDLED_PREFIX}")

        if(HASE_MPI_ENABLED)
            set(HASE_PROVIDER_USE_MPI ON)
        else()
            set(HASE_PROVIDER_USE_MPI OFF)
        endif()

        if(HASE_OPENPMD_USE_HDF5 AND HASE_OPENPMD_FETCH_HDF5)
            hase_openpmd_run_provider_stage(hdf5 HaseHdf5Provider.cmake.in)
        endif()
        if(HASE_OPENPMD_USE_ADIOS2 AND HASE_OPENPMD_FETCH_ADIOS2)
            hase_openpmd_run_provider_stage(adios2 HaseAdios2Provider.cmake.in)
        endif()
        hase_openpmd_run_provider_stage(openpmd HaseOpenPmdProvider.cmake.in)

        file(
            WRITE
            "${HASE_OPENPMD_BUNDLED_STAMP}"
            "${HASE_OPENPMD_EXPECTED_STAMP}"
        )
        hase_openpmd_find_bundled_config(HASE_OPENPMD_EXISTING_CONFIG_DIR)
    endif()

    if(NOT HASE_OPENPMD_EXISTING_CONFIG_DIR)
        message(
            FATAL_ERROR
            "The HASE-managed bundled openPMD provider did not install openPMDConfig.cmake under "
            "${HASE_OPENPMD_BUNDLED_PREFIX}."
        )
    endif()

    set(openPMD_DIR
        "${HASE_OPENPMD_EXISTING_CONFIG_DIR}"
        CACHE PATH
        "openPMD CMake config directory provided by the HASE-managed bundled provider"
        FORCE
    )
    # Keep the bundled provider prefix visible after this function returns so
    # openPMDConfig.cmake can resolve transitive find_dependency() calls such as
    # ADIOS2 and HDF5 during the outer HASE find_package(openPMD).
    set(CMAKE_PREFIX_PATH
        "${HASE_OPENPMD_BUNDLED_PREFIX};${CMAKE_PREFIX_PATH}"
        PARENT_SCOPE
    )

    if(
        HASE_OPENPMD_BUILD_PYTHON_BINDINGS
        AND NOT HASE_OPENPMD_PYTHON_PACKAGE_DIR
    )
        hase_openpmd_find_bundled_python(HASE_OPENPMD_BUNDLED_PYTHON_DIR)
        if(HASE_OPENPMD_BUNDLED_PYTHON_DIR)
            set(HASE_OPENPMD_PYTHON_PACKAGE_DIR
                "${HASE_OPENPMD_BUNDLED_PYTHON_DIR}"
                CACHE PATH
                "Directory containing the openpmd_api Python package matching openPMD::openPMD"
                FORCE
            )
        else()
            message(
                WARNING
                "HASE_OPENPMD_BUILD_PYTHON_BINDINGS=ON, but no installed openpmd_api package "
                "was found under ${HASE_OPENPMD_BUNDLED_PREFIX}. Runtime Python must provide "
                "a matching openpmd_api module."
            )
        endif()
    endif()
endfunction()

set(HASE_OPENPMD_ACTUAL_PROVIDER "")
if(HASE_OPENPMD_PROVIDER STREQUAL "system")
    message(STATUS "Using system openPMD-api for the HASE openPMD transport")
    find_package(openPMD CONFIG QUIET)
    if(NOT openPMD_FOUND)
        message(
            FATAL_ERROR
            "HASE_OPENPMD_PROVIDER=system requires an installed openPMD-api C++ package, "
            "but CMake could not find openPMD. Point CMake to it with openPMD_DIR or "
            "CMAKE_PREFIX_PATH, or use HASE_OPENPMD_PROVIDER=bundled."
        )
    endif()
    hase_openpmd_validate_found("system")
    set(HASE_OPENPMD_ACTUAL_PROVIDER "system")
elseif(HASE_OPENPMD_PROVIDER STREQUAL "auto")
    message(
        STATUS
        "HASE_OPENPMD_PROVIDER=auto: probing for a system openPMD-api provider"
    )
    find_package(openPMD CONFIG QUIET)
    if(openPMD_FOUND AND TARGET openPMD::openPMD)
        message(
            STATUS
            "HASE_OPENPMD_PROVIDER=auto: using system openPMD-api provider"
        )
        set(HASE_OPENPMD_ACTUAL_PROVIDER "system")
    else()
        message(
            STATUS
            "HASE_OPENPMD_PROVIDER=auto: no system provider found; using bundled provider"
        )
        hase_openpmd_bootstrap_bundled_provider()
        find_package(openPMD CONFIG REQUIRED)
        hase_openpmd_validate_found("bundled")
        set(HASE_OPENPMD_ACTUAL_PROVIDER "bundled")
    endif()
else()
    hase_openpmd_bootstrap_bundled_provider()
    find_package(openPMD CONFIG REQUIRED)
    hase_openpmd_validate_found("bundled")
    set(HASE_OPENPMD_ACTUAL_PROVIDER "bundled")
endif()

if(HASE_OPENPMD_ACTUAL_PROVIDER STREQUAL "system")
    set(HASE_USE_SYSTEM_OPENPMD TRUE)
    set(HASE_BUILD_OPENPMD_FROM_SOURCE
        OFF
        CACHE BOOL
        "Deprecated compatibility alias; use HASE_OPENPMD_PROVIDER"
        FORCE
    )
else()
    set(HASE_USE_SYSTEM_OPENPMD FALSE)
    set(HASE_BUILD_OPENPMD_FROM_SOURCE
        ON
        CACHE BOOL
        "Deprecated compatibility alias; use HASE_OPENPMD_PROVIDER"
        FORCE
    )
endif()

hase_openpmd_configure_python("${HASE_OPENPMD_ACTUAL_PROVIDER}")
