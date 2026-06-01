# -----------------------
# Detect CUDA
# -----------------------
set(HASE_HAVE_CUDA OFF)
if(DEFINED ENV{CUDACXX} AND EXISTS "$ENV{CUDACXX}")
    set(HASE_HAVE_CUDA ON)
    set(alpaka_DEP_CUDA ON CACHE BOOL "Enable Alpaka CUDA backend" FORCE)
    message(
        STATUS
        "Found CUDA via CUDACXX=$ENV{CUDACXX} -> enabling Alpaka CUDA backend"
    )
else()
    find_program(NVCC_EXECUTABLE nvcc)
    if(NVCC_EXECUTABLE)
        set(HASE_HAVE_CUDA ON)
        set(alpaka_DEP_CUDA ON CACHE BOOL "Enable Alpaka CUDA backend" FORCE)
        message(
            STATUS
            "Found CUDA via nvcc=${NVCC_EXECUTABLE} -> enabling Alpaka CUDA backend"
        )
    endif()
endif()
if(NOT HASE_HAVE_CUDA)
    set(alpaka_DEP_CUDA OFF CACHE BOOL "Enable Alpaka CUDA backend" FORCE)
    message(STATUS "CUDA not found -> disabling Alpaka CUDA backend")
endif()

# -----------------------
# Detect HIP
# -----------------------
set(HASE_HAVE_HIP OFF)
if(DEFINED ENV{HIPCXX} AND EXISTS "$ENV{HIPCXX}")
    set(HASE_HAVE_HIP ON)
    set(alpaka_DEP_HIP ON CACHE BOOL "Enable Alpaka HIP backend" FORCE)
    message(
        STATUS
        "Found HIP via HIPCXX=$ENV{HIPCXX} -> enabling Alpaka HIP backend"
    )
else()
    find_program(HIPCC_EXECUTABLE hipcc)
    if(HIPCC_EXECUTABLE)
        set(HASE_HAVE_HIP ON)
        set(alpaka_DEP_HIP ON CACHE BOOL "Enable Alpaka HIP backend" FORCE)
        message(
            STATUS
            "Found HIP via hipcc=${HIPCC_EXECUTABLE} -> enabling Alpaka HIP backend"
        )
    elseif(DEFINED ENV{ROCM_PATH} AND IS_DIRECTORY "$ENV{ROCM_PATH}")
        set(HASE_HAVE_HIP ON)
        set(alpaka_DEP_HIP ON CACHE BOOL "Enable Alpaka HIP backend" FORCE)
        message(
            STATUS
            "Found HIP via ROCM_PATH=$ENV{ROCM_PATH} -> enabling Alpaka HIP backend"
        )
    endif()
endif()
if(NOT HASE_HAVE_HIP)
    set(alpaka_DEP_HIP OFF CACHE BOOL "Enable Alpaka HIP backend" FORCE)
    message(STATUS "HIP not found -> disabling Alpaka HIP backend")
endif()

# -----------------------
# Detect OpenMP
# -----------------------
set(HASE_HAVE_OMP OFF)
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    set(HASE_HAVE_OMP ON)
    set(alpaka_DEP_OMP ON CACHE BOOL "Enable Alpaka OpenMP backend" FORCE)
    message(STATUS "Found OpenMP -> enabling Alpaka OpenMP backend")
endif()
if(NOT HASE_HAVE_OMP)
    set(alpaka_DEP_OMP OFF CACHE BOOL "Enable Alpaka OpenMP backend" FORCE)
    message(STATUS "OpenMP not found -> disabling Alpaka OpenMP backend")
endif()

# -----------------------
# Detect TBB (oneTBB or classic)
# -----------------------
set(HASE_HAVE_TBB OFF)
find_package(TBB QUIET)
if(TBB_FOUND)
    set(HASE_HAVE_TBB ON)
    set(alpaka_DEP_TBB ON CACHE BOOL "Enable Alpaka TBB backend" FORCE)
    message(STATUS "Found TBB -> enabling Alpaka TBB backend")
else()
    find_package(oneTBB QUIET)
    if(oneTBB_FOUND)
        set(HASE_HAVE_TBB ON)
        set(alpaka_DEP_TBB ON CACHE BOOL "Enable Alpaka TBB backend" FORCE)
        message(STATUS "Found oneTBB -> enabling Alpaka TBB backend")
    endif()
endif()
if(NOT HASE_HAVE_TBB)
    set(alpaka_DEP_TBB OFF CACHE BOOL "Enable Alpaka TBB backend" FORCE)
    message(STATUS "TBB not found -> disabling Alpaka TBB backend")
endif()
# -----------------------
# Backend conflict check
# -----------------------
if(HASE_HAVE_HIP AND HASE_HAVE_CUDA)
    message(
        FATAL_ERROR
        "Currently HIP and CUDA backends are mutually exclusive.\n"
        "Either unload one of the dependencies or target a backend more specifically using:\n"
        "  HASE_SELECT_BACKEND_ALPAKA=ON\n"
        "which allows you to select backends using Alpaka's CMake knobs."
    )
endif()
