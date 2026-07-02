include_guard(GLOBAL)

set(HASE_ALPAKA_GIT_REPOSITORY
    "https://github.com/alpaka-group/alpaka3.git"
    CACHE STRING
    "Git repository used when fetching alpaka"
)
set(HASE_ALPAKA_GIT_TAG
    "58065f933b3d46406c4abb7c74111567938a0ee6"
    CACHE STRING
    "Git tag or commit used when fetching alpaka"
)

if(HASE_USE_SYSTEM_ALPAKA)
    message(
        STATUS
        "HASE_USE_SYSTEM_ALPAKA=ON -> using find_package(alpaka CONFIG REQUIRED)"
    )
    find_package(alpaka CONFIG REQUIRED)
else()
    message(
        STATUS
        "HASE_USE_SYSTEM_ALPAKA=OFF -> fetching pinned alpaka into the CMake build tree"
    )
    include(FetchContent)
    FetchContent_Declare(
        alpaka
        GIT_REPOSITORY "${HASE_ALPAKA_GIT_REPOSITORY}"
        GIT_TAG "${HASE_ALPAKA_GIT_TAG}"
        ${HASE_FETCHCONTENT_EXCLUDE_FROM_ALL}
    )
    hase_fetchcontent_make_available(alpaka)
endif()
