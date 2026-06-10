add_library(
    HaseAlpakaBackendNames
    SHARED
    "${CMAKE_CURRENT_SOURCE_DIR}/src/alpakaUtils/backendNamesLibrary.cpp"
)
set_target_properties(
    HaseAlpakaBackendNames
    PROPERTIES
        BUILD_RPATH_USE_ORIGIN ON
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "${HASE_CUDA_ARCHITECTURES}"
        INSTALL_RPATH "$ORIGIN"
        INSTALL_RPATH_USE_LINK_PATH ON
        POSITION_INDEPENDENT_CODE ON
)
target_link_libraries(HaseAlpakaBackendNames PRIVATE hase::core)
alpaka_finalize(HaseAlpakaBackendNames)
