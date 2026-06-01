set_target_properties(hase PROPERTIES POSITION_INDEPENDENT_CODE ON)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG QUIET)

if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found -> fetching pybind11")

    include(FetchContent)

    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v3.0.4
    )

    FetchContent_MakeAvailable(pybind11)
endif()
pybind11_add_module(HASEonGPU HASEonGPU_Bindings/module.cpp)
set(HASE_PYTHON_RUNTIME_DIR "${CMAKE_BINARY_DIR}/python/HASEonGPU_Bindings")
file(MAKE_DIRECTORY "${HASE_PYTHON_RUNTIME_DIR}")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/HASEonGPU_Bindings/__init__.py"
    "${HASE_PYTHON_RUNTIME_DIR}/__init__.py"
    COPYONLY
)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
set_target_properties(
    HASEonGPU
    PROPERTIES
        BUILD_RPATH_USE_ORIGIN ON
        CUDA_SEPARABLE_COMPILATION ON
        INSTALL_RPATH_USE_LINK_PATH ON
)
target_link_libraries(HASEonGPU PRIVATE hase)
alpaka_finalize(HASEonGPU)
add_custom_command(
    TARGET HASEonGPU
    POST_BUILD
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:HASEonGPU>"
        "${HASE_PYTHON_RUNTIME_DIR}/$<TARGET_FILE_NAME:HASEonGPU>"
    COMMENT
        "Copying Python extension to build package ${HASE_PYTHON_RUNTIME_DIR}"
    VERBATIM
)
add_custom_command(
    TARGET HaseAlpakaBackendNames
    POST_BUILD
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:HaseAlpakaBackendNames>"
        "${HASE_PYTHON_RUNTIME_DIR}/$<TARGET_FILE_NAME:HaseAlpakaBackendNames>"
    COMMENT
        "Copying Alpaka backend-name library to build package ${HASE_PYTHON_RUNTIME_DIR}"
    VERBATIM
)
if(HASE_BUILD_PhiAse)
    add_custom_command(
        TARGET calcPhiASE
        POST_BUILD
        COMMAND
            ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:calcPhiASE>"
            "${HASE_PYTHON_RUNTIME_DIR}/$<TARGET_FILE_NAME:calcPhiASE>"
        COMMENT
            "Copying calcPhiASE executable to build package ${HASE_PYTHON_RUNTIME_DIR}"
        VERBATIM
    )
endif()
install(TARGETS HASEonGPU LIBRARY DESTINATION HASEonGPU_Bindings)
if(HASE_BUILD_PhiAse)
    install(TARGETS calcPhiASE RUNTIME DESTINATION HASEonGPU_Bindings)
endif()
install(TARGETS HaseAlpakaBackendNames LIBRARY DESTINATION HASEonGPU_Bindings)
install(FILES HASEonGPU.py DESTINATION .)
install(
    DIRECTORY HASEonGPU_Bindings
    DESTINATION .
    FILES_MATCHING
    PATTERN "__init__.py"
)
