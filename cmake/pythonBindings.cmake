set_target_properties(hase PROPERTIES POSITION_INDEPENDENT_CODE ON)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG QUIET)

if(NOT pybind11_FOUND AND NOT TARGET pybind11::module)
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
if(TARGET hase_openpmd_python)
    add_dependencies(HASEonGPU hase_openpmd_python)
    add_dependencies(calcPhiASE hase_openpmd_python)
endif()
add_dependencies(HASEonGPU HaseOpenPmdBackendProbe)
if(HASE_FORWARD_LOGGING)
    set(HASE_FORWARD_LOGGING_PY True)
else()
    set(HASE_FORWARD_LOGGING_PY False)
endif()
# HASE wheels do not vendor openPMD-api runtime libraries or Python bindings.
# Even when CMake fetches openPMD-api for the build, the installed package uses
# the runtime environment's openpmd_api module and dynamic libraries.
set(HASE_USE_SYSTEM_OPENPMD_PY True)
string(
    REPLACE
    "\\"
    "\\\\"
    HASE_OPENPMD_PYTHON_PACKAGE_DIR_ESCAPED
    "${HASE_OPENPMD_PYTHON_PACKAGE_DIR}"
)
string(
    REPLACE
    "\""
    "\\\""
    HASE_OPENPMD_PYTHON_PACKAGE_DIR_ESCAPED
    "${HASE_OPENPMD_PYTHON_PACKAGE_DIR_ESCAPED}"
)
file(MAKE_DIRECTORY "${HASE_PYTHON_RUNTIME_DIR}")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/HASEonGPU_Bindings/__init__.py"
    "${HASE_PYTHON_RUNTIME_DIR}/__init__.py"
    COPYONLY
)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/HASEonGPU_Bindings/_config.py.in"
    "${HASE_PYTHON_RUNTIME_DIR}/_config.py"
    @ONLY
)
set_target_properties(
    HASEonGPU
    PROPERTIES
        BUILD_RPATH_USE_ORIGIN ON
        BUILD_RPATH "${HASE_BUILD_TREE_RPATH}"
        CUDA_SEPARABLE_COMPILATION ON
        INSTALL_RPATH "${HASE_INSTALL_RPATH}"
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
add_custom_command(
    TARGET HaseOpenPmdBackendProbe
    POST_BUILD
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:HaseOpenPmdBackendProbe>"
        "${HASE_PYTHON_RUNTIME_DIR}/$<TARGET_FILE_NAME:HaseOpenPmdBackendProbe>"
    COMMENT
        "Copying openPMD backend-probe library to build package ${HASE_PYTHON_RUNTIME_DIR}"
    VERBATIM
)
install(TARGETS HASEonGPU LIBRARY DESTINATION HASEonGPU_Bindings)
install(TARGETS calcPhiASE RUNTIME DESTINATION HASEonGPU_Bindings)
install(TARGETS HaseAlpakaBackendNames LIBRARY DESTINATION HASEonGPU_Bindings)
install(TARGETS HaseOpenPmdBackendProbe LIBRARY DESTINATION HASEonGPU_Bindings)
install(
    FILES "${HASE_PYTHON_RUNTIME_DIR}/_config.py"
    DESTINATION HASEonGPU_Bindings
)
install(FILES HASEonGPU.py DESTINATION .)
install(
    DIRECTORY HASEonGPU_Bindings
    DESTINATION .
    FILES_MATCHING
    PATTERN "__init__.py"
)
