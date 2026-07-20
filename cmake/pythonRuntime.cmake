# Install the Python frontend's private C++ runtime files. The public Python
# API communicates with calcPhiASE exclusively through openPMD.
set(HASE_PYTHON_RUNTIME_DIR "${CMAKE_BINARY_DIR}/python/pyInclude/_runtime")

if(TARGET hase_openpmd_python)
    add_dependencies(calcPhiASE hase_openpmd_python)
endif()
add_dependencies(calcPhiASE HaseAlpakaBackendNames HaseOpenPmdBackendProbe)

# Wheels do not vendor openPMD-api runtime libraries or generated Python
# bindings. The frontend uses the matching provider from its runtime environment.
set(HASE_USE_SYSTEM_OPENPMD_PY True)
string(REPLACE "\\" "\\\\" HASE_OPENPMD_PYTHON_PACKAGE_DIR_ESCAPED "${HASE_OPENPMD_PYTHON_PACKAGE_DIR}")
string(REPLACE "\"" "\\\"" HASE_OPENPMD_PYTHON_PACKAGE_DIR_ESCAPED "${HASE_OPENPMD_PYTHON_PACKAGE_DIR_ESCAPED}")

file(MAKE_DIRECTORY "${HASE_PYTHON_RUNTIME_DIR}")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/pythonRuntimeConfig.py.in"
    "${HASE_PYTHON_RUNTIME_DIR}/_config.py"
    @ONLY
)

foreach(HASE_RUNTIME_TARGET IN ITEMS calcPhiASE HaseAlpakaBackendNames HaseOpenPmdBackendProbe)
    add_custom_command(
        TARGET ${HASE_RUNTIME_TARGET}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:${HASE_RUNTIME_TARGET}>"
            "${HASE_PYTHON_RUNTIME_DIR}/$<TARGET_FILE_NAME:${HASE_RUNTIME_TARGET}>"
        COMMENT "Copying ${HASE_RUNTIME_TARGET} to the private HASE Python runtime"
        VERBATIM
    )
endforeach()

install(TARGETS calcPhiASE RUNTIME DESTINATION pyInclude/_runtime)
install(TARGETS HaseAlpakaBackendNames HaseOpenPmdBackendProbe LIBRARY DESTINATION pyInclude/_runtime)
install(FILES "${HASE_PYTHON_RUNTIME_DIR}/_config.py" DESTINATION pyInclude/_runtime)
install(FILES HASEonGPU.py DESTINATION .)
