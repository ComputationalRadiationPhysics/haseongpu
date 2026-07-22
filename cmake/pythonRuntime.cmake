# Generate the runtime metadata consumed by the thin Python frontend. Native
# artifacts stay in the durable CMake build selected by HASE_RUNTIME_DIR.
set(HASE_PYTHON_RUNTIME_DIR "${CMAKE_BINARY_DIR}/python/pyInclude/_runtime")

if(TARGET hase_openpmd_python AND TARGET calcPhiASE)
    add_dependencies(calcPhiASE hase_openpmd_python)
endif()
if(TARGET calcPhiASE)
    add_dependencies(calcPhiASE HaseAlpakaBackendNames HaseOpenPmdBackendProbe)
endif()

# Wheels do not vendor openPMD-api runtime libraries or generated Python
# bindings. The frontend uses the matching provider from its runtime environment.
set(HASE_USE_SYSTEM_OPENPMD_PY True)
string(REPLACE "\\" "\\\\" HASE_RUNTIME_DIR_ESCAPED "${HASE_RUNTIME_DIR}")
string(
    REPLACE
    "\""
    "\\\""
    HASE_RUNTIME_DIR_ESCAPED
    "${HASE_RUNTIME_DIR_ESCAPED}"
)
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
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/pythonRuntimeConfig.py.in"
    "${HASE_PYTHON_RUNTIME_DIR}/_config.py"
    @ONLY
)

foreach(
    HASE_RUNTIME_TARGET
    IN
    ITEMS calcPhiASE HaseAlpakaBackendNames HaseOpenPmdBackendProbe
)
    if(TARGET ${HASE_RUNTIME_TARGET})
        add_custom_command(
            TARGET ${HASE_RUNTIME_TARGET}
            POST_BUILD
            COMMAND
                ${CMAKE_COMMAND} -E copy_if_different
                "$<TARGET_FILE:${HASE_RUNTIME_TARGET}>"
                "${HASE_PYTHON_RUNTIME_DIR}/$<TARGET_FILE_NAME:${HASE_RUNTIME_TARGET}>"
            COMMENT
                "Copying ${HASE_RUNTIME_TARGET} to the HASE runtime metadata directory"
            VERBATIM
        )
    endif()
endforeach()

if(TARGET calcPhiASE)
    install(
        TARGETS calcPhiASE
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT runtime
    )
endif()
if(TARGET HaseAlpakaBackendNames AND TARGET HaseOpenPmdBackendProbe)
    install(
        TARGETS HaseAlpakaBackendNames HaseOpenPmdBackendProbe
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT runtime
    )
endif()
install(
    FILES "${HASE_PYTHON_RUNTIME_DIR}/_config.py"
    DESTINATION pyInclude/_runtime
    COMPONENT python
)
install(FILES HASEonGPU.py DESTINATION . COMPONENT python)
