cmake_minimum_required(VERSION 3.24)

include(FetchContent)

# Catch2 v3
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.6.0
)

FetchContent_MakeAvailable(Catch2)
include(CTest)
include(Catch)
file(GLOB test_folders CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*")

foreach(path IN LISTS test_folders)
    if(NOT IS_DIRECTORY "${path}")
        list(REMOVE_ITEM test_folders "${path}")
    endif()
endforeach()

list(APPEND test_folders "${CMAKE_CURRENT_SOURCE_DIR}")

foreach(test_folder ${test_folders})
    message(DEBUG "test folder: ${test_folder}")

    file(GLOB test_sources "${test_folder}/**")

    # iterate over all sub folders
    foreach(source_file ${test_sources})
        # the include folder is a special case, therefore ignore it
        if(source_file STREQUAL "${test_folder}/include")
            continue()
        endif()

        get_filename_component(source_file_extension ${source_file} EXT)
        get_filename_component(source_file_name ${source_file} NAME_WE)
        # get absolute path of the folder of the source file
        get_filename_component(test_folder_name ${source_file} DIRECTORY)
        # get name of the folder
        get_filename_component(test_folder_name ${test_folder_name} NAME)

        if("${source_file_extension}" STREQUAL ".cpp")
            set(_test_target_name
                "${target_prefix}${test_folder_name}_${source_file_name}"
            )
            message(
                DEBUG
                "  add test target ${_test_target_name}: ${source_file}"
            )

            add_executable(${_test_target_name} ${source_file})
            get_filename_component(src_dir ${source_file} DIRECTORY)
            target_include_directories(${_test_target_name} PRIVATE ${src_dir})
            if(EXISTS "${test_folder}/include")
                message(STATUS "  add include folder: ${test_folder}/include")
                target_include_directories(
                    ${_test_target_name}
                    PRIVATE ${test_folder}/include
                )
            endif()
            set_target_properties(
                ${_test_target_name}
                PROPERTIES CUDA_SEPARABLE_COMPILATION ON
            )
            target_link_libraries(
                ${_test_target_name}
                PUBLIC hase::hase Catch2::Catch2WithMain
            )
            target_compile_definitions(
                ${_test_target_name}
                PRIVATE HASE_WORKDIR="${PROJECT_SOURCE_DIR}"
            )
            alpaka_finalize(${_test_target_name})
            catch_discover_tests(
                    ${_test_target_name}
                    PROPERTIES
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            )
            add_custom_command(
                TARGET ${_test_target_name}
                POST_BUILD
                COMMAND
                    ${CMAKE_COMMAND} -E copy_directory
                    "${CMAKE_CURRENT_SOURCE_DIR}/data"
                    "${CMAKE_BINARY_DIR}/data"
            )
        endif()
    endforeach()
endforeach()
# ------------------------------------------------------------
# Python test dependencies
# ------------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter REQUIRED)

execute_process(
    COMMAND "${Python3_EXECUTABLE}" -m pytest --version
    RESULT_VARIABLE PYTEST_CHECK_RESULT
    OUTPUT_QUIET
    ERROR_QUIET
)

if(NOT PYTEST_CHECK_RESULT EQUAL 0)
    message(
        WARNING
        "pytest is not installed for Python interpreter:\n"
        "  ${Python3_EXECUTABLE}\n\n"
        "Python tests will not be available.\n"
        "Install pytest manually with:\n"
        "  ${Python3_EXECUTABLE} -m pip install --user pytest"
    )
    set(HASE_PYTEST_AVAILABLE OFF)
else()
    set(HASE_PYTEST_AVAILABLE ON)
endif()
if(HASE_PYTEST_AVAILABLE)
    add_test(
        NAME python_tests
        COMMAND "${Python3_EXECUTABLE}" -m pytest "${CMAKE_CURRENT_SOURCE_DIR}"
    )
endif()
