cmake_minimum_required(VERSION 3.12.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type to debug.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
message(STATUS "Build type = " ${CMAKE_BUILD_TYPE})

# Check compile flags availability
include(CheckCXXCompilerFlag)
# list(APPEND custom_opt_flags "-march=native" "-ffast-math")

set(asan_flag -fsanitize=address)
set(CMAKE_REQUIRED_LIBRARIES ${asan_flag}) # This is also a linker flag
check_cxx_compiler_flag(${asan_flag} compiler_support_asan)
if(${compiler_support_asan})
    add_compile_options($<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:${asan_flag}>)
    add_link_options($<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:${asan_flag}>)
endif()

foreach(custom_flag ${custom_opt_flags})
    string(REPLACE "-" "" flag_without_hyphen ${custom_flag})
    check_cxx_compiler_flag(${custom_flag} compiler_support_${flag_without_hyphen})
    if(${compiler_support_${flag_without_hyphen}})
        add_compile_options($<$<CONFIG:Release>:${custom_flag}>)
    endif()
endforeach(custom_flag ${custom_opt_flags})

# Add options according to compiler.
if(MSVC)
    set(CMAKE_CXX_FLAGS_INIT "/W4 /EHsc")
else()
    if(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
        set(CMAKE_CXX_FLAGS_INIT "-Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-exit-time-destructors")
    else()
        set(CMAKE_CXX_FLAGS_INIT "-Wall -Wextra")
    endif()
endif()

add_compile_definitions($<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:INTP_DEBUG>
                        $<$<CONFIG:Release>:INTP_MULTITHREAD>)

add_compile_definitions(INTP_PERIODIC_NO_DUMMY_POINT
                        INTP_CELL_LAYOUT)

#add_compile_definitions(INTP_TRACE)

# Specify tests
list(APPEND regular_tests "util-test" "mesh-test" "band-matrix-and-solver-test" "bspline-test" "interpolation-test")
list(APPEND speed_tests "interpolation-speed-test" "interpolation-template-test" "interpolation-eval-proxy-test")

list(LENGTH regular_tests test_num)
message(STATUS)
message(STATUS  ${test_num} " tests found: ")
enable_testing()
foreach(test ${regular_tests})
    string(JOIN "" test_src_file src/ ${test} .cpp)
    add_executable(${test} ${test_src_file})
    target_include_directories(
        ${test} PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/include>)

    string(JOIN "-" test_name intp ${test})
    add_test(${test_name} ${test})
    set_tests_properties(${test_name} PROPERTIES LABELS intp-regular)
    message(STATUS ${test} " from source file " ${test_src_file})
endforeach(test ${regular_tests})

if(CMAKE_BUILD_TYPE STREQUAL Release OR CMAKE_BUILD_TYPE STREQUAL RelWithDebInfo)
    message(STATUS)
    list(LENGTH speed_tests test_num)
    message(STATUS  ${test_num} " speed tests found: ")
    foreach(test ${speed_tests})
        string(JOIN "" test_src_file src/ ${test} .cpp)
        add_executable(${test} ${test_src_file} src/Timer.cpp)
        target_include_directories(
            ${test} PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/include>)

        string(JOIN "-" test_name intp ${test})
        add_test(${test_name} ${test})
        set_tests_properties(${test_name} PROPERTIES LABELS intp-speed)
        message(STATUS ${test} " from source file " ${test_src_file})
    endforeach(test ${speed_tests})
endif()


# target_compile_options(interpolation-speed-test PRIVATE $<$<CONFIG:Release>:-march=native>)

message(STATUS)
