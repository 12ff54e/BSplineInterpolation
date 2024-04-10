cmake_minimum_required(VERSION 3.12.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type to debug.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
message(STATUS "Build type = " ${CMAKE_BUILD_TYPE})

# Add options according to compiler.
if(MSVC)
    set(CMAKE_CXX_FLAGS "/W4 /EHsc")
else()
    set(debug_flags " -fsanitize=address")
    if(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
        set(CMAKE_CXX_FLAGS "-Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-exit-time-destructors")
    else()
        set(CMAKE_CXX_FLAGS "-Wall -Wextra")
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL Debug)
        string(APPEND CMAKE_CXX_FLAGS ${debug_flags})
    endif()
endif()

add_compile_definitions($<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:INTP_DEBUG>
                        $<$<CONFIG:Release>:INTP_MULTITHREAD>)

add_compile_definitions(INTP_PERIODIC_NO_DUMMY_POINT
                        INTP_CELL_LAYOUT)

#add_compile_definitions(INTP_TRACE)

# Specify tests
list(APPEND tests "util-test" "mesh-test" "band-matrix-and-solver-test" "bspline-test" "interpolation-test" "interpolation-speed-test" "interpolation-template-test")

list(LENGTH tests test_num)
message(STATUS)
message(STATUS  ${test_num} " tests found: ")
enable_testing()
foreach(test ${tests})
    string(JOIN "" test_src_file src/ ${test} .cpp)
    add_executable(${test} ${test_src_file})
    target_compile_features(${test} PRIVATE cxx_std_17)
    target_include_directories(
        ${test} PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/include>)

    string(JOIN "-" test_name intp ${test})
    add_test(${test_name} ${test})
    message(STATUS ${test} " from source file " ${test_src_file})
endforeach(test ${tests})

# target_compile_options(interpolation-speed-test PRIVATE $<$<CONFIG:Release>:-march=native>)

message(STATUS)
