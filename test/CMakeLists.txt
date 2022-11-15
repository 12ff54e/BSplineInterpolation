cmake_minimum_required(VERSION 3.12.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

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

add_compile_definitions($<$<CONFIG:Debug>:_DEBUG>)
# add_compile_definitions($<$<CONFIG:Debug>:_TRACE>)

add_compile_definitions($<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>:_MULTITHREAD>)

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
