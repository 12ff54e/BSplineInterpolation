cmake_minimum_required(VERSION 3.12.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_definitions($<$<CONFIG:Debug>:_DEBUG>)

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
message(STATUS)
