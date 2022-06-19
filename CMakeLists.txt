cmake_minimum_required(VERSION 3.12.0)
project(BSplineInterpolation VERSION 0.1.0)
include(FetchContent)
include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall)
add_compile_options(-g)
add_definitions(-D_DEBUG)

find_package(Eigen3 3.0 CONFIG)

message(STATUS)
if(TARGET Eigen3::Eigen)
    message(STATUS "Eigen3 v${EIGEN3_VERSION_STRING} found in ${EIGEN3_INCLUDE_DIR}")
    else()
    message(STATUS "Download Eigen 3.4.0 from its repository ...")
    FetchContent_Declare(
        eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen
        GIT_TAG 34e5f34b391e9b964917d54139aae1ecf522d4e5 # 3.4.0
    )
    FetchContent_MakeAvailable(eigen3)
    include_directories(${eigen3_SOURCE_DIR})
endif()

list(APPEND tests util-test mesh-test bspline-test interpolation-test interpolation-speed-test interpolation-template-test)

list(LENGTH tests test_num)
message(STATUS)
message(STATUS  ${test_num} " tests found: ")
enable_testing()
foreach(test ${tests})
    string(JOIN "" test_src_file test/ ${test} .cpp)
    add_executable(${test} ${test_src_file})
    target_link_libraries(${test} PRIVATE Eigen3::Eigen)

    string(JOIN "-" test_name intp ${test})
    add_test(${test_name} ${test})
    message(STATUS ${test} " from source file " ${test_src_file})
endforeach(test ${tests})
message(STATUS)

target_compile_options(interpolation-speed-test PRIVATE -O2)
target_compile_options(interpolation-template-test PRIVATE -O2)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
