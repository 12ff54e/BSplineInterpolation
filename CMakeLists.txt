cmake_minimum_required(VERSION 3.11.0)
project(BSplineInterpolation VERSION 0.1.0)
include(FetchContent)

add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-std=c++11>)
add_compile_options(-Wall)
add_compile_options(-g)
add_definitions(-D_DEBUG)

FetchContent_Declare(
    eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen
    GIT_TAG 34e5f34b391e9b964917d54139aae1ecf522d4e5 # 3.4.0
)
FetchContent_MakeAvailable(eigen3)
set(BUILD_TESTING OFF)
include_directories(${eigen3_SOURCE_DIR})

add_executable(bsp_test test/bspline-test.cpp)
add_executable(itp_test test/interpolation-test.cpp)
add_executable(mesh_test test/mesh-test.cpp)
add_executable(util_test test/util-test.cpp)

include(CTest)
enable_testing()
add_test(intp-bspline bsp_test)
add_test(intp-interp itp_test)
add_test(intp-mesh mesh_test)
add_test(intp-util util_test)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
