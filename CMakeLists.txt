cmake_minimum_required(VERSION 3.0.0)
project(BSplineInterpolation VERSION 0.1.0)

add_compile_options(-std=c++11)
add_compile_options(-Wall)
add_compile_options(-g)
add_definitions(-D_DEBUG)
include_directories(/usr/include/eigen3)

add_executable(bsp_test test/bspline-test.cpp)
add_executable(itp_test test/interpolation-test.cpp)
add_executable(mesh_test test/mesh-test.cpp)
add_executable(util_test test/util-test.cpp)

include(CTest)
enable_testing()
add_test(bspline bsp_test)
add_test(interp itp_test)
add_test(mesh mesh_test)
add_test(util util_test)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
