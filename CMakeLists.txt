cmake_minimum_required(VERSION 3.0.0)
project(mct VERSION 0.1.0)

add_compile_options(-std=c++11)

add_executable(mct src/main.cpp)
add_executable(bmt_test test/band-matrix-test.cpp)
add_executable(bsp_test test/bspline-test.cpp)
add_executable(itp_test test/interpolation-test.cpp)
add_executable(mesh_test test/mesh-test.cpp)
add_executable(util_test test/util-test.cpp)

include(CTest)
enable_testing()
add_test(band-matrix bmt_test)
add_test(bspline bsp_test)
add_test(interp itp_test)
add_test(mesh mesh_test)
add_test(util util_test)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
