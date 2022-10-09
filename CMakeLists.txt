cmake_minimum_required(VERSION 3.12.0)
project(
    BSplineInterpolation
    VERSION 1.2.0
    HOMEPAGE_URL "https://github.com/12ff54e/BSplineInterpolation.git")
include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# aux vars
set(ConfigPostfix Config)
set(VerPostfix ${ConfigPostfix}Version)

# Add target for installation, use interface to indicate it's header-only.
add_library(BSplineInterpolation INTERFACE)
target_include_directories(
    BSplineInterpolation INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Version management boilerplate
write_basic_package_version_file(
    ${PROJECT_NAME}${VerPostfix}.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}${ConfigPostfix}.cmake.in
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}${ConfigPostfix}.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

# Specify install target and location
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets)
install(
    DIRECTORY src/include/
    TYPE DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
install(
    EXPORT ${PROJECT_NAME}Targets
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
    FILE ${PROJECT_NAME}${ConfigPostfix}.cmake)
install(
    FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}${VerPostfix}.cmake
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

# Set default build type to release.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Requires c++17 and no compiler extension.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add options according to compiler.
if(MSVC)
    set(CMAKE_CXX_FLAGS "/W4 /EHsc")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    set(CMAKE_CXX_FLAGS "-Weverything -Wno-c++98-compat -Wno-padded")
else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra")
endif()

add_compile_definitions($<$<CONFIG:Debug>:_DEBUG>)
# add_compile_definitions(_TRACE)

# Specify tests
list(APPEND tests util-test mesh-test band-matrix-and-solver-test bspline-test interpolation-test interpolation-speed-test interpolation-template-test)

list(LENGTH tests test_num)
message(STATUS)
message(STATUS  ${test_num} " tests found: ")
enable_testing()
foreach(test ${tests})
    string(JOIN "" test_src_file test/ ${test} .cpp)
    add_executable(${test} ${test_src_file})

    string(JOIN "-" test_name intp ${test})
    add_test(${test_name} ${test})
    message(STATUS ${test} " from source file " ${test_src_file})
endforeach(test ${tests})
message(STATUS)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
