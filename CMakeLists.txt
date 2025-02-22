cmake_minimum_required(VERSION 3.12.0)
project(
    BSplineInterpolation
    VERSION 2.0.0
    HOMEPAGE_URL "https://github.com/12ff54e/BSplineInterpolation.git"
    LANGUAGES CXX)
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

# Requires c++11, do not need any compiler extension.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

message(STATUS "Entering test directory ...")
add_subdirectory(test)
message(STATUS "Leaving test directory ...")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
