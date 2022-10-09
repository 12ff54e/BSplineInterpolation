[![CMake](https://github.com/12ff54e/BSplineInterpolation/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/12ff54e/BSplineInterpolation/actions/workflows/cmake.yml)

# BSpline Interpolation Lib

This library use B-spline to interpolate given data of any dimension on a Cartesian mesh grid.

## Perquisite

This library is standard C++11. CMake is required only for build and run tests.

## How to use

Since the library is header only, include the file include/Interpolation.hpp and one can use all its functionality under the namespace intp. Check test/interpolation-test.cpp and test/interpolation-template-test.cpp for more info.

### Packaging using CMake

This library can be installed using CMake by executing the following commands:
```bash
git clone https://github.com/12ff54e/BSplineInterpolation.git
cd BSplineInterpolation
git checkout latest

cmake -B build -DCMAKE_INSTALL_PREFIX=/folder/you/want/to/install/in
cmake --install build
```
If you do not add `CMAKE_INSTALL_PREFIX` in configuration, you may need sudo to perform the installation. Then in the CMakeLists.txt file you can do:
```cmake
find_package(BSplineInterpolation) # can also specify version here

add_executable(main main.cpp)
target_link_libraries(main BSplineInterpolation)
```
Note: this project follows [Semantic Version 2.0.0](https://semver.org/) so the interface will be compatible within one major version.

## Note

- **Main branch is under-development and is unstable, checkout the latest tag.**
- Limited compiler support: only tested using g++ 10.0 and clang++ 11.0 (on WSL2 Ubuntu).
- The Interpolation Template feature is not efficient in dimension 2 or higher.

## Known Issues

- No

## Future Plan

- Support multithreading.
- Interpolation with boundary derivative specified.
