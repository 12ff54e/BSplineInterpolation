[![CMake](https://github.com/12ff54e/BSplineInterpolation/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/12ff54e/BSplineInterpolation/actions/workflows/cmake.yml)

# BSpline Interpolation Lib

This library use B-spline to interpolate given data of any dimension on a Cartesian mesh grid.

## Perquisite

This library requires C++11. CMake is required only for build and run tests.

## How to use

Since the library is header only, include the file include/Interpolation.hpp and one can use all its functionality under the namespace intp. Check test/interpolation-test.cpp and test/interpolation-template-test.cpp for more info.

### Packaging using CMake

This library can be installed using CMake by executing the following commands:
```bash
git clone https://github.com/12ff54e/BSplineInterpolation.git
cd BSplineInterpolation

cmake -B build -DCMAKE_INSTALL_PREFIX=/folder/you/want/to/install/in
cmake --install build
```
If you do not add `CMAKE_INSTALL_PREFIX` in configuration, you may need sudo to perform the installation. Then in the CMakeLists.txt file you can do:
```cmake
find_package(BSplineInterpolation) # can also specify version here

add_executable(main main.cpp)
target_link_libraries(main BSplineInterpolation)
```
and in the source file:
```cpp
#include <BSplineInterpolation/Interpolation.hpp>

using namespace intp;
InterpolationFunction<double, 2> func(
    3, // interpolation function order
    z_mesh, // mesh storing z values, an object of type intp::Mesh<double, 2>
    std::make_pair(x_min, x_max), // x range
    std::make_pair(y_min, y_max)); // y range
```
Note: this project follows [Semantic Version 2.0.0](https://semver.org/) so the interface will be compatible within one major version.

### Matlab Interface

Check matlab/Example.m for instructions. Please notice that matlab/bspline.mexw64 is compiled for windows platform, you may need to compile your own MEX file from bspline.cpp (basically a wrapper) on other platforms.

## Note

- Limited compiler support: only tested using g++ 11.3 and clang++ 14.0 (on Ubuntu).
- The Interpolation Template feature is not efficient in dimension 2 or higher.

## Known Issues

- No

## Future Plan

[x] Support multithreading.
[ ] Interpolation with boundary derivative specified.
