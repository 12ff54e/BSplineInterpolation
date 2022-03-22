# BSpline Interpolation Lib

This library use B-spline to interpolate given data of any dimension on a Cartesian mesh grid.

## Perquisite

This library is standard C++11 and relies on [Eigen3](https://eigen.tuxfamily.org/). CMake is required only for build and run tests.
- The eigen3 library should be in system path such that `#include <Eigen/SparseLU>` works properly.

## How to use

Since the library is header only, include BSplineInterpolation and one can use all its functionality under namespace intp. Check test/interpolation-test.cpp for more info.

## Disclaimer

- Limited compiler support: only tested using g++ 9.3 and clang++ 10.0 (on WSL2 Ubuntu).
- The algorithm is un-optimized and has low efficiency when constructing interpolation function on mesh of dimension more than 2.