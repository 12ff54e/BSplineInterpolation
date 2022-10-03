[![CMake](https://github.com/12ff54e/BSplineInterpolation/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/12ff54e/BSplineInterpolation/actions/workflows/cmake.yml)

# BSpline Interpolation Lib

This library use B-spline to interpolate given data of any dimension on a Cartesian mesh grid.

## Perquisite

This library is standard C++11. CMake is required only for build and run tests.

## How to use

Since the library is header only, include BSplineInterpolation and one can use all its functionality under namespace intp. Check test/interpolation-test.cpp and test/interpolation-template-test.cpp for more info.

PS: main branch is under-development and is unstable, checkout the latest tag.

## Note

- Limited compiler support: only tested using g++ 10.0 and clang++ 11.0 (on WSL2 Ubuntu).
- The Interpolation Template feature is not efficient in dimension 2 or higher.

## Known Issues

- No