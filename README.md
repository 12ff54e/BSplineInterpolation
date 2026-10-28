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

// template arguments are type of interpolated value, dimension,
//  interpolation order and type of coordinate, respectively
InterpolationFunction<double, 2, 3, double> func(
    // mesh storing z values, an object of type intp::Mesh<double, 2>
    z_mesh,
    // x range
    std::make_pair(x_min, x_max),
    // y range
    std::make_pair(y_min, y_max));
```
Note: this project follows [Semantic Version 2.0.0](https://semver.org/) so the interface will be compatible within one major version.

### Optional polynomial coefficients

`BSpline` and `InterpolationFunction` evaluate from control points by default.
Call `generate_polynomial_form()` after construction to compute and store
polynomial coefficients for each cell:

```cpp
func.generate_polynomial_form();
```

Values and derivatives then use these coefficients. This can speed up evaluation
at the cost of setup time and memory. It supports both control-point layouts,
including `INTP_CELL_LAYOUT`. `spline.has_polynomial()` reports whether a
`BSpline` stores polynomial coefficients; for an interpolation function, use
`func.spline().has_polynomial()`.

Loading new knots or control points clears the polynomial coefficients. Call
`generate_polynomial_form()` again to rebuild them. Scaling control points also
scales any stored polynomial coefficients.

After generating the polynomial form, release control-point storage if only
polynomial evaluation is needed:

```cpp
func.generate_polynomial_form();
func.release_control_points();  // Also available directly on BSpline.
```

This frees the control-point allocation, rather than only setting its size to
zero. Values, derivatives and scaling still work. Repeating either call keeps
the polynomial coefficients. Releasing nonempty controls without polynomial
coefficients throws `std::logic_error`; releasing an empty spline does nothing.
Control-point evaluation proxies cannot evaluate a spline after release. To
change its knots, first load new control points with `load_ctrlPts()`; a rejected
knot change leaves the spline intact.

### Modify Library's Behavior

You can control the bahavior of this library by defining the following macros:

- `INTP_TRACE`, `INTP_DEBUG`: Defines levels of logging, with decreasing verbosity.
- `INTP_MULTITHREAD`: Parallel control-point solving, cell-layout construction,
  and polynomial coefficient generation.
- `INTP_CUSTOM_PARALLEL_BACKEND`: The consumer supplies `intp::parallel_submit`
  and `intp::parallel_wait` instead of the default worker-pool wrappers.
- `INTP_ENABLE_ASSERTION`: Add assertion to ensure data to be interplolated is valid.
- `INTP_PERIODIC_NO_DUMMY_POINT`: Whether to accept dummy point when interpolating periodic data. If this macro is defined, number of data point is one less than the specified dimension, since the last point is implicitly set as the same of the first one.
- `INTP_STACK_ALLOCATOR`: Use stack allocator in some small dynamic memory allocation scenario.
- `INTP_CELL_LAYOUT`: Store weights in cell layout redundantly to speed up evaluation.

By default all of the above macros is not defined.

### Custom parallel backend

Define both `INTP_MULTITHREAD` and `INTP_CUSTOM_PARALLEL_BACKEND`, then define
these hooks in namespace `intp` before including the spline headers:

```cpp
template <typename Func>
std::future<void> parallel_submit(Func&& task);

void parallel_wait(std::future<void>& task);
```

For example, an adapter header can provide the definitions:

```cpp
#define INTP_MULTITHREAD
#define INTP_CUSTOM_PARALLEL_BACKEND
#include "my_parallel_backend.hpp"
#include <bsintp/Interpolation.hpp>
```

Use the same configuration and hook definitions in every translation unit.
Mark a header-defined `parallel_wait` as `inline`. An adapter may include
`bsintp/DedicatedThreadPool.hpp` before defining the hooks; the macro suppresses
only the default wrappers. The pool class remains available, and including it
does not create workers.

`parallel_submit` must return a valid `std::future<void>`, or throw without
leaving untracked work. `parallel_wait` must finish that task before returning
or propagating its error. It may consume the future. For nested work on a
shared pool, the wait hook must let the scheduler make progress, for example
by executing queued tasks. Waiting for the whole pool to become idle can
include the waiting parent itself and deadlock.

The library waits for all accepted tasks before releasing their buffers,
including when submission or work on the calling thread throws. Polynomial
construction uses at most 256 chunks, independent of the backend's worker
count. This is a task cap, not a concurrency limit.

Without `INTP_CUSTOM_PARALLEL_BACKEND`, the hooks use
`DedicatedThreadPool<void>::get_instance(8)` and `future.get()`. This preserves
the default pool's blocking waits; it does not make that backend safe for
nested waits on its own workers. Defining only `INTP_CUSTOM_PARALLEL_BACKEND`
does not enable parallel execution.

### Matlab Interface

Check matlab/Example.m for instructions. Please notice that matlab/bspline.mexw64 is compiled for windows platform, you may need to compile your own MEX file from bspline.cpp (basically a wrapper) on other platforms.

## Note

- The Interpolation Template feature is not efficient in dimension 2 or higher.

## Known Issues

- Periodic interpolation of high orders (>=6, maybe) are extremely slow.

## Future Plan

- [ ] Interpolation with boundary derivative specified.
