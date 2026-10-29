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
- `INTP_CUSTOM_PARALLEL_BACKEND`: The consumer supplies `intp::parallel_for`
  instead of the default worker-pool wrapper.
- `INTP_CUSTOM_ALLOCATOR`: The consumer supplies `intp::BufferAllocator<T>`
  for control-point and polynomial coefficient buffers.
- `INTP_ENABLE_ASSERTION`: Add assertion to ensure data to be interplolated is valid.
- `INTP_PERIODIC_NO_DUMMY_POINT`: Whether to accept dummy point when interpolating periodic data. If this macro is defined, number of data point is one less than the specified dimension, since the last point is implicitly set as the same of the first one.
- `INTP_STACK_ALLOCATOR`: Use stack allocator in some small dynamic memory allocation scenario.
- `INTP_CELL_LAYOUT`: Store weights in cell layout redundantly to speed up evaluation.

By default all of the above macros is not defined.

### Custom parallel backend

Define both `INTP_MULTITHREAD` and `INTP_CUSTOM_PARALLEL_BACKEND`, then define
this hook in namespace `intp` before including the spline headers:

```cpp
template <typename Func>
void parallel_for(std::size_t count, std::size_t grain, const Func& body);
```

For example, an adapter header can provide the definition:

```cpp
#define INTP_MULTITHREAD
#define INTP_CUSTOM_PARALLEL_BACKEND
#include "my_parallel_backend.hpp"
#include <bsintp/Interpolation.hpp>
```

Use the same configuration and hook definition in every translation unit.
An adapter may include `bsintp/DedicatedThreadPool.hpp` before defining the
hook. The pool class remains available; including it does not create workers.

The hook calls `body(begin, end)` on disjoint, nonempty subranges of
`[0, count)`. A successful call covers every index exactly once. An empty
range invokes no callbacks. `grain` must be positive; it is a hint for the
number of indices per chunk. A backend may choose other chunk sizes or run
the whole range inline. The body must support concurrent calls.

The hook is synchronous: all accepted work must finish before it returns or
throws, including after submission or callback failure. This keeps borrowed
buffers and the body alive until the last callback finishes. For nested work
on a shared pool, waiting threads must let the scheduler make progress, for
example by executing queued tasks. Waiting for the whole pool to become idle
can include the waiting parent itself and deadlock.

The spline supplies its existing grain hints: square-root block counts for
control-point solving and cell-layout construction, and at most 256 chunks
with roughly 4096 or more coefficients per chunk for polynomial construction.
Custom backends control the actual partitioning and execution.

Without `INTP_CUSTOM_PARALLEL_BACKEND`, `parallel_for` uses those hints,
`DedicatedThreadPool<void>::get_instance(8)`, and blocking `future.get()`
waits. It computes the first chunk on the caller and submits the rest. This
default backend remains unsuitable for nested waits on its own workers.
Defining only `INTP_CUSTOM_PARALLEL_BACKEND` does not enable parallel execution.

### Custom buffer allocator

Define `INTP_CUSTOM_ALLOCATOR` and the allocator alias before including the
spline headers, using the same definition in every translation unit:

```cpp
#define INTP_CUSTOM_ALLOCATOR
#include "my_allocator.hpp"
namespace intp {
template <typename T>
using BufferAllocator = MyAllocator<T>;
}
#include <bsintp/Interpolation.hpp>
```

Without this macro, `BufferAllocator<T>` is
`AlignedAllocator<T, Alignment::AVX>`, preserving the default container types.
The hook does not depend on `INTP_MULTITHREAD` or change the template parameters
of `BSpline`, `InterpolationFunction`, or `InterpolationFunctionTemplate`.

The allocator supplies the two control-point solve buffers, cell-layout
storage, and polynomial coefficient storage. It does not replace allocators
for input meshes, knots, solver matrices, or scheduler tasks. Small fixed-size
scratch arrays remain on the stack. The existing initialization behavior and
`util::default_init_allocator` adaptor remain unchanged.

Supply a default-constructible standard allocator with copy, rebind, equality,
and propagation behavior appropriate to its ownership. Preserve the default
32-byte alignment (or a greater alignment required by the value type).
Allocation may occur concurrently in separate spline constructions.
Deallocation must not throw and must use the resource captured by the
allocator, even on another thread or after the current resource changes.

These are not all temporary allocations: the solve result becomes the retained
control points in row layout; cell controls and polynomial coefficients also
outlive construction. A stateful allocator must keep its resource valid until
the last owning container releases it. A constructor-local arena cannot be
reset while any of these buffers remain live. Both solve buffers use the hook
because either one can become the result, depending on dimension.

### Matlab Interface

Check matlab/Example.m for instructions. Please notice that matlab/bspline.mexw64 is compiled for windows platform, you may need to compile your own MEX file from bspline.cpp (basically a wrapper) on other platforms.

## Note

- The Interpolation Template feature is not efficient in dimension 2 or higher.

## Known Issues

- Periodic interpolation of high orders (>=6, maybe) are extremely slow.

## Future Plan

- [ ] Interpolation with boundary derivative specified.
