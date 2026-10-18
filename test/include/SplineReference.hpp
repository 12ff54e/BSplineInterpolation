#ifndef INTP_TEST_SPLINE_REFERENCE
#define INTP_TEST_SPLINE_REFERENCE

#include "Assertion.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <string_view>

namespace intp_test {

struct ReferenceValue {
    long double value{};
    long double scale{};
    long double polynomial_scale{};
};

template <typename T>
void check_close(Assertion& assertion,
                 T actual,
                 ReferenceValue expected,
                 std::string_view name,
                 std::size_t degree,
                 T x,
                 std::size_t derivative = 0,
                 long double ulps = 128) {
    const auto tol =
        derivative > degree || ulps == 0
            ? 0.L
            : std::numeric_limits<T>::epsilon() *
                  (ulps * std::max(
                              {1.L, std::abs(expected.value), expected.scale}) +
                   // Horner multiply/add steps and coordinate scaling.
                   (2 * (degree + 1) + 4) * expected.polynomial_scale);
    const bool passed =
        std::isfinite(actual) && std::isfinite(expected.value) &&
        std::isfinite(expected.scale) &&
        std::isfinite(expected.polynomial_scale) && std::isfinite(tol) &&
        std::abs(actual - expected.value) <= tol;
    assertion(passed);
    if (!passed) {
        assertion(false,
                  std::format("{} ({}), degree {}, x {}, derivative {}: "
                              "expected {}, got {}, tolerance {}",
                              name, sizeof(T) == 4 ? "float" : "double", degree,
                              x, derivative, expected.value, actual, tol));
    }
}

// Evaluate the chosen interval's polynomial, including its endpoint and
// extrapolation. Derivatives at an interior knot use the right-hand interval.
// This recursive reference does not call the production basis evaluator.
template <typename Iter>
long double basis_derivative(Iter knots,
                             std::size_t index,
                             std::size_t degree,
                             std::size_t derivative,
                             std::size_t cell,
                             long double x) {
    if (derivative > degree) { return 0; }
    if (degree == 0) { return index == cell ? 1 : 0; }
    const long double left =
        static_cast<long double>(knots[index + degree]) - knots[index];
    const long double right =
        static_cast<long double>(knots[index + degree + 1]) - knots[index + 1];
    const auto next_derivative = derivative == 0 ? 0 : derivative - 1;
    long double value{};
    if (left != 0) {
        value += (derivative == 0 ? x - knots[index] : degree) / left *
                 basis_derivative(knots, index, degree - 1, next_derivative,
                                  cell, x);
    }
    if (right != 0) {
        value += (derivative == 0 ? knots[index + degree + 1] - x
                                  : -static_cast<long double>(degree)) /
                 right *
                 basis_derivative(knots, index + 1, degree - 1, next_derivative,
                                  cell, x);
    }
    return value;
}

template <typename Spline>
std::size_t reference_cell(const Spline& spline,
                           typename Spline::knot_type& x) {
    using T = typename Spline::knot_type;
    const auto [left, right] = spline.range(0);
    if (spline.periodicity(0)) {
        // Scalar reduction isolates knot lookup from basis evaluation. Exact
        // binary-period identities test wrapping separately in the fixtures.
        x = left + std::fmod(x - left, right - left) +
            (x < left ? right - left : T{});
        if (x >= right) { x = left; }
    }
    const auto knots = spline.knots_begin(0);
    const auto first = knots + Spline::order;
    const auto end =
        std::lower_bound(first + 1, spline.knots_end(0) - Spline::order, right);
    // Select the last interval at the nonperiodic upper endpoint and retain
    // the first interval when extrapolating below the lower endpoint.
    return std::upper_bound(first + 1, end, x) - knots - 1;
}

template <typename Spline, typename ControlPoint>
ReferenceValue reference_value(const Spline& spline,
                               typename Spline::knot_type x,
                               std::size_t derivative,
                               ControlPoint control_point,
                               bool polynomial_condition = false) {
    const auto cell = reference_cell(spline, x);
    ReferenceValue result{};
    for (std::size_t j = 0; j <= Spline::order; ++j) {
        const auto index = cell - Spline::order + j;
        const auto term = control_point(index) *
                          basis_derivative(spline.knots_begin(0), index,
                                           Spline::order, derivative, cell, x);
        result.value += term;
        result.scale += std::abs(term);
    }
    if (polynomial_condition) {
        // Derive Taylor coefficients independently, without using the cached
        // production polynomial. Their absolute sum measures cancellation in
        // conversion to monomials and Horner evaluation near cell boundaries.
        const auto knots = spline.knots_begin(0);
        const long double left = knots[cell];
        long double weight = 1;
        for (std::size_t d = derivative; d <= Spline::order; ++d) {
            long double coefficient{};
            for (std::size_t j = 0; j <= Spline::order; ++j) {
                const auto index = cell - Spline::order + j;
                coefficient += control_point(index) *
                               basis_derivative(knots, index, Spline::order, d,
                                                cell, left);
            }
            result.polynomial_scale += std::abs(coefficient) * weight;
            weight *= std::abs(x - left) / (d - derivative + 1);
        }
    }
    return result;
}

}  // namespace intp_test

#endif
