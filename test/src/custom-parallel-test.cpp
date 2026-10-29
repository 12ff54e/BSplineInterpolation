#include "CustomParallelBackend.hpp"
#include "bsintp/Interpolation.hpp"

#include "Assertion.hpp"

#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#ifndef _MSC_VER
static_assert(__cplusplus == 201103L, "This test must compile in C++11 mode");
#endif

parallel_test::Backend* check_other_translation_unit();

namespace {

using Interpolation = intp::InterpolationFunction<double, 2, 3>;

double polynomial(double x, double y) {
    return x * x * x + 2. * y * y + x * y;
}

Interpolation make_interpolation() {
    intp::Mesh<double, 2> data(25, 25);
    for (std::size_t i = 0; i < 25; ++i) {
        for (std::size_t j = 0; j < 25; ++j) {
            data(i, j) = polynomial(i / 24., j / 24.);
        }
    }
    return Interpolation(data, std::make_pair(0., 1.), std::make_pair(0., 1.));
}

void check_values(Interpolation& interpolation, Assertion& assertion) {
    for (bool use_polynomial : {false, true}) {
        if (use_polynomial) { interpolation.generate_polynomial_form(); }
        for (double x : {.07, .37, .83}) {
            for (double y : {.11, .61, .93}) {
                assertion(
                    std::abs(interpolation(x, y) - polynomial(x, y)) < 1e-11,
                    "Custom backend must reproduce a cubic polynomial");
                assertion(
                    std::abs(interpolation.derivative_at({{x, y}}, {{1, 0}}) -
                             (3. * x * x + y)) < 1e-10,
                    "Custom backend must preserve analytic derivatives");
                assertion(
                    std::abs(interpolation.derivative_at({{x, y}}, {{1, 1}}) -
                             1.) < 1e-10,
                    "Custom backend must preserve mixed derivatives");
            }
        }
    }
}

void check_nested_construction(Assertion& assertion) {
    auto& backend = parallel_test::backend();
    backend.reset();
    auto parent = backend.submit([&] {
        const auto before = backend.submitted;
        auto interpolation = make_interpolation();
        assertion(backend.submitted > before,
                  "Interpolation must submit work through the custom hook");
        check_values(interpolation, assertion);
        assertion(check_other_translation_unit() == &backend,
                  "Both translation units must use the same backend");
    });
    backend.wait(parent);
    assertion(backend.submitted == backend.completed &&
                  backend.waited == backend.submitted && backend.empty(),
              "Nested construction must finish and join every task");
}

template <typename Func>
void check_failure(Func&& operation,
                   const std::string& message,
                   Assertion& assertion) {
    bool caught = false;
    try {
        operation();
    } catch (const std::runtime_error& error) {
        caught = error.what() == message;
    }
    auto& backend = parallel_test::backend();
    assertion(caught, "Parallel construction must propagate its failure");
    assertion(backend.submitted > 0 && backend.submitted == backend.completed &&
                  backend.waited == backend.submitted && backend.empty(),
              "Failure must finish every accepted task before unwinding");
}

void check_construction_failures(Assertion& assertion) {
    auto& backend = parallel_test::backend();
    backend.reset();
    backend.fail_submission = 2;
    check_failure([] { make_interpolation(); }, "submission failure",
                  assertion);

    backend.reset();
    backend.fail_execution = 1;
    check_failure([] { make_interpolation(); }, "task failure", assertion);

    // A failed coefficient build must leave a usable control-point spline.
    backend.reset();
    auto interpolation = make_interpolation();
    backend.reset();
    backend.fail_execution = 1;
    check_failure([&] { interpolation.generate_polynomial_form(); },
                  "task failure", assertion);
    assertion(
        !interpolation.spline().has_polynomial() &&
            std::abs(interpolation(.37, .61) - polynomial(.37, .61)) < 1e-11,
        "A failed coefficient build must preserve control-point values");
    backend.reset();
    check_values(interpolation, assertion);
}

void check_ranges(Assertion& assertion) {
    auto& backend = parallel_test::backend();
    for (std::size_t count : {0u, 1u, 17u, 1003u}) {
        backend.reset();
        std::vector<int> visits(count);
        intp::parallel_for(count, 64, [&](std::size_t begin, std::size_t end) {
            if (begin >= end || end > count) {
                throw std::logic_error("invalid callback range");
            }
            for (auto i = begin; i < end; ++i) { ++visits[i]; }
        });
        assertion(std::all_of(visits.begin(), visits.end(),
                              [](int n) { return n == 1; }),
                  "The backend must visit every index exactly once");
        assertion(backend.submitted == backend.completed && backend.empty(),
                  "The range must finish before returning");
        if (count == 0) {
            assertion(backend.submitted == 0,
                      "An empty range must not invoke the callback");
        }
    }
}

void check_callback_failure(Assertion& assertion) {
    auto& backend = parallel_test::backend();
    backend.reset();
    std::size_t completed_with_live_buffer = 0;
    check_failure(
        [&] {
            const auto buffer = std::make_shared<int>(42);
            const std::weak_ptr<int> borrowed = buffer;
            intp::parallel_for(
                17, 64, [&, borrowed](std::size_t begin, std::size_t end) {
                    if (begin == 0) {
                        throw std::runtime_error("callback failure");
                    }
                    const auto live = borrowed.lock();
                    if (live && *live == 42) {
                        completed_with_live_buffer += end - begin;
                    }
                });
        },
        "callback failure", assertion);
    assertion(completed_with_live_buffer == 10,
              "Callback failure must join work before destroying buffers");
}

}  // namespace

int main() {
    Assertion assertion;
    check_ranges(assertion);
    check_nested_construction(assertion);
    check_construction_failures(assertion);
    check_callback_failure(assertion);
    return assertion.status();
}
