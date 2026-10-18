#include "Assertion.hpp"
#include "bsintp/Interpolation.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#ifndef _MSC_VER
static_assert(__cplusplus == 201103L, "This test must compile in C++11 mode");
#endif

namespace {

double cubic(double x, std::size_t derivative) {
    switch (derivative) {
        case 0:
            return x * x * x;
        case 1:
            return 3. * x * x;
        case 2:
            return 6. * x;
        case 3:
            return 6.;
        default:
            return 0.;
    }
}

void check_interpolation(bool uniform, Assertion& assertion) {
    std::vector<double> coords{0., .125, .25, .5, .75, .875, 1.};
    std::vector<double> values(coords.size());
    for (std::size_t i = 0; i < coords.size(); ++i) {
        if (uniform) {
            coords[i] = static_cast<double>(i) / (coords.size() - 1);
        }
        values[i] = cubic(coords[i], 0);
    }
    using Template = intp::InterpolationFunctionTemplate1D<3>;
    const auto interpolation_template =
        uniform ? Template(std::make_pair(0., 1.), coords.size())
                : Template(intp::util::get_range(coords), coords.size());
    auto interpolation =
        interpolation_template.interpolate(intp::util::get_range(values));
    for (const bool polynomial : {false, true}) {
        if (polynomial) { interpolation.generate_polynomial_form(); }
        for (const auto x : {0., .17, .5, .93, 1.}) {
            assertion(std::abs(interpolation(x) - cubic(x, 0)) < 1e-11,
                      "C++11 interpolation must reproduce a cubic");
#ifdef INTP_CELL_LAYOUT
            const auto value_proxy = interpolation.eval_proxy({{x}});
            const auto template_value_proxy =
                interpolation_template.eval_proxy(x);
            assertion(std::abs(value_proxy(interpolation.spline()) -
                               cubic(x, 0)) < 1e-11 &&
                          std::abs(template_value_proxy(interpolation) -
                                   cubic(x, 0)) < 1e-11,
                      "C++11 value proxies must reproduce a cubic");
#endif
            for (std::size_t derivative = 0; derivative <= 4; ++derivative) {
                const auto expected = cubic(x, derivative);
                assertion(std::abs(interpolation.derivative_at({{x}},
                                                               {{derivative}}) -
                                   expected) < 1e-10,
                          "C++11 derivatives must match the analytic cubic");
#ifdef INTP_CELL_LAYOUT
                const auto derivative_proxy =
                    interpolation.derivative_eval_proxy({{x}}, {{derivative}});
                const auto template_derivative_proxy =
                    interpolation_template.derivative_eval_proxy(
                        {{x}}, {{derivative}});
                assertion(
                    std::abs(derivative_proxy(interpolation.spline()) -
                             expected) < 1e-10 &&
                        std::abs(template_derivative_proxy(interpolation) -
                                 expected) < 1e-10,
                    "C++11 derivative proxies must match the analytic "
                    "cubic");
#endif
            }
        }
    }
}

void check_complex_derivatives(Assertion& assertion) {
    using Complex = std::complex<double>;
    using Spline1D = intp::BSpline<Complex, 1, 3>;
    using Spline2D = intp::BSpline<Complex, 2, 3>;
    Spline1D::ControlPointContainer cp1(4);
    Spline2D::ControlPointContainer cp2(4, 4);
    std::fill(cp1.data(), cp1.data() + cp1.size(), Complex{});
    std::fill(cp2.data(), cp2.data() + cp2.size(), Complex{});
    const Complex scale(1., 2.);
    cp1(3) = scale;
    cp2(3, 3) = scale;
    const std::vector<double> knots{0., 0., 0., 0., 1., 1., 1., 1.};
    Spline1D spline1(cp1, intp::util::get_range(knots));
    Spline2D spline2(cp2, intp::util::get_range(knots),
                     intp::util::get_range(knots));
    // These Bezier control points represent scale*x^3 and scale*x^3*y^3.
    for (const bool polynomial : {false, true}) {
        if (polynomial) {
            spline1.generate_polynomial_form();
            spline2.generate_polynomial_form();
        }
        const double x = .37;
        const double y = .61;
        for (std::size_t dx = 0; dx <= 4; ++dx) {
            const Spline1D::DimArray<std::pair<double, std::size_t>> d1{
                {{x, dx}}};
            assertion(
                std::abs(spline1.derivative_at(d1) - scale * cubic(x, dx)) <
                    1e-10,
                "C++11 complex derivatives must match the analytic cubic");
            for (std::size_t dy = 0; dy <= 4; ++dy) {
                const Spline2D::DimArray<std::pair<double, std::size_t>> d2{
                    {{x, dx}, {y, dy}}};
                assertion(std::abs(spline2.derivative_at(d2) -
                                   scale * cubic(x, dx) * cubic(y, dy)) < 1e-10,
                          "C++11 complex mixed derivatives must match the "
                          "analytic polynomial");
            }
        }
    }
}

}  // namespace

int main() {
    Assertion assertion;
    check_interpolation(false, assertion);
    check_interpolation(true, assertion);
    check_complex_derivatives(assertion);
    return assertion.status();
}
