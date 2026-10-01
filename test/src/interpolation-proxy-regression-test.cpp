#include "Assertion.hpp"
#include "bsintp/Interpolation.hpp"

#include <array>
#include <cmath>
#include <format>
#include <stdexcept>
#include <vector>

namespace {

double polynomial(double x, std::size_t derivative) {
    switch (derivative) {
        case 0:
            return 1. - 2. * x + .5 * x * x + x * x * x;
        case 1:
            return -2. + x + 3. * x * x;
        case 2:
            return 1. + 6. * x;
        case 3:
            return 6.;
        default:
            return 0.;
    }
}

template <typename Template>
void check_polynomial(const Template& interpolation_template,
                      const std::vector<double>& coords,
                      Assertion& assertion) {
    std::vector<double> values;
    for (const auto x : coords) { values.push_back(polynomial(x, 0)); }
    const std::vector<double> other_values(8, 1.);
    const auto other = intp::InterpolationFunction1D<3>(
        std::make_pair(0., 1.), intp::util::get_range(other_values));

    for (const auto x : {0., .17, .5, .93, 1.}) {
        // A template has knots but no control points when these are created.
        const auto value_proxy = interpolation_template.eval_proxy(x);
        for (std::size_t derivative = 0; derivative <= 4; ++derivative) {
            const auto derivative_proxy =
                interpolation_template.derivative_eval_proxy({{x}},
                                                             {{derivative}});
            for (const auto scale : {1., -2.}) {
                auto scaled = values;
                for (auto& value : scaled) { value *= scale; }
                const auto interpolation = interpolation_template.interpolate(
                    intp::util::get_range(scaled));

                const auto value = value_proxy(interpolation);
                const auto deriv = derivative_proxy(interpolation);
                const auto real_value = scale * polynomial(x, 0);
                const auto real_deriv = scale * polynomial(x, derivative);
                assertion(std::abs(value - real_value) < 1e-11,
                          std::format(
                              "Template value proxy failed, expect {}, get {}",
                              value, real_value));
                assertion(
                    std::abs(deriv - real_deriv) < 1e-10,
                    std::format("Template derivative proxy failed, derivative "
                                "order {}, expect {}, get {}",
                                derivative, real_deriv, deriv));
            }
        }
        bool rejected = false;
        try {
            value_proxy(other);
        } catch (const std::invalid_argument&) { rejected = true; }
        assertion(rejected, "Value proxy must reject a different grid extent");
    }
}

}  // namespace

int main() {
    Assertion assertion;
    const std::vector<double> uniform{0.,      1. / 6., 2. / 6., .5,
                                      4. / 6., 5. / 6., 1.};
    const std::vector<double> nonuniform{0., .08, .21, .48, .62, .87, 1.};
    check_polynomial(intp::InterpolationFunctionTemplate1D<3>(
                         std::make_pair(0., 1.), uniform.size()),
                     uniform, assertion);
    check_polynomial(intp::InterpolationFunctionTemplate1D<3>(
                         intp::util::get_range(nonuniform), nonuniform.size()),
                     nonuniform, assertion);

    // Mixed periodic/nonuniform axes also exercise the even-order cell layout.
    const std::array<double, 5> y{{0., .12, .43, .76, 1.}};
    intp::Mesh<double, 2> values(8, y.size());
    for (std::size_t i = 0; i < values.dim_size(0); ++i) {
        for (std::size_t j = 0; j < y.size(); ++j) {
            values(i, j) =
                std::sin(2. * std::acos(-1.) * i / 8.) * (1. + y[j] * y[j]);
        }
    }
    const intp::InterpolationFunctionTemplate<double, 2, 2> template_2d(
        {{true, false}}, values.dimension(), std::make_pair(0., 1.),
        intp::util::get_range(y));
    for (const auto coord : {std::array<double, 2>{{0., 0.}},
                             {{.31, .37}},
                             {{.99, 1.}},
                             {{1.23, .72}}}) {
        const auto value_proxy = template_2d.eval_proxy(coord);
        for (const auto derivatives : {std::array<std::size_t, 2>{{0, 0}},
                                       {{1, 0}},
                                       {{0, 1}},
                                       {{1, 1}},
                                       {{3, 0}}}) {
            const auto derivative_proxy =
                template_2d.derivative_eval_proxy(coord, derivatives);
            const auto interpolation = template_2d.interpolate(values);
            assertion(std::abs(value_proxy(interpolation) -
                               interpolation(coord)) < 1e-12,
                      "2D template value proxy must match direct evaluation");
            assertion(std::abs(derivative_proxy(interpolation) -
                               interpolation.derivative_at(
                                   coord, derivatives)) < 1e-10,
                      "2D template derivative proxy must match direct "
                      "evaluation, including periodic wrapping");
        }
    }
    return assertion.status();
}
