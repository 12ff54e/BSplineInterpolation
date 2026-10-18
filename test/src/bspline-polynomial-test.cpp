#include "Assertion.hpp"
#include "bsintp/BSpline.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace {

void check_close(double actual, double expected, Assertion& assertion) {
    assertion(std::abs(actual - expected) < 2e-10 * (1. + std::abs(expected)),
              "Polynomial evaluation must match the B-spline basis");
}

template <typename Spline>
void check_splines(const Spline& actual,
                   const Spline& expected,
                   Assertion& assertion) {
    const std::array<double, 5> fractions{{-.125, 0., .375, 1., 1.125}};
    const intp::MeshDimension<Spline::dim> samples(fractions.size());
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto indices = samples.dimwise_indices(i);
        typename Spline::template DimArray<double> coords;
        typename Spline::template DimArray<std::pair<double, std::size_t>>
            derivatives;
        for (std::size_t d = 0; d < Spline::dim; ++d) {
            const auto range = expected.range(d);
            coords[d] = range.first +
                        fractions[indices[d]] * (range.second - range.first);
            derivatives[d] = {coords[d], 0};
        }
        check_close(actual(coords), expected(coords), assertion);
        for (std::size_t d = 0; d < Spline::dim; ++d) {
            for (std::size_t n = 1; n <= Spline::order + 1; ++n) {
                derivatives[d].second = n;
                check_close(actual.derivative_at(derivatives),
                            expected.derivative_at(derivatives), assertion);
            }
            derivatives[d].second = 0;
        }
        for (auto& derivative : derivatives) { derivative.second = 1; }
        check_close(actual.derivative_at(derivatives),
                    expected.derivative_at(derivatives), assertion);
    }
}

std::vector<double> make_knots(std::size_t count,
                               std::size_t order,
                               bool periodic,
                               bool uniform) {
    const auto cells = periodic ? count : count - order;
    std::vector<double> breaks(cells + 1);
    for (std::size_t i = 0; i < cells; ++i) {
        breaks[i + 1] = breaks[i] + (uniform ? .5 : .25 + .125 * (i % 3));
    }
    std::vector<double> knots;
    if (periodic) {
        knots.resize(count + 2 * order + 1);
        for (std::size_t i = 0; i < knots.size(); ++i) {
            const auto offset = static_cast<int>(i) - static_cast<int>(order);
            const auto n = static_cast<int>(count);
            const auto cycle = offset < 0 ? (offset - n + 1) / n : offset / n;
            knots[i] = cycle * breaks.back() + breaks[offset - cycle * n];
        }
    } else {
        knots.assign(order + 1, breaks.front());
        knots.insert(knots.end(), breaks.begin() + 1, breaks.end() - 1);
        knots.insert(knots.end(), order + 1, breaks.back());
    }
    return knots;
}

template <typename Spline, std::size_t... d>
Spline make_spline(
    const typename Spline::ControlPointContainer& cp,
    const std::array<typename Spline::KnotContainer, Spline::dim>& knots,
    typename Spline::template DimArray<bool> periodicity,
    intp::util::index_sequence<d...>) {
    return Spline(periodicity, cp, intp::util::get_range(knots[d])...);
}

template <std::size_t dim, std::size_t order>
void check_layout(Assertion& assertion) {
    using Spline = intp::BSpline<double, dim, order>;
    intp::MeshDimension<dim> dimensions;
    for (std::size_t d = 0; d < dim; ++d) {
        dimensions.dim_size(d) = order + 4 + d;
    }
    typename Spline::ControlPointContainer cp(dimensions);
    for (std::size_t i = 0; i < cp.size(); ++i) {
        cp.data()[i] = std::sin(.7 * i) + .25 * (i % 7) - .13 * (i % 11);
    }
    for (const bool uniform : {false, true}) {
        for (std::size_t mask = 0; mask < (std::size_t{1} << dim); ++mask) {
            typename Spline::template DimArray<bool> periodicity;
            std::array<typename Spline::KnotContainer, dim> knots;
            for (std::size_t d = 0; d < dim; ++d) {
                periodicity[d] = (mask & (std::size_t{1} << d)) != 0;
                knots[d] =
                    make_knots(cp.dim_size(d), order, periodicity[d], uniform);
            }
            auto spline = make_spline<Spline>(
                cp, knots, periodicity, intp::util::make_index_sequence<dim>{});
            if (uniform) {
                for (std::size_t d = 0; d < dim; ++d) { spline.set_uniform(d); }
            }
            assertion(!spline.has_polynomial(),
                      "Construction must leave polynomial coefficients empty");
            const auto reference = spline;
            spline.generate_polynomial_form();
            assertion(
                spline.has_polynomial(),
                "Explicit calculation must store polynomial coefficients");
            for (std::size_t d = 0; d < dim; ++d) {
                assertion(spline.polynomial_coef().dim_size(d) ==
                              cp.dim_size(d) - (periodicity[d] ? 0 : order),
                          "Coefficient grid must contain one block per cell");
            }
            assertion(spline.polynomial_coef().dim_size(dim) ==
                          intp::util::pow(order + 1, dim),
                      "Each cell must store every polynomial coefficient");
            check_splines(spline, reference, assertion);
        }
    }
}

void check_updates(Assertion& assertion) {
    using Spline = intp::BSpline<double, 1, 3>;
    Spline empty;
    empty.generate_polynomial_form();
    assertion(!empty.has_polynomial(),
              "An empty spline must not acquire polynomial coefficients");

    Spline::ControlPointContainer cp(16);
    for (std::size_t i = 0; i < cp.size(); ++i) { cp(i) = std::sin(.5 * i); }
    auto knots = make_knots(cp.size(), Spline::order, false, true);
    Spline spline(cp, intp::util::get_range(knots));
    spline.set_uniform(0);
    auto reference = spline;
    spline.generate_polynomial_form();
    spline.generate_polynomial_form();
    check_splines(spline, reference, assertion);

    spline.scale_control_points(-2.);
    reference.scale_control_points(-2.);
    assertion(spline.has_polynomial(),
              "Scaling must preserve the polynomial representation");
    check_splines(spline, reference, assertion);

    for (std::size_t i = 0; i < cp.size(); ++i) { cp(i) = std::cos(.3 * i); }
    spline.load_ctrlPts(cp);
    assertion(!spline.has_polynomial(),
              "Replacing control points must clear polynomial coefficients");
    reference = Spline(cp, intp::util::get_range(knots));
    check_splines(spline, reference, assertion);
    spline.generate_polynomial_form();
    check_splines(spline, reference, assertion);

    knots = make_knots(cp.size(), Spline::order, false, false);
    spline.load_knots(0, knots);
    assertion(!spline.has_polynomial(),
              "Replacing knots must clear polynomial coefficients");
    assertion(!spline.uniformity(0, 7) && spline.knot_spacing_inv(0) == 0.,
              "Replacing knots must clear uniform spacing metadata");
    reference = Spline(cp, intp::util::get_range(knots));
    check_splines(spline, reference, assertion);
    spline.generate_polynomial_form();
    check_splines(spline, reference, assertion);
}

}  // namespace

int main() {
    Assertion assertion;
    check_layout<1, 0>(assertion);
    check_layout<1, 3>(assertion);
    check_layout<2, 0>(assertion);
    check_layout<2, 2>(assertion);
    check_layout<2, 3>(assertion);
    check_layout<3, 2>(assertion);
    check_layout<3, 3>(assertion);
    check_updates(assertion);
    return assertion.status();
}
