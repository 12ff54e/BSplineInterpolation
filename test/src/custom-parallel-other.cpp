#include "CustomParallelBackend.hpp"
#include "bsintp/BSpline.hpp"

#include <cmath>
#include <vector>

// Instantiate the same hooks and spline operations in a second translation
// unit, and return its backend address to check that the adapter is shared.
parallel_test::Backend* check_other_translation_unit() {
    using Spline = intp::BSpline<double, 2, 3>;
    Spline::ControlPointContainer cp(25, 25);
    std::fill(cp.data(), cp.data() + cp.size(), 2.);
    std::vector<double> knots(29);
    for (std::size_t i = 0; i < knots.size(); ++i) {
        knots[i] = (static_cast<double>(i) - 3.) / 22.;
    }
    auto& backend = parallel_test::backend();
    const auto before_layout = backend.submitted;
    Spline spline(cp, intp::util::get_range(knots),
                  intp::util::get_range(knots));
#ifdef INTP_CELL_LAYOUT
    if (backend.submitted == before_layout) {
        throw std::runtime_error("cell layout bypassed the custom backend");
    }
#else
    (void)before_layout;
#endif
    const auto before_polynomial = backend.submitted;
    spline.generate_polynomial_form();
    if (backend.submitted == before_polynomial) {
        throw std::runtime_error("coefficients bypassed the custom backend");
    }
    if (std::abs(spline(Spline::DimArray<double>{{.37, .61}}) - 2.) > 1e-11) {
        throw std::runtime_error("constant spline changed across backends");
    }
    return &backend;
}
