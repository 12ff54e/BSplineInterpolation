#include "bsintp/BSpline.hpp"
#include "Assertion.hpp"
#include "SplineReference.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

#ifdef _DEBUG
#include <iomanip>
#endif

template <typename Func, typename InputIterPt, typename InputIterVal>
double rel_err(const Func& interp,
               std::pair<InputIterPt, InputIterPt> pts,
               std::pair<InputIterVal, InputIterVal> vals) {
    double err{}, l2{};
#ifdef _DEBUG
    std::cout.precision(17);
    std::cout << "\n[DEBUG] Spline Value           \tExpected\n";
#endif
    auto pt_it = pts.first;
    auto val_it = vals.first;
    for (; pt_it != pts.second && val_it != vals.second; ++pt_it, ++val_it) {
        double f = interp(*pt_it);
        err += (f - *val_it) * (f - *val_it);
        l2 += (*val_it) * (*val_it);
#ifdef _DEBUG
        std::cout << "[DEBUG] " << std::setw(20) << f << ",\t" << std::setw(20)
                  << *val_it << '\n';
#endif
    }
    return std::sqrt(err / l2);
}

template <std::size_t D>
using AlignedMesh = typename intp::BSpline<double, D, 3>::ControlPointContainer;
// order will not change type of control point container

namespace {

template <typename T, std::size_t degree, typename ControlPoint>
void check_spline_edges(const intp::BSpline<T, 1, degree, T>& spline,
                        const std::vector<T>& points,
                        ControlPoint control_point,
                        Assertion& assertion) {
    using namespace intp_test;
    const auto last = spline.knots_num(0) - degree - 2;
    for (const auto x : points) {
        auto expected_x = x;
        const auto cell = reference_cell(spline, expected_x);
        // A deliberately stale hint must still find the stored knot interval.
        for (const auto hint : {degree, last}) {
            auto actual_x = x;
            const auto iter = spline.get_knot_iter(0, actual_x, hint);
            assertion(iter - spline.knots_begin(0) ==
                              static_cast<std::ptrdiff_t>(cell) &&
                          actual_x == expected_x,
                      std::format("Knot lookup, degree {}, x {}, hint {}",
                                  degree, x, hint));
        }
        for (std::size_t derivative = 0; derivative <= degree + 1;
             ++derivative) {
            const auto expected =
                reference_value(spline, x, derivative, control_point);
            if (derivative <= degree) {
                const auto basis =
                    spline.base_spline_value(0, spline.knots_begin(0) + cell,
                                             expected_x, degree - derivative);
                ReferenceValue sum{};
                for (std::size_t j = 0; j <= degree; ++j) {
                    const auto value = basis_derivative(
                        spline.knots_begin(0), cell - degree + j, degree,
                        derivative, cell, expected_x);
                    check_close(assertion, basis[j], {value, std::abs(value)},
                                "Individual B-spline basis", degree, x,
                                derivative);
                    sum.value += basis[j];
                    sum.scale += std::abs(basis[j]);
                }
                check_close(assertion, static_cast<T>(sum.value),
                            {derivative == 0 ? 1.L : 0.L, sum.scale},
                            "Basis partition of unity", degree, x, derivative);
            }
            check_close(assertion, spline.derivative_at({{{x, derivative}}}),
                        expected, "B-spline derivative", degree, x, derivative);
            if (derivative == 0) {
                check_close(assertion, spline({{{x, degree}}}), expected,
                            "B-spline value", degree, x);
            }
        }
    }
#ifdef INTP_CELL_LAYOUT
    // Keep proxy failures separate so they do not hide the direct checks.
    try {
        for (const auto x : points) {
            for (std::size_t derivative = 0; derivative <= degree + 1;
                 ++derivative) {
                const auto proxy =
                    spline.pre_calc_coef_helper({{{x, degree, derivative}}});
                check_close(
                    assertion, proxy(spline),
                    reference_value(spline, x, derivative, control_point),
                    "B-spline proxy", degree, x, derivative);
            }
        }
    } catch (const std::exception& e) {
        assertion(false,
                  std::format("B-spline proxy ({}), degree {}, periodic {}, "
                              "uniform {}: {}",
                              sizeof(T) == 4 ? "float" : "double", degree,
                              spline.periodicity(0),
                              spline.uniformity(0, degree), e.what()));
    }
#endif
}

template <typename T, std::size_t degree>
void check_bspline_edges(Assertion& assertion) {
    using Spline = intp::BSpline<T, 1, degree, T>;
    const auto infinity = std::numeric_limits<T>::infinity();
    const auto add_neighbors = [&](std::vector<T>& points, T x) {
        points.insert(points.end(), {std::nextafter(x, -infinity), x,
                                     std::nextafter(x, infinity)});
    };

    // Alternating Bezier control points give the exact polynomial (1-2x)^p.
    std::vector<T> knots(2 * (degree + 1), T{});
    std::fill(knots.begin() + degree + 1, knots.end(), T{1});
    typename Spline::ControlPointContainer cp(degree + 1);
    for (std::size_t i = 0; i <= degree; ++i) { cp(i) = i % 2 ? T{-1} : T{1}; }
    Spline bezier(cp, intp::util::get_range(knots));
    std::vector<T> points{T{-.25}, T{.25}, T{.5}, T{.75}, T{1.25}};
    add_neighbors(points, T{});
    add_neighbors(points, T{1});
    for (const auto x : points) {
        for (std::size_t derivative = 0; derivative <= degree + 1;
             ++derivative) {
            long double expected = derivative > degree ? 0 : 1;
            if (derivative <= degree) {
                for (std::size_t i = 0; i < derivative; ++i) {
                    expected *= -2.L * (degree - i);
                }
                for (std::size_t i = derivative; i < degree; ++i) {
                    expected *= 1.L - 2.L * x;
                }
            }
            const auto reference = intp_test::reference_value(
                bezier, x, derivative, [&](std::size_t i) { return cp(i); });
            // The derivative can vanish while its basis terms are large.
            intp_test::check_close(
                assertion, bezier.derivative_at({{{x, derivative}}}),
                {expected, reference.scale}, "Bezier analytic derivative",
                degree, x, derivative);
            assertion(std::abs(reference.value - expected) <=
                          128 * std::numeric_limits<long double>::epsilon() *
                              std::max(1.L, reference.scale),
                      "Recursive reference must match the Bezier polynomial");
        }
    }
    check_spline_edges(
        bezier, points, [&](std::size_t i) { return cp(i); }, assertion);
    bezier.generate_polynomial_form();
    check_spline_edges(
        bezier, points, [&](std::size_t i) { return cp(i); }, assertion);

    if constexpr (degree > 0) {
        // Multiplicity p leaves only C0 continuity at the interior knot.
        knots.assign(degree + 1, T{});
        knots.insert(knots.end(), degree, T{.5});
        knots.insert(knots.end(), degree + 1, T{1});
        cp = typename Spline::ControlPointContainer(2 * degree + 1);
        for (std::size_t i = 0; i < cp.size(); ++i) {
            cp(i) = static_cast<T>(i % 2 ? -2 : 1);
        }
        Spline repeated(cp, intp::util::get_range(knots));
        points.clear();
        for (const T x : {T{}, T{.5}, T{1}}) { add_neighbors(points, x); }
        check_spline_edges(
            repeated, points, [&](std::size_t i) { return cp(i); }, assertion);
        repeated.generate_polynomial_form();
        check_spline_edges(
            repeated, points, [&](std::size_t i) { return cp(i); }, assertion);
    }

    constexpr std::size_t count = 7;
    for (const bool uniform : {false, true}) {
        std::array<T, count + 1> breaks{};
        for (std::size_t i = 0; i < count; ++i) {
            breaks[i + 1] =
                breaks[i] + (uniform ? T{.5} : (i % 2 ? T{.75} : T{1.25}));
        }
        const auto period = breaks.back();
        knots.resize(count + 2 * degree + 1);
        for (std::size_t i = 0; i < knots.size(); ++i) {
            const auto offset = static_cast<int>(i) - static_cast<int>(degree);
            const auto cycle = offset < 0
                                   ? (offset - static_cast<int>(count) + 1) /
                                         static_cast<int>(count)
                                   : offset / static_cast<int>(count);
            knots[i] = T{-2.25} + cycle * period +
                       breaks[offset - cycle * static_cast<int>(count)];
        }
        cp = typename Spline::ControlPointContainer(count);
        for (std::size_t i = 0; i < count; ++i) {
            cp(i) = static_cast<T>((i % 2 ? -1 : 2) + .25 * i);
        }
        Spline periodic({true}, cp, intp::util::get_range(knots));
        if (uniform) { periodic.set_uniform(0); }
        points.clear();
        for (std::size_t i = degree; i <= degree + count; ++i) {
            add_neighbors(points, knots[i]);
            for (const int shift : {-3, -1, 1, 3}) {
                points.push_back(knots[i] + shift * period);
            }
        }
        for (const T x : {T{-2}, T{-.25}}) {
            for (const int shift : {-3, -1, 0, 1, 3}) {
                intp_test::check_close(
                    assertion, periodic({{{x + shift * period, degree}}}),
                    {periodic({{{x, degree}}}), 1}, "Exact binary period shift",
                    degree, x);
            }
        }
        check_spline_edges(
            periodic, points, [&](std::size_t i) { return cp(i % count); },
            assertion);
        periodic.generate_polynomial_form();
        check_spline_edges(
            periodic, points, [&](std::size_t i) { return cp(i % count); },
            assertion);
    }
}

template <typename T, std::size_t... degree>
void check_bspline_orders(Assertion& assertion,
                          std::index_sequence<degree...>) {
    (check_bspline_edges<T, degree>(assertion), ...);
}

}  // namespace

int main() {
    using namespace std;
    using namespace intp;

    Assertion assertion;
    constexpr double tol = 1e-14;

    check_bspline_orders<float>(assertion, std::make_index_sequence<6>{});
    check_bspline_orders<double>(assertion, std::make_index_sequence<6>{});

    // 1D uniform test
    {
        std::array<double, 16> knots;
        std::array<double, 12> cp;
        for (std::size_t i = 0; i < knots.size(); ++i) {
            knots[i] = static_cast<double>(i < 4 ? 0 : i > 11 ? 9 : i - 3);
        }
        cp.fill(1.);

        BSpline<double, 1, 3> spline_1d_3(
            AlignedMesh<1>(std::make_pair(cp.begin(), cp.end())),
            std::make_pair(knots.begin(), knots.end()));

        spline_1d_3.set_uniform(0);

        std::array<double, 91> coords_1d;
        std::array<double, 91> vals_1d;
        for (std::size_t i = 0; i < coords_1d.size(); ++i) {
            coords_1d[i] = static_cast<double>(i) / 10.;
            vals_1d[i] = 1.;
        }
        std::cout << "\nB-Spline Base Test:\n";
        double d = rel_err([&](double x) { return spline_1d_3({x}); },
                           std::make_pair(coords_1d.begin(), coords_1d.end()),
                           std::make_pair(vals_1d.begin(), vals_1d.end()));

        assertion(d < tol);
        std::cout << "\nBase test "
                  << (assertion.last_status() == 0 ? "succeed" : "failed")
                  << '\n';
        std::cout << "Relative Error = " << d << '\n';
    }
    // 1D B-Spline

    auto knots = {0., 0., 0., 0., .5, 1., 1., 1., 1.};
    auto cp = {1., 16. / 3, -16. / 3, 4. / 3, 3.};

    BSpline<double, 1, 3> spline_1d_3(
        AlignedMesh<1>(std::make_pair(cp.begin(), cp.end())),
        std::make_pair(knots.begin(), knots.end()));

    // some random points
    auto coords_1d = {0.4560373422725581,  0.8888069703323336,
                      0.04461098353789361, 0.5439456019541871,
                      0.532230619202863,   0.5153105562793239,
                      0.21346152704976218, 0.08993837286502737,
                      0.20928777742014026, 0.11139302362905723};
    // Values pre-computed by MMA
    auto vals_1d = {-0.6452391969741351, 1.706781268175239,
                    1.9389729439006869,  -1.1673240687437252,
                    -1.1409725647301816, -1.079846126512429,
                    2.3537378042387784,  2.4815644189239574,
                    2.387237177311428,   2.6116510526724888};

    std::cout << "\n1D B-Spline Test:\n";
    double d = rel_err([&](double x) { return spline_1d_3({x}); },
                       std::make_pair(coords_1d.begin(), coords_1d.end()),
                       std::make_pair(vals_1d.begin(), vals_1d.end()));

    assertion(d < tol);
    std::cout << "\n1D test "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

    // 2D B-Spline

    // random 5x5 control points
    constexpr array<array<double, 5>, 5> cp2{{{-4, -4, 4, 1, 1},
                                              {0, 2, -2, 0, 5},
                                              {5, -4, -4, -3, -5},
                                              {-2, 4, 2, 2, -3},
                                              {1, 0, 5, -3, 5}}};

    AlignedMesh<2> cp2d{5, 5};
    for (unsigned i = 0; i < 5; ++i) {
        for (unsigned j = 0; j < 5; ++j) { cp2d(i, j) = cp2[i][j]; }
    }

    BSpline<double, 2, 3> spline_2d_3(cp2d,
                                      make_pair(knots.begin(), knots.end()),
                                      make_pair(knots.begin(), knots.end()));

    std::cout << "\n2D B-Spline Test:\n";

    // some random points
    constexpr array<array<double, 2>, 10> coords_2d{
        {{0.049800826265867126, 0.8612651106123255},
         {0.575715550141984, 0.5153782982448287},
         {0.051389608765929795, 0.20636530700031974},
         {0.3406193462245328, 0.6559088769326575},
         {0.6562156301169633, 0.5229827075599092},
         {0.7065591320791091, 0.18017930762978818},
         {0.5929532206588517, 0.015110486908892717},
         {0.35737449937918164, 0.2771152465740585},
         {0.04848381891207554, 0.9244449256192726},
         {0.9740876351532088, 0.9174673912595357}}};

    // and their spline value, pre-computed by MMA
    constexpr array<double, 10> vals_2d{
        {1.3002562089978034, -0.9049532347293099, -1.5715077611939692,
         -1.3555951852638821, -0.19537752208751952, 0.7893147505912226,
         1.3459399154210203, -0.7670721516764836, 1.495128187689829,
         1.5162241194068322}};

    d = rel_err(
        [&](const std::array<double, 2>& coord) { return spline_2d_3(coord); },
        std::make_pair(coords_2d.begin(), coords_2d.end()),
        std::make_pair(vals_2d.begin(), vals_2d.end()));
    assertion(d < tol);
    std::cout << "\n2D test "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

    // 3D Spline

    // random 5x5x5 control points
    constexpr array<array<array<double, 5>, 5>, 5> cp3{{{{{-4, -4, -4, -2, -5},
                                                          {-4, 0, -5, 5, 4},
                                                          {4, -5, -3, 5, -4},
                                                          {3, 4, 0, -3, -1},
                                                          {-5, 2, 1, -5, -2}}},
                                                        {{{5, 4, -3, -5, 5},
                                                          {2, -3, 0, 1, 4},
                                                          {5, 0, -4, 3, -1},
                                                          {2, -2, 3, 0, 2},
                                                          {-1, 1, -4, -3, -2}}},
                                                        {{{-1, 3, -1, 2, 2},
                                                          {-1, -5, -3, -5, -1},
                                                          {3, -3, -5, 2, 0},
                                                          {2, 1, 1, 3, 3},
                                                          {2, 2, 4, 4, -1}}},
                                                        {{{-4, 0, 2, -5, 2},
                                                          {-2, -3, 1, 2, -2},
                                                          {4, -4, 5, -3, 0},
                                                          {-4, -5, 3, -4, 2},
                                                          {2, -4, -2, 1, 5}}},
                                                        {{{-3, 2, 1, 0, -2},
                                                          {2, 3, 2, -5, 5},
                                                          {-1, 2, 3, 0, -5},
                                                          {1, 2, -1, 1, -1},
                                                          {0, -1, 5, 5, 4}}}}};

    AlignedMesh<3> cp3d{5, 5, 5};
    for (unsigned i = 0; i < 5; ++i) {
        for (unsigned j = 0; j < 5; ++j) {
            for (unsigned k = 0; k < 5; ++k) { cp3d(i, j, k) = cp3[i][j][k]; }
        }
    }

    BSpline<double, 3, 3> spline_3d_3(cp3d,
                                      make_pair(knots.begin(), knots.end()),
                                      make_pair(knots.begin(), knots.end()),
                                      make_pair(knots.begin(), knots.end()));

    std::cout << "\n3D B-Spline Test:\n";

    // some random points
    constexpr array<array<double, 3>, 20> coords_3d{
        {{0.6771169667756889, 0.9215782737838114, 0.4073093854149159},
         {0.5587048417879672, 0.27660875852915034, 0.480940182622974},
         {0.22166100919838327, 0.20423034083657687, 0.6828162120660075},
         {0.946117300279044, 0.6987537482071375, 0.988288195761676},
         {0.6657135834502097, 0.6362329447266486, 0.6098637234378605},
         {0.4831720672055162, 0.04449885480903437, 0.10794187260224075},
         {0.01970883141191182, 0.707507573262028, 0.3842127548987244},
         {0.2708700179218677, 0.9968009386232979, 0.7784159674453355},
         {0.34911366851070724, 0.37498001215861887, 0.5313335523213891},
         {0.7831337922295372, 0.3523807977438491, 0.5616443755053859},
         {0.5651837447297534, 0.3521803904983014, 0.3766317887491841},
         {0.42938234954177146, 0.03216063762167898, 0.2257110836573426},
         {0.5320732583223098, 0.2563335748568327, 0.4284182446681921},
         {0.27509099660894876, 0.3977734281932175, 0.7406749938043244},
         {0.2645390059243993, 0.8727863468285368, 0.8960358519530291},
         {0.8338988630403117, 0.2793369171836766, 0.6474522893994397},
         {0.01796573121296885, 0.8821327182314185, 0.023056755634403236},
         {0.29524574329222886, 0.7820003453674633, 0.1953107714234823},
         {0.27122873705368566, 0.050574347888162174, 0.8570032085001116},
         {0.3328181108995292, 0.22075349315084103, 0.7955065113795101}}};

    // and their spline value, pre-computed by MMA
    constexpr array<double, 20> vals_3d{
        {-0.1680401055752912,  -1.5184595156006049,  -0.6181298890131312,
         -0.7634699499420325,  -0.11071964593272185, 0.36484761815950045,
         -0.13606044311702267, -1.066384805143148,   -1.3030396801520705,
         -0.1366609994525347,  -1.845924629958203,   0.7750141662857806,
         -1.7374305394795615,  0.1278576545692187,   0.12587781246729338,
         -0.16274434441451202, -0.04873351906943385, 0.11126350145310085,
         -0.24228756896767534, -0.35239496660939207}};

    d = rel_err(
        [&](const std::array<double, 3>& coord) {
            return spline_3d_3({coord[0], coord[1], coord[2]});
        },
        std::make_pair(coords_3d.begin(), coords_3d.end()),
        std::make_pair(vals_3d.begin(), vals_3d.end()));
    assertion(d < tol);
    std::cout << "\n3D test "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

    // 2D with one dimension being periodic

    std::cout << "\n2D B-Spline with periodic boundary Test:\n";

    auto knots2 = {-0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6};
    BSpline<double, 2, 3> spline_2d_3_periodic(
        {false, true}, cp2d, make_pair(knots.begin(), knots.end()),
        make_pair(knots2.begin(), knots2.end()));

    constexpr array<double, 10> vals_2d_periodic{
        -2.4226917831367354, -1.6666970964438268, 1.292431192225407,
        0.7418964927133986,  -1.6357978631810504, 0.32021767652671507,
        -0.0198531979338478, -1.6728476167273418, -2.414215484638793,
        0.7421859678077125};

    d = rel_err(
        [&](const array<double, 2>& coord) {
            return spline_2d_3_periodic(coord);
        },
        std::make_pair(coords_2d.begin(), coords_2d.end()),
        std::make_pair(vals_2d_periodic.begin(), vals_2d_periodic.end()));
    assertion(d < tol);
    std::cout << "\n2D test with periodic boundary "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

    // 1D derivative

    std::cout << "\n1D B-Spline derivative Test:\n";

    auto vals_1d_derivative_1 = {-9.92272407143533,   12.667321048785523,
                                 16.318938809554005,  -1.7077524908452282,
                                 -2.7976046481470465, -4.4414576381861774,
                                 -8.212957569324931,  7.852171562286172,
                                 -7.837502772245249,  4.326042837701884};
    auto vals_1d_derivative_2 = {74.45709400715904,   -7.976407455712078,
                                 -202.02141906253547, 91.34366663719412,
                                 94.71758166957548,   99.59055979155474,
                                 -88.5538538225598,   -171.56141343470156,
                                 -91.35861357366572,  -157.14388812127353};
    d = rel_err(
        [&](double x) {
            return spline_1d_3.derivative_at({std::make_pair(x, 0)});
        },
        std::make_pair(coords_1d.begin(), coords_1d.end()),
        std::make_pair(vals_1d.begin(), vals_1d.end()));

    assertion(d < tol);
    std::cout << "\n1D test of derivative 0 (fallback to spline value) "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

    d = rel_err(
        [&](double x) {
            return spline_1d_3.derivative_at({std::make_pair(x, 1)});
        },
        std::make_pair(coords_1d.begin(), coords_1d.end()),
        std::make_pair(vals_1d_derivative_1.begin(),
                       vals_1d_derivative_1.end()));
    assertion(d < tol);
    std::cout << "\n1D test of derivative 1 "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

    d = rel_err(
        [&](double x) {
            return spline_1d_3.derivative_at({std::make_pair(x, 2)});
        },
        std::make_pair(coords_1d.begin(), coords_1d.end()),
        std::make_pair(vals_1d_derivative_2.begin(),
                       vals_1d_derivative_2.end()));
    assertion(d < tol);
    std::cout << "\n1D test of derivative 2 "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

    // 2D partial derivative

    std::cout << "\n2D B-Spline derivative Test:\n";

    auto vals_2d_derivative_x2_y0 = {-64.18140090609639,  33.606993827386646,
                                     -97.40128619523041,  40.27561307446061,
                                     9.28340750029671,    -7.934058734360374,
                                     -24.125345223515264, 14.846578817272864,
                                     -113.9333016161396,  30.79197201901215};
    auto vals_2d_derivative_x1_y1 = {88.74503910110498,  -4.238681068935503,
                                     -81.71582304860588, -7.40376361081399,
                                     -8.066246965615212, 47.312250959237815,
                                     91.39055060903914,  -10.267520769902886,
                                     106.54723295116875, 253.25136807731457};

    d = rel_err(
        [&](array<double, 2> coord) {
            return spline_2d_3.derivative_at(
                {std::make_pair(coord[0], 2), std::make_pair(coord[1], 0)});
        },
        std::make_pair(coords_2d.begin(), coords_2d.end()),
        std::make_pair(vals_2d_derivative_x2_y0.begin(),
                       vals_2d_derivative_x2_y0.end()));

    assertion(d < tol);
    std::cout << "\n2D test of derivative 2,0 "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

    d = rel_err(
        [&](array<double, 2> coord) {
            return spline_2d_3.derivative_at(
                {std::make_pair(coord[0], 1), std::make_pair(coord[1], 1)});
        },
        std::make_pair(coords_2d.begin(), coords_2d.end()),
        std::make_pair(vals_2d_derivative_x1_y1.begin(),
                       vals_2d_derivative_x1_y1.end()));

    assertion(d < tol);
    std::cout << "\n2D test of derivative 1,1 "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

    std::cout << "\n2D B-Spline with periodic boundary derivative Test:\n";

    auto vals_2d_periodic_derivative_x1_y1 = {
        10.419194997422988, -44.893624896638144, -67.96918096244443,
        64.25725539946846,  -24.772206659332948, -5.667565671444081,
        39.05723974033617,  -5.748203566140304,  -36.99426057331604,
        -104.22161375710107};

    d = rel_err(
        [&](array<double, 2> coord) {
            return spline_2d_3_periodic.derivative_at(
                {std::make_pair(coord[0], 1), std::make_pair(coord[1], 1)});
        },
        std::make_pair(coords_2d.begin(), coords_2d.end()),
        std::make_pair(vals_2d_periodic_derivative_x1_y1.begin(),
                       vals_2d_periodic_derivative_x1_y1.end()));

    assertion(d < tol);
    std::cout << "\n2D test of derivative 1,1 with periodic boundary "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

#ifdef INTP_CELL_LAYOUT
    d = rel_err(
        [&](double x) {
            auto evaluator = spline_1d_3.pre_calc_coef({{{x, 3u}}});
            return evaluator(spline_1d_3);
        },
        std::make_pair(coords_1d.begin(), coords_1d.end()),
        std::make_pair(vals_1d.begin(), vals_1d.end()));
#endif
    assertion(d < tol);
    std::cout << "\nProxy test "
              << (assertion.last_status() == 0 ? "succeed" : "failed") << '\n';
    std::cout << "Relative Error = " << d << '\n';

    return assertion.status();
}
