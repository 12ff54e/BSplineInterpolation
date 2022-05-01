#include "../src/include/Interpolation.hpp"
#include "Assertion.hpp"

#include <chrono>
#include <iostream>

int main(int argc, char const* argv[]) {
    using namespace std::chrono;
    using namespace intp;

    Assertion assertion;

    const auto t_start_1d = high_resolution_clock::now();

    constexpr size_t len = 1024;
    constexpr double dx = 2 * M_PI / (len * len);
    std::vector<double> trig_vec{};
    trig_vec.reserve(len * len + 1);
    for (size_t i = 0; i < len * len + 1; ++i) {
        trig_vec.emplace_back(std::sin(7 * i * dx));
    }

    const auto t_after_vec = high_resolution_clock::now();

    InterpolationFunctionTemplate<double, 1> interp1d_template(
        3, true, trig_vec.size(), std::make_pair(0., 2 * M_PI));

    const auto t_after_template_1d = high_resolution_clock::now();

    auto interp1d = interp1d_template.interpolate(
        std::make_pair(trig_vec.begin(), trig_vec.end()));

    const auto t_after_interpolation_1d = high_resolution_clock::now();

    std::cout << "\nInterpolation on a 1D mesh consisting " << trig_vec.size()
              << " points.\n\n";
    std::cout << "Phase\t\t\tTime\n";
    std::cout << "Mesh\t\t\t"
              << duration_cast<microseconds>(t_after_vec - t_start_1d).count() /
                     1000.
              << "ms\n";
    std::cout << "Template\t\t"
              << duration_cast<microseconds>(t_after_template_1d - t_after_vec)
                         .count() /
                     1000.
              << "ms\n";
    std::cout << "Interpolate\t\t"
              << duration_cast<microseconds>(t_after_interpolation_1d -
                                             t_after_template_1d)
                         .count() /
                     1000.
              << "ms\n\n";

    double err_1d = std::abs(interp1d(M_PI_2 / 7) - 1.);
    assertion(err_1d < std::sqrt(std::numeric_limits<double>::epsilon()));
    if (assertion.last_status() == 0) {
        std::cout << "Interpolation 1d trigonometric function accuracy is "
                     "within sqrt(epsilon).\n";
    } else {
        std::cout << "Interpolation 1d trigonometric function is not "
                     "accurate, error1 = "
                  << err_1d << '\n';
    }

    const auto t_start_2d = high_resolution_clock::now();

    // In order to test the speed of bspline construction, use a 1024x1024 mesh
    constexpr double dt = 2 * M_PI / len;

    Mesh<double, 2> trig_mesh_2d_1(len + 1);
    Mesh<double, 2> trig_mesh_2d_2(len + 1);
    for (size_t i = 0; i <= len; ++i) {
        for (size_t j = 0; j <= len; ++j) {
            trig_mesh_2d_1(i, j) =
                std::sin(i * dt - M_PI) * std::cos(j * dt - M_PI);
            trig_mesh_2d_2(i, j) =
                std::sin(5 * (i * dt - M_PI)) * std::cos(5 * (j * dt - M_PI));
        }
    }

    const auto t_after_mesh = high_resolution_clock::now();

    InterpolationFunctionTemplate<double, 2> interp2d_template(
        3, {true, true}, trig_mesh_2d_1.dimension(),
        std::make_pair(-M_PI, M_PI), std::make_pair(-M_PI, M_PI));

    const auto t_after_template = high_resolution_clock::now();

    auto interp2d_1 = interp2d_template.interpolate(trig_mesh_2d_1);
    auto interp2d_2 = interp2d_template.interpolate(trig_mesh_2d_2);

    const auto t_after_interpolation = high_resolution_clock::now();

    std::cout << "\nInterpolation on a 2D mesh consisting "
              << trig_mesh_2d_1.size() << " points.\n\n";
    std::cout << "Phase\t\t\tTime\n";
    std::cout
        << "Mesh\t\t\t"
        << duration_cast<microseconds>(t_after_mesh - t_start_2d).count() /
               1000.
        << "ms\n";
    std::cout << "Template\t\t"
              << duration_cast<microseconds>(t_after_template - t_after_mesh)
                         .count() /
                     1000.
              << "ms\n";
    std::cout << "Interpolate(2)\t\t"
              << duration_cast<microseconds>(t_after_interpolation -
                                             t_after_template)
                         .count() /
                     1000.
              << "ms\n\n";

    double err1 = std::abs(interp2d_1(M_PI_2, 0.) - 1);
    double err2 = std::abs(interp2d_2(M_PI_4, M_PI_4) - .5);
    assertion(std::max(err1, err2) <
              std::sqrt(std::numeric_limits<double>::epsilon()));
    if (assertion.last_status() == 0) {
        std::cout << "Interpolation 2d trigonometric function accuracy is "
                     "within sqrt(epsilon).\n";
    } else {
        std::cout << "Interpolation 2d trigonometric function is not accurate, "
                     "error1 = "
                  << err1 << '\n';
        std::cout << "Interpolation 2d trigonometric function is not accurate, "
                     "error2 = "
                  << err2 << '\n';
    }

    return assertion.status();
}
