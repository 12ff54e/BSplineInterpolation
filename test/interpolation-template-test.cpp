#include "../src/include/Interpolation.hpp"
#include "Assertion.hpp"

#include <chrono>
#include <iostream>

int main(int argc, char const* argv[]) {
    using namespace std::chrono;
    using namespace intp;

    Assertion assertion;

    const auto t_start = high_resolution_clock::now();

    // In order to test the speed of bspline construction, use a 256x256 mesh
    constexpr size_t len = 256;
    constexpr double dt = 2 * M_PI / len;

    Mesh<double, 2> trig_mesh_1(len + 1);
    Mesh<double, 2> trig_mesh_2(len + 1);
    for (size_t i = 0; i <= len; ++i) {
        for (size_t j = 0; j <= len; ++j) {
            trig_mesh_1(i, j) =
                std::sin(i * dt - M_PI) * std::cos(j * dt - M_PI);
            trig_mesh_2(i, j) =
                std::sin(5 * (i * dt - M_PI)) * std::cos(5 * (j * dt - M_PI));
        }
    }

    const auto t_after_mesh = high_resolution_clock::now();

    InterpolationFunctionTemplate<double, 2> interp_template(
        3, {true, true}, trig_mesh_1.dimension(), std::make_pair(-M_PI, M_PI),
        std::make_pair(-M_PI, M_PI));

    const auto t_after_template = high_resolution_clock::now();

    auto interp1 = interp_template.interpolate(trig_mesh_1);
    auto interp2 = interp_template.interpolate(trig_mesh_2);

    const auto t_after_interpolation = high_resolution_clock::now();

    std::cout << "Interpolation on a 2D mesh consisting " << trig_mesh_1.size()
              << " points.\n\n";
    std::cout << "Phase\t\t\tTime\n";
    std::cout << "Mesh\t\t\t"
              << duration_cast<microseconds>(t_after_mesh - t_start).count() /
                     1000.
              << "ms\n";
    std::cout << "Template\t\t"
              << duration_cast<microseconds>(t_after_template - t_after_mesh)
                         .count() /
                     1000.
              << "ms\n";
    std::cout << "Interpolate\t\t"
              << duration_cast<microseconds>(t_after_interpolation -
                                             t_after_template)
                         .count() /
                     1000.
              << "ms\n\n";

    double err1 = std::abs(interp1(M_PI_2, 0.) - 1);
    double err2 = std::abs(interp2(M_PI_4, M_PI_4) + .5);
    assertion(err1 < std::sqrt(std::numeric_limits<double>::epsilon()));
    if (assertion.last_status() == 0) {
        std::cout << "Interpolation accuracy is within sqrt(epsilon).\n";
    } else {
        std::cout << "Interpolation is not accurate, error1 = " << err1 << '\n';
        std::cout << "Interpolation is not accurate, error2 = " << err2 << '\n';
    }

    return assertion.status();
}
