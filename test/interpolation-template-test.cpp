#include "../src/include/Interpolation.hpp"
#include "./rel_err.hpp"
#include "Assertion.hpp"

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

// M_PI is not part of the standard
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    using namespace std::chrono;
    using namespace intp;

    // random sample points

    std::vector<double> coord_1d;
    std::vector<std::array<double, 2>> coord_2d;
    std::vector<std::array<double, 3>> coord_3d;
    {
        std::mt19937 rand_gen(static_cast<unsigned int>(
            high_resolution_clock::now().time_since_epoch().count()));
        std::uniform_real_distribution<> rand_dist(-M_PI, M_PI);
        std::uniform_real_distribution<> rand_dist2(-.5, .5);

        constexpr size_t coord_num = 100;
        coord_1d.reserve(coord_num);
        coord_2d.reserve(coord_num);
        coord_3d.reserve(coord_num);
        for (size_t i = 0; i < coord_num; ++i) {
            coord_1d.push_back(rand_dist(rand_gen));
            coord_2d.push_back({rand_dist(rand_gen), rand_dist(rand_gen)});
            coord_3d.push_back({rand_dist(rand_gen), rand_dist(rand_gen),
                                rand_dist2(rand_gen)});
        }
    }

    Assertion assertion;
    constexpr double eps = 1e-4;  // TODO: use a more reasonable error tolerance
                                  //  according to Taylor series expansion

    const auto t_start_1d = high_resolution_clock::now();

    constexpr size_t len = 1024;
    constexpr double dx = 2 * M_PI / (len * len);
    std::vector<double> trig_vec{};
    trig_vec.reserve(len * len + 1);
    for (size_t i = 0; i < len * len + 1; ++i) {
        trig_vec.emplace_back(std::sin(7 * (i * dx - M_PI)));
    }

    std::vector<double> vals_1d;
    vals_1d.reserve(coord_1d.size());
    for (auto& x : coord_1d) { vals_1d.push_back(std::sin(7 * x)); }

    const auto t_after_vec = high_resolution_clock::now();

    InterpolationFunctionTemplate1D<> interp1d_template(
        std::make_pair(-M_PI, M_PI), trig_vec.size(), 3, true);

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

    double err_1d =
        rel_err(interp1d, std::make_pair(coord_1d.begin(), coord_1d.end()),
                std::make_pair(vals_1d.begin(), vals_1d.end()));
    assertion(err_1d < eps);
    if (assertion.last_status() == 0) {
        std::cout << "Interpolation 1d trigonometric function accuracy is "
                     "within volume size.\n";
    } else {
        std::cout << "Interpolation 1d trigonometric function is not "
                     "accurate, error1 = "
                  << err_1d << '\n';
    }

    const auto t_start_2d = high_resolution_clock::now();

    constexpr double dt = 2 * M_PI / len;

    Mesh<double, 2> trig_mesh_2d_1(len + 1);
    Mesh<double, 2> trig_mesh_2d_2(len + 1);
    for (size_t i = 0; i <= len; ++i) {
        for (size_t j = 0; j <= len; ++j) {
            trig_mesh_2d_1(i, j) =
                std::sin(i * dt - M_PI) * std::cos(j * dt - M_PI);
            trig_mesh_2d_2(i, j) =
                std::cos(2 * (i * dt - M_PI)) * std::sin(2 * (j * dt - M_PI));
        }
    }

    std::vector<double> vals_2d_1;
    std::vector<double> vals_2d_2;
    vals_2d_1.reserve(coord_2d.size());
    vals_2d_2.reserve(coord_2d.size());

    for (auto& pt : coord_2d) {
        vals_2d_1.push_back(std::sin(pt[0]) * std::cos(pt[1]));
        vals_2d_2.push_back(std::cos(2 * pt[0]) * std::sin(2 * pt[1]));
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

    double err1 =
        rel_err(interp2d_1, std::make_pair(coord_2d.begin(), coord_2d.end()),
                std::make_pair(vals_2d_1.begin(), vals_2d_1.end()));
    double err2 =
        rel_err(interp2d_2, std::make_pair(coord_2d.begin(), coord_2d.end()),
                std::make_pair(vals_2d_2.begin(), vals_2d_2.end()));
    assertion(std::max(err1, err2) < eps);
    if (assertion.last_status() == 0) {
        std::cout << "Interpolation 2d trigonometric function accuracy is "
                     "within volume size.\n";
    } else {
        std::cout << "Interpolation 2d trigonometric function is not accurate, "
                     "error1 = "
                  << err1 << '\n';
        std::cout << "Interpolation 2d trigonometric function is not accurate, "
                     "error2 = "
                  << err2 << '\n';
    }

    const auto t_start_3d = high_resolution_clock::now();

    constexpr size_t lt = 256;
    constexpr size_t lz = 256;
    constexpr double dt_3d = 2 * M_PI / lt;
    constexpr double dz = 1. / (lz - 1);

    Mesh<double, 3> mesh_3d_1({lt + 1, lt + 1, lz});
    Mesh<double, 3> mesh_3d_2({lt + 1, lt + 1, lz});
    for (size_t i = 0; i <= lt; ++i) {
        for (size_t j = 0; j <= lt; ++j) {
            for (size_t k = 0; k < lz; ++k) {
                mesh_3d_1(i, j, k) = std::sin(i * dt_3d - M_PI) *
                                     std::cos(j * dt_3d - M_PI) *
                                     std::exp(-std::pow(k * dz - .5, 2));
                mesh_3d_2(i, j, k) = std::cos(2 * (i * dt_3d - M_PI)) *
                                     std::sin(2 * (j * dt_3d - M_PI)) *
                                     std::exp(-std::pow(k * dz - .5, 3));
            }
        }
    }

    std::vector<double> vals_3d_1;
    std::vector<double> vals_3d_2;
    vals_3d_1.reserve(coord_3d.size());
    vals_3d_2.reserve(coord_3d.size());

    for (auto& pt : coord_3d) {
        vals_3d_1.push_back(std::sin(pt[0]) * std::cos(pt[1]) *
                            std::exp(-std::pow(pt[2], 2)));
        vals_3d_2.push_back(std::cos(2 * pt[0]) * std::sin(2 * pt[1]) *
                            std::exp(-std::pow(pt[2], 3)));
    }

    const auto t_after_mesh_3d = high_resolution_clock::now();

    InterpolationFunctionTemplate<double, 3> interp3d_template(
        3, {true, true, false}, mesh_3d_1.dimension(),
        std::make_pair(-M_PI, M_PI), std::make_pair(-M_PI, M_PI),
        std::make_pair(-.5, .5));

    const auto t_after_template_3d = high_resolution_clock::now();

    auto interp3d_1 = interp3d_template.interpolate(mesh_3d_1);
    auto interp3d_2 = interp3d_template.interpolate(mesh_3d_2);

    const auto t_after_interpolation_3d = high_resolution_clock::now();

    std::cout << "\nInterpolation on a 3D mesh consisting " << mesh_3d_1.size()
              << " points.\n\n";
    std::cout << "Phase\t\t\tTime\n";
    std::cout
        << "Mesh\t\t\t"
        << duration_cast<microseconds>(t_after_mesh_3d - t_start_3d).count() /
               1000.
        << "ms\n";
    std::cout << "Template\t\t"
              << duration_cast<microseconds>(t_after_template_3d -
                                             t_after_mesh_3d)
                         .count() /
                     1000.
              << "ms\n";
    std::cout << "Interpolate(2)\t\t"
              << duration_cast<microseconds>(t_after_interpolation_3d -
                                             t_after_template_3d)
                         .count() /
                     1000.
              << "ms\n\n";

    double err1_3d =
        rel_err(interp3d_1, std::make_pair(coord_3d.begin(), coord_3d.end()),
                std::make_pair(vals_3d_1.begin(), vals_3d_1.end()));
    double err2_3d =
        rel_err(interp3d_2, std::make_pair(coord_3d.begin(), coord_3d.end()),
                std::make_pair(vals_3d_2.begin(), vals_3d_2.end()));
    assertion(std::max(err1_3d, err2_3d) < eps);
    if (assertion.last_status() == 0) {
        std::cout << "Interpolation 3d trigonometric function accuracy is "
                     "within volume size.\n";
    } else {
        std::cout << "Interpolation 3d trigonometric function is not accurate, "
                     "error1 = "
                  << err1_3d << '\n';
        std::cout << "Interpolation 3d trigonometric function is not accurate, "
                     "error2 = "
                  << err2_3d << '\n';
    }

    return assertion.status();
}
