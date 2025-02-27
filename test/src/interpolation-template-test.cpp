#include <Interpolation.hpp>
#include "include/Assertion.hpp"
#include "include/Timer.h"
#include "include/rel_err.hpp"

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

// M_PI is not part of the standard
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    using namespace intp;
    auto& timer = Timer::get_timer();

    // random sample points

    std::vector<double> coord_1d;
    std::vector<std::array<double, 2>> coord_2d;
    std::vector<std::array<double, 3>> coord_3d;
    {
        std::mt19937 rand_gen(
            static_cast<unsigned int>(std::chrono::high_resolution_clock::now()
                                          .time_since_epoch()
                                          .count()));
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

    constexpr size_t len = 1024;
    constexpr double dx = 2 * M_PI / (len * len);
    std::vector<double> trig_vec{};

    timer.start("1D Mesh");
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
    trig_vec.reserve(len * len);
    for (size_t i = 0; i < len * len; ++i) {
        trig_vec.emplace_back(
            std::sin(7 * (static_cast<double>(i) * dx - M_PI)));
    }
#else
    trig_vec.reserve(len * len + 1);
    for (size_t i = 0; i < len * len + 1; ++i) {
        trig_vec.emplace_back(
            std::sin(7 * (static_cast<double>(i) * dx - M_PI)));
    }
#endif

    std::vector<double> vals_1d;
    vals_1d.reserve(coord_1d.size());
    for (auto& x : coord_1d) { vals_1d.push_back(std::sin(7 * x)); }

    timer.pause_and_start("1D Interpolation Template");

    InterpolationFunctionTemplate1D<3> interp1d_template(
        std::make_pair(-M_PI, M_PI), trig_vec.size(), true);

    timer.pause_and_start("1D Interpolation");

    auto interp1d = interp1d_template.interpolate(
        std::make_pair(trig_vec.begin(), trig_vec.end()));

    timer.pause();

    std::cout << "Interpolation on a 1D mesh consisting " << trig_vec.size()
              << " points.\n";
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
    timer.print();
    timer.reset();
    std::cout << '\n';

    constexpr double dt = 2 * M_PI / len;

    timer.start("2D Mesh (2)");
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
    Mesh<double, 2> trig_mesh_2d_1(len);
    Mesh<double, 2> trig_mesh_2d_2(len);
    for (size_t i = 0; i < len; ++i) {
        for (size_t j = 0; j < len; ++j) {
            trig_mesh_2d_1(i, j) =
                std::sin(static_cast<double>(i) * dt - M_PI) *
                std::cos(static_cast<double>(j) * dt - M_PI);
            trig_mesh_2d_2(i, j) =
                std::cos(2 * (static_cast<double>(i) * dt - M_PI)) *
                std::sin(2 * (static_cast<double>(j) * dt - M_PI));
        }
    }
#else
    Mesh<double, 2> trig_mesh_2d_1(len + 1);
    Mesh<double, 2> trig_mesh_2d_2(len + 1);
    for (size_t i = 0; i <= len; ++i) {
        for (size_t j = 0; j <= len; ++j) {
            trig_mesh_2d_1(i, j) =
                std::sin(static_cast<double>(i) * dt - M_PI) *
                std::cos(static_cast<double>(j) * dt - M_PI);
            trig_mesh_2d_2(i, j) =
                std::cos(2 * (static_cast<double>(i) * dt - M_PI)) *
                std::sin(2 * (static_cast<double>(j) * dt - M_PI));
        }
    }
#endif
    std::vector<double> vals_2d_1;
    std::vector<double> vals_2d_2;
    vals_2d_1.reserve(coord_2d.size());
    vals_2d_2.reserve(coord_2d.size());

    for (auto& pt : coord_2d) {
        vals_2d_1.push_back(std::sin(pt[0]) * std::cos(pt[1]));
        vals_2d_2.push_back(std::cos(2 * pt[0]) * std::sin(2 * pt[1]));
    }

    timer.pause_and_start("2D Interpolation Template");

    InterpolationFunctionTemplate<double, 2, 3> interp2d_template(
        {true, true}, trig_mesh_2d_1.dimension(), std::make_pair(-M_PI, M_PI),
        std::make_pair(-M_PI, M_PI));

    timer.pause_and_start("2D Interpolation (2)");

    auto interp2d_1 = interp2d_template.interpolate(trig_mesh_2d_1);
    auto interp2d_2 = interp2d_template.interpolate(trig_mesh_2d_2);

    timer.pause();

    std::cout << "Interpolation on a 2D mesh consisting "
              << trig_mesh_2d_1.size() << " points.\n";

    double err1_2d =
        rel_err(interp2d_1, std::make_pair(coord_2d.begin(), coord_2d.end()),
                std::make_pair(vals_2d_1.begin(), vals_2d_1.end()));
    double err2_2d =
        rel_err(interp2d_2, std::make_pair(coord_2d.begin(), coord_2d.end()),
                std::make_pair(vals_2d_2.begin(), vals_2d_2.end()));
    assertion(std::max(err1_2d, err2_2d) < eps);
    if (assertion.last_status() == 0) {
        std::cout << "Interpolation 2d trigonometric function accuracy is "
                     "within volume size.\n";
    } else {
        std::cout << "Interpolation 2d trigonometric function is not accurate, "
                     "error1 = "
                  << err1_2d << '\n';
        std::cout << "Interpolation 2d trigonometric function is not accurate, "
                     "error2 = "
                  << err2_2d << '\n';
    }
    timer.print();
    timer.reset();
    std::cout << '\n';

    constexpr size_t lt = 256;
    constexpr size_t lz = 256;
    constexpr double dt_3d = 2 * M_PI / lt;
    constexpr double dz = 1. / (lz - 1);

    timer.start("3D Mesh (2)");
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
    Mesh<double, 3> mesh_3d_1(lt, lt, lz);
    Mesh<double, 3> mesh_3d_2(lt, lt, lz);
    for (size_t i = 0; i < lt; ++i) {
        for (size_t j = 0; j < lt; ++j) {
            for (size_t k = 0; k < lz; ++k) {
                mesh_3d_1(i, j, k) =
                    std::sin(static_cast<double>(i) * dt_3d - M_PI) *
                    std::cos(static_cast<double>(j) * dt_3d - M_PI) *
                    std::exp(-std::pow(static_cast<double>(k) * dz - .5, 2));
                mesh_3d_2(i, j, k) =
                    std::cos(2 * (static_cast<double>(i) * dt_3d - M_PI)) *
                    std::sin(2 * (static_cast<double>(j) * dt_3d - M_PI)) *
                    std::exp(-std::pow(static_cast<double>(k) * dz - .5, 3));
            }
        }
    }
#else
    Mesh<double, 3> mesh_3d_1(lt + 1, lt + 1, lz);
    Mesh<double, 3> mesh_3d_2(lt + 1, lt + 1, lz);
    for (size_t i = 0; i <= lt; ++i) {
        for (size_t j = 0; j <= lt; ++j) {
            for (size_t k = 0; k < lz; ++k) {
                mesh_3d_1(i, j, k) =
                    std::sin(static_cast<double>(i) * dt_3d - M_PI) *
                    std::cos(static_cast<double>(j) * dt_3d - M_PI) *
                    std::exp(-std::pow(static_cast<double>(k) * dz - .5, 2));
                mesh_3d_2(i, j, k) =
                    std::cos(2 * (static_cast<double>(i) * dt_3d - M_PI)) *
                    std::sin(2 * (static_cast<double>(j) * dt_3d - M_PI)) *
                    std::exp(-std::pow(static_cast<double>(k) * dz - .5, 3));
            }
        }
    }
#endif

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

    timer.pause_and_start("3D Interpolation Template");

    InterpolationFunctionTemplate<double, 3, 2> interp3d_template(
        {true, true, false}, mesh_3d_1.dimension(), std::make_pair(-M_PI, M_PI),
        std::make_pair(-M_PI, M_PI), std::make_pair(-.5, .5));

    timer.pause_and_start("3D Interpolation (2)");

    // InterplationFunction is default constructible
    decltype(interp3d_template)::function_type interp3d_1, interp3d_2;
    interp3d_template.interpolate(interp3d_1, mesh_3d_1);
    interp3d_template.interpolate(interp3d_2, mesh_3d_2);

    timer.pause();

    std::cout << "Interpolation on a 3D mesh consisting " << mesh_3d_1.size()
              << " points.\n";

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
    timer.print();

    return assertion.status();
}
