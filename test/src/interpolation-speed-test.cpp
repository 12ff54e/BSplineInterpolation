#include <Interpolation.hpp>
#include "include/Assertion.hpp"
#include "include/Timer.h"
#include "include/rel_err.hpp"

#include <algorithm>
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

    std::mt19937 rand_gen(static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<> rand_dist(-M_PI, M_PI);
    std::uniform_real_distribution<> rand_dist2(-.5, .5);
    {
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
    constexpr size_t len_power = 6 * 4;
    constexpr size_t eval_count = 1 << 20;

    std::vector<double> eval_coord_1d;
    std::vector<std::array<double, 2>> eval_coord_2d;
    std::vector<std::array<double, 3>> eval_coord_3d;
    {
        eval_coord_1d.reserve(eval_count);
        eval_coord_2d.reserve(eval_count);
        eval_coord_3d.reserve(eval_count);
        for (size_t i = 0; i < eval_count; ++i) {
            eval_coord_1d.push_back(rand_dist(rand_gen));
            eval_coord_2d.push_back({rand_dist(rand_gen), rand_dist(rand_gen)});
            eval_coord_3d.push_back({rand_dist(rand_gen), rand_dist(rand_gen),
                                     rand_dist2(rand_gen)});
        }
    }
    // 1D case
    {
        constexpr size_t len_1d = 1 << len_power;
        constexpr double dx = 2 * M_PI / (len_1d);
        std::vector<double> vec_1d{};

        timer.start("1D Mesh");

#ifdef INTP_PERIODIC_NO_DUMMY_POINT
        vec_1d.reserve(len_1d);
        for (size_t i = 0; i < len_1d; ++i) {
            vec_1d.emplace_back(
                std::sin(13 * (static_cast<double>(i) * dx - M_PI)));
        }
#else
        vec_1d.reserve(len_1d + 1);
        for (size_t i = 0; i <= len_1d; ++i) {
            vec_1d.emplace_back(
                std::sin(13 * (static_cast<double>(i) * dx - M_PI)));
        }
#endif
        std::vector<double> vals_1d;
        vals_1d.reserve(coord_1d.size());
        for (auto& x : coord_1d) { vals_1d.push_back(std::sin(13 * x)); }

        timer.pause_and_start("1D Interpolation");

        InterpolationFunction1D<> interp1d(
            std::make_pair(-M_PI, M_PI),
            std::make_pair(vec_1d.begin(), vec_1d.end()), 3, true);

        timer.pause_and_start("1D Evaluation (Random)");

        for (size_t i = 0; i < eval_count; ++i) { interp1d(eval_coord_1d[i]); }
        // for (auto& x : eval_coord_1d) { interp1d(x); }

        timer.pause();

        std::sort(eval_coord_1d.begin(), eval_coord_1d.end());

        timer.start("1D Evaluation (Sequential)");

        for (size_t i = 0; i < eval_count; ++i) { interp1d(eval_coord_1d[i]); }
        // for (auto& x : eval_coord_1d) { interp1d(x); }

        timer.pause();

        double err_1d =
            rel_err(interp1d, std::make_pair(coord_1d.begin(), coord_1d.end()),
                    std::make_pair(vals_1d.begin(), vals_1d.end()));

        constexpr double eps = .5e-14;
        assertion(err_1d < eps);
        std::cout << "Interpolation 1d trigonometric function with err = "
                  << err_1d << '\n';

        std::cout << "Interpolation on a 1D mesh consisting " << vec_1d.size()
                  << " points. Then evaluate the function on " << eval_count
                  << " points.\n";

        timer.print();
        timer.reset();

        std::cout << '\n';
    }

    // 2d case
    {
        constexpr size_t len_2d = 1 << (len_power / 2);
        constexpr double dt = 2 * M_PI / static_cast<double>(len_2d);

        timer.start("2D Mesh");
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
        Mesh<double, 2> trig_mesh_2d(len_2d);
        for (size_t i = 0; i < len_2d; ++i) {
            for (size_t j = 0; j < len_2d; ++j) {
                trig_mesh_2d(i, j) =
                    std::sin(static_cast<double>(i) * dt - M_PI) *
                    std::cos(static_cast<double>(j) * dt - M_PI);
            }
        }
#else
        Mesh<double, 2> trig_mesh_2d(len_2d + 1);
        for (size_t i = 0; i <= len_2d; ++i) {
            for (size_t j = 0; j <= len_2d; ++j) {
                trig_mesh_2d(i, j) =
                    std::sin(static_cast<double>(i) * dt - M_PI) *
                    std::cos(static_cast<double>(j) * dt - M_PI);
            }
        }
#endif

        std::vector<double> vals_2d;
        vals_2d.reserve(coord_2d.size());

        for (auto& pt : coord_2d) {
            vals_2d.push_back(std::sin(pt[0]) * std::cos(pt[1]));
        }

        timer.pause_and_start("2D Interpolation");

        InterpolationFunction<double, 2> interp2d(3, {true, true}, trig_mesh_2d,
                                                  std::make_pair(-M_PI, M_PI),
                                                  std::make_pair(-M_PI, M_PI));

        timer.pause_and_start("2D Evaluation (Random)");

        for (size_t i = 0; i < eval_count; ++i) { interp2d(eval_coord_2d[i]); }

        timer.pause();

        std::sort(eval_coord_2d.begin(), eval_coord_2d.end(),
                  [](const std::array<double, 2>& p1,
                     const std::array<double, 2>& p2) {
                      const auto x_1 =
                          static_cast<int>(std::floor((p1[0] + M_PI) / dt));
                      const auto x_2 =
                          static_cast<int>(std::floor((p2[0] + M_PI) / dt));
                      return x_1 < x_2 || (x_1 == x_2 &&
                                           std::floor((p1[1] + M_PI) / dt) <=
                                               std::floor((p2[1] + M_PI) / dt));
                  });

        timer.start("2D Evaluation (Sequential)");

        for (size_t i = 0; i < eval_count; ++i) { interp2d(eval_coord_2d[i]); }

        timer.pause();

        double err_2d =
            rel_err(interp2d, std::make_pair(coord_2d.begin(), coord_2d.end()),
                    std::make_pair(vals_2d.begin(), vals_2d.end()));

        const double eps =
            std::pow(2 * M_PI / static_cast<double>(len_2d), 4) / 4;
        assertion(err_2d < eps);
        std::cout << "Interpolation 2d trigonometric function with err = "
                  << err_2d << '\n';

        std::cout << "Interpolation on a 2D mesh consisting "
                  << trig_mesh_2d.size()
                  << " points. Then evaluate the function on " << eval_count
                  << " points.\n";

        timer.print();
        timer.reset();

        std::cout << '\n';
    }

    // 3d case
    {
        constexpr size_t len_3d = 1 << (len_power / 3);
        constexpr double dt_3d = 2 * M_PI / len_3d;
        constexpr double dt_3d_aperiodic = 1. / (len_3d - 1);

        timer.start("3D Mesh");
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
        Mesh<double, 3> mesh_3d{len_3d, len_3d, len_3d};
        for (size_t i = 0; i < len_3d; ++i) {
            for (size_t j = 0; j < len_3d; ++j) {
                for (size_t k = 0; k < len_3d; ++k) {
                    mesh_3d(i, j, k) =
                        std::sin(static_cast<double>(i) * dt_3d - M_PI) *
                        std::cos(static_cast<double>(j) * dt_3d - M_PI) *
                        std::exp(-std::pow(
                            static_cast<double>(k) * dt_3d_aperiodic - .5, 2));
                }
            }
        }
#else
        Mesh<double, 3> mesh_3d{len_3d + 1, len_3d + 1, len_3d};
        for (size_t i = 0; i <= len_3d; ++i) {
            for (size_t j = 0; j <= len_3d; ++j) {
                for (size_t k = 0; k < len_3d; ++k) {
                    mesh_3d(i, j, k) =
                        std::sin(static_cast<double>(i) * dt_3d - M_PI) *
                        std::cos(static_cast<double>(j) * dt_3d - M_PI) *
                        std::exp(-std::pow(
                            static_cast<double>(k) * dt_3d_aperiodic - .5, 2));
                }
            }
        }
#endif

        std::vector<double> vals_3d;
        vals_3d.reserve(coord_3d.size());

        for (auto& pt : coord_3d) {
            vals_3d.push_back(std::sin(pt[0]) * std::cos(pt[1]) *
                              std::exp(-std::pow(pt[2], 2)));
        }

        timer.pause_and_start("3D Interpolation");

        InterpolationFunction<double, 3> interp3d(
            3, {true, true, false}, mesh_3d, std::make_pair(-M_PI, M_PI),
            std::make_pair(-M_PI, M_PI), std::make_pair(-.5, .5));

        timer.pause_and_start("3D Evaluation (Random)");

        for (size_t i = 0; i < eval_count; ++i) { interp3d(eval_coord_3d[i]); }

        timer.pause();

        std::sort(eval_coord_3d.begin(), eval_coord_3d.end(),
                  [](const std::array<double, 3>& p1,
                     const std::array<double, 3>& p2) {
                      const auto x_1 =
                          static_cast<int>(std::floor((p1[0] + M_PI) / dt_3d));
                      const auto x_2 =
                          static_cast<int>(std::floor((p2[0] + M_PI) / dt_3d));
                      const auto y_1 =
                          static_cast<int>(std::floor((p1[1] + M_PI) / dt_3d));
                      const auto y_2 =
                          static_cast<int>(std::floor((p2[1] + M_PI) / dt_3d));
                      return x_1 < x_2 || (x_1 == x_2 && y_1 < y_2) ||
                             (x_1 == x_2 && y_1 == y_2 &&
                              std::floor((p1[2] + .5) / dt_3d_aperiodic) <=
                                  std::floor((p2[2] + .5) / dt_3d_aperiodic));
                  });

        timer.start("3D Evaluation (Sequential)");

        for (size_t i = 0; i < eval_count; ++i) { interp3d(eval_coord_3d[i]); }

        timer.pause();

        double err_3d =
            rel_err(interp3d, std::make_pair(coord_3d.begin(), coord_3d.end()),
                    std::make_pair(vals_3d.begin(), vals_3d.end()));

        const double eps = std::pow(2 * M_PI / len_3d, 4) / 2;
        assertion(err_3d < eps);
        std::cout << "Interpolation 3d trig-exp function with err = " << err_3d
                  << '\n';

        std::cout << "Interpolation on a 3D mesh consisting " << mesh_3d.size()
                  << " points. Then evaluate the function on " << eval_count
                  << " points.\n";

        timer.print();
    }

    return assertion.status();
}
