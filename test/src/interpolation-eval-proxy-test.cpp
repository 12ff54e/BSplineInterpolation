#include "Interpolation.hpp"
#include "include/Assertion.hpp"
#include "include/Timer.h"
#include "include/rel_err.hpp"

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

// M_PI is not part of the standard
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    Assertion assertion;
#ifndef INTP_CELL_LAYOUT
    assertion(0);
    std::cout << "Evaluation proxy needs `INTP_CELL_LAYOUT`.\n";
#else

    using namespace intp;
    auto& timer = Timer::get_timer();

    // random sample points

    std::mt19937 rand_gen(static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<> rand_dist(-M_PI, M_PI);
    std::uniform_real_distribution<> rand_dist2(-.5, .5);

    constexpr size_t len_power = 6 * 2;
    constexpr size_t eval_count = 1 << 10;
    constexpr size_t repeat_time = 1 << 10;
    constexpr double delta = 1.e-2 * M_PI;

    constexpr size_t len_1d = 1 << len_power;
    constexpr size_t len_2d = 1 << (len_power / 2);
    constexpr size_t len_3d = 1 << (len_power / 3);

    constexpr double dx = 2 * M_PI / (len_1d);
    constexpr double dt_2d = 2 * M_PI / static_cast<double>(len_2d);
    constexpr double dt_3d = 2 * M_PI / len_3d;
    constexpr double dt_3d_aperiodic = 1. / (len_3d - 1);

    std::vector<double> eval_coord_1d;
    std::vector<std::array<double, 2> > eval_coord_2d;
    std::vector<std::array<double, 3> > eval_coord_3d;
    eval_coord_1d.reserve(eval_count);
    eval_coord_2d.reserve(eval_count);
    eval_coord_3d.reserve(eval_count);
    for (size_t i = 0; i < eval_count; ++i) {
        eval_coord_1d.push_back(rand_dist(rand_gen));
        eval_coord_2d.push_back({rand_dist(rand_gen), rand_dist(rand_gen)});
        eval_coord_3d.push_back(
            {rand_dist(rand_gen), rand_dist(rand_gen), rand_dist2(rand_gen)});
    }
    std::sort(eval_coord_1d.begin(), eval_coord_1d.end());
    std::sort(
        eval_coord_2d.begin(), eval_coord_2d.end(),
        [](const std::array<double, 2>& p1, const std::array<double, 2>& p2) {
            const auto x_1 =
                static_cast<int>(std::floor((p1[0] + M_PI) / dt_2d));
            const auto x_2 =
                static_cast<int>(std::floor((p2[0] + M_PI) / dt_2d));
            return x_1 < x_2 ||
                   (x_1 == x_2 && std::floor((p1[1] + M_PI) / dt_2d) <=
                                      std::floor((p2[1] + M_PI) / dt_2d));
        });
    std::sort(
        eval_coord_3d.begin(), eval_coord_3d.end(),
        [](const std::array<double, 3>& p1, const std::array<double, 3>& p2) {
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

    std::vector<double> diff(eval_count);

    // 1D case
    {
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
        constexpr std::size_t data_len = len_1d;
#else
        constexpr std::size_t data_len = len_1d + 1;
#endif
        std::vector<double> vec_1d(data_len);

        InterpolationFunctionTemplate1D<> interp1d_template(
            std::make_pair(-M_PI, M_PI), vec_1d.size(), 3, true);
        std::vector<decltype(interp1d_template.eval_proxy(0.))> evaluators;

        timer.start("1D Evaluator");
        for (std::size_t i = 0; i < eval_count; ++i) {
            evaluators.push_back(
                interp1d_template.eval_proxy(eval_coord_1d[i]));
        }
        timer.pause();

        for (std::size_t t = 0; t < repeat_time; ++t) {
            timer.start("1D Mesh");
            for (std::size_t i = 0; i < data_len; ++i) {
                vec_1d[i] = std::sin(13 * (static_cast<double>(i) * dx - M_PI) +
                                     static_cast<double>(t) * delta);
            }

            timer.pause_and_start("1D Interpolation");
            auto interp1d = interp1d_template.interpolate(
                std::make_pair(vec_1d.begin(), vec_1d.end()));

            timer.pause_and_start("1D Evaluation (Direct)");
            for (std::size_t i = 0; i < eval_count; ++i) {
                diff[i] = interp1d(eval_coord_1d[i]);
            }
            timer.pause_and_start("1D Evaluation (Use proxy)");
            for (std::size_t i = 0; i < eval_count; ++i) {
                diff[i] -= evaluators[i](interp1d);
            }
            timer.pause();
        }
        std::cout << "Interpolation on a 1D mesh consisting " << vec_1d.size()
                  << " points. Then evaluate the function on " << eval_count
                  << " points. Repeat for " << repeat_time << " times.\n";
        double diff_L2{};
        for (auto v : diff) { diff_L2 += v * v; }
        diff_L2 = std::sqrt(diff_L2 / eval_count);
        std::cout << "Difference between two methods: " << diff_L2 << '\n';
        assertion(diff_L2 < std::numeric_limits<double>::epsilon());

        timer.print();
        timer.reset();
        std::cout << '\n';
    }

    // 2d case
    {
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
        constexpr std::size_t data_len = len_2d;
#else
        constexpr std::size_t data_len = len_2d + 1;
#endif
        Mesh<double, 2> trig_mesh_2d(data_len);

        InterpolationFunctionTemplate<double, 2> interp2d_template(
            3, {true, true}, trig_mesh_2d.dimension(),
            std::make_pair(-M_PI, M_PI), std::make_pair(-M_PI, M_PI));

        timer.start("2D Evaluator");
        std::vector<decltype(interp2d_template.eval_proxy({0., 0.}))>
            evaluators;

        for (std::size_t i = 0; i < eval_count; ++i) {
            evaluators.push_back(
                interp2d_template.eval_proxy(eval_coord_2d[i]));
        }
        timer.pause();

        for (std::size_t t = 0; t < repeat_time; ++t) {
            timer.start("2D Mesh");
            for (std::size_t i = 0; i < data_len; ++i) {
                for (std::size_t j = 0; j < data_len; ++j) {
                    trig_mesh_2d(i, j) =
                        std::sin(static_cast<double>(i) * dt_2d - M_PI +
                                 static_cast<double>(t) * delta) *
                        std::cos(static_cast<double>(j) * dt_2d - M_PI);
                }
            }

            timer.pause_and_start("2D Interpolation");
            auto interp2d = interp2d_template.interpolate(trig_mesh_2d);

            timer.pause_and_start("2D Evaluation (Direct)");
            for (std::size_t i = 0; i < eval_count; ++i) {
                diff[i] = interp2d(eval_coord_2d[i]);
            }

            timer.pause_and_start("2D Evaluation (Use proxy)");
            for (std::size_t i = 0; i < eval_count; ++i) {
                diff[i] -= evaluators[i](interp2d);
            }
            timer.pause();
        }
        std::cout << "Interpolation on a 2D mesh consisting "
                  << trig_mesh_2d.size()
                  << " points. Then evaluate the function on " << eval_count
                  << " points. Repeat for " << repeat_time << " times.\n";
        double diff_L2{};
        for (auto v : diff) { diff_L2 += v * v; }
        diff_L2 = std::sqrt(diff_L2 / eval_count);
        std::cout << "Difference between two methods: " << diff_L2 << '\n';
        assertion(diff_L2 < std::numeric_limits<double>::epsilon());

        timer.print();
        timer.reset();
        std::cout << '\n';
    }

    // 3d case
    {
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
        constexpr std::size_t data_len = len_3d;
#else
        constexpr std::size_t data_len = len_3d + 1;
#endif

        Mesh<double, 3> mesh_3d{data_len, data_len, len_3d};
        InterpolationFunctionTemplate<double, 3> interp3d_template(
            3, {true, true, false}, mesh_3d.dimension(),
            std::make_pair(-M_PI, M_PI), std::make_pair(-M_PI, M_PI),
            std::make_pair(-.5, .5));

        timer.start("3D Evaluator");
        std::vector<decltype(interp3d_template.eval_proxy({0., 0., 0.}))>
            evaluators;
        for (std::size_t i = 0; i < eval_count; ++i) {
            evaluators.push_back(
                interp3d_template.eval_proxy(eval_coord_3d[i]));
        }
        timer.pause();

        for (std::size_t t = 0; t < repeat_time; ++t) {
            timer.start("3D Mesh");
            for (std::size_t i = 0; i < data_len; ++i) {
                for (std::size_t j = 0; j < data_len; ++j) {
                    for (std::size_t k = 0; k < len_3d; ++k) {
                        mesh_3d(i, j, k) =
                            std::sin(static_cast<double>(i) * dt_3d - M_PI) *
                            std::cos(static_cast<double>(j) * dt_3d - M_PI +
                                     static_cast<double>(t) * delta) *
                            std::exp(-std::pow(
                                static_cast<double>(k) * dt_3d_aperiodic - .5,
                                2));
                    }
                }
            }

            timer.pause_and_start("3D Interpolate");
            auto interp3d = interp3d_template.interpolate(mesh_3d);

            timer.pause_and_start("3D Evaluation (Direct)");
            for (std::size_t i = 0; i < eval_count; ++i) {
                diff[i] = interp3d(eval_coord_3d[i]);
            }

            timer.pause_and_start("3D Evaluation (Use proxy)");
            for (std::size_t i = 0; i < eval_count; ++i) {
                diff[i] -= evaluators[i](interp3d);
            }
            timer.pause();
        }
        std::cout << "Interpolation on a 3D mesh consisting " << mesh_3d.size()
                  << " points. Then evaluate the function on " << eval_count
                  << " points. Repeat for " << repeat_time << " times.\n";
        double diff_L2{};
        for (auto v : diff) { diff_L2 += v * v; }
        diff_L2 = std::sqrt(diff_L2 / eval_count);
        std::cout << "Difference between two methods: " << diff_L2 << '\n';
        assertion(diff_L2 < std::numeric_limits<double>::epsilon());

        timer.print();
        timer.reset();
        std::cout << '\n';
    }
#endif  // INTP_CELL_LAYOUT
    return assertion.status();
}
