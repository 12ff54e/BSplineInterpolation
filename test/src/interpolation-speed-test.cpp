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

    std::mt19937 rand_gen(static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<> rand_dist(-M_PI, M_PI);
    std::uniform_real_distribution<> rand_dist2(-.5, .5);

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
    std::vector<double> eval_coord_1d_sorted(eval_coord_1d);
    std::vector<std::array<double, 2>> eval_coord_2d_sorted(eval_coord_2d);
    std::vector<std::array<double, 3>> eval_coord_3d_sorted(eval_coord_3d);

    std::sort(eval_coord_1d_sorted.begin(), eval_coord_1d_sorted.end());

    constexpr size_t len_2d_max = 1 << (len_power / 2);
    constexpr size_t len_3d_max = 1 << (len_power / 3);
    constexpr double dt_min = 2 * M_PI / static_cast<double>(len_2d_max);
    constexpr double dt_3d_min = 2 * M_PI / len_3d_max;
    constexpr double dt_3d_aperiodic_min = 1. / (len_3d_max - 1);

    std::sort(
        eval_coord_2d_sorted.begin(), eval_coord_2d_sorted.end(),
        [](const std::array<double, 2>& p1, const std::array<double, 2>& p2) {
            const auto x_1 =
                static_cast<int>(std::floor((p1[0] + M_PI) / dt_min));
            const auto x_2 =
                static_cast<int>(std::floor((p2[0] + M_PI) / dt_min));
            return x_1 < x_2 ||
                   (x_1 == x_2 && std::floor((p1[1] + M_PI) / dt_min) <
                                      std::floor((p2[1] + M_PI) / dt_min));
        });

    std::sort(
        eval_coord_3d_sorted.begin(), eval_coord_3d_sorted.end(),
        [](const std::array<double, 3>& p1, const std::array<double, 3>& p2) {
            const auto x_1 =
                static_cast<int>(std::floor((p1[0] + M_PI) / dt_3d_min));
            const auto x_2 =
                static_cast<int>(std::floor((p2[0] + M_PI) / dt_3d_min));
            const auto y_1 =
                static_cast<int>(std::floor((p1[1] + M_PI) / dt_3d_min));
            const auto y_2 =
                static_cast<int>(std::floor((p2[1] + M_PI) / dt_3d_min));
            return x_1 < x_2 || (x_1 == x_2 && y_1 < y_2) ||
                   (x_1 == x_2 && y_1 == y_2 &&
                    std::floor((p1[2] + .5) / dt_3d_aperiodic_min) <
                        std::floor((p2[2] + .5) / dt_3d_aperiodic_min));
        });

    // 1D case
    auto test_1D = [&eval_coord_1d, &eval_coord_1d_sorted](auto order) {
        constexpr size_t len_1d = 1 << len_power;
        constexpr double dx = 2 * M_PI / (len_1d);
        std::vector<double> vec_1d{};

        auto& timer = Timer::get_timer();
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
        timer.pause_and_start("1D Interpolation");

        InterpolationFunction1D<decltype(order)::value> interp1d(
            std::make_pair(-M_PI, M_PI),
            std::make_pair(vec_1d.begin(), vec_1d.end()), true);

        timer.pause_and_start("1D Evaluation (Random)");

        double diff{};  // prevent loops being optimized out
        for (auto& x : eval_coord_1d) { diff += interp1d(x); }

        timer.pause_and_start("1D Evaluation (Sequential)");

        for (auto& x : eval_coord_1d_sorted) { diff -= interp1d(x); }

        timer.pause();

        std::cout << "1D mesh(" << vec_1d.size() << "), order " << order.value
                  << ", evaluate " << eval_count
                  << " times, unsorted and sorted.\nThe diffreence is " << diff
                  << ". (Due to float point arithmetic error if it is not 0)\n";

        timer.print();
        timer.reset();

        std::cout << '\n';
    };

    // 2d case
    auto test_2D = [&eval_coord_2d, &eval_coord_2d_sorted](auto order) {
        constexpr size_t len_2d = 1 << (len_power / 2);
        constexpr double dt = 2 * M_PI / static_cast<double>(len_2d);

        auto& timer = Timer::get_timer();
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

        timer.pause_and_start("2D Interpolation");

        InterpolationFunction<double, 2, decltype(order)::value> interp2d(
            {true, true}, trig_mesh_2d, std::make_pair(-M_PI, M_PI),
            std::make_pair(-M_PI, M_PI));

        timer.pause_and_start("2D Evaluation (Random)");

        double diff{};
        for (auto& x : eval_coord_2d) { diff += interp2d(x); }

        timer.pause_and_start("2D Evaluation (Sequential)");

        for (auto& x : eval_coord_2d_sorted) { diff -= interp2d(x); }

        timer.pause();

        std::cout << "Interpolation on a 2D mesh consisting "
                  << trig_mesh_2d.size()
                  << " points. Then evaluate the function on " << eval_count
                  << " points, unsorted and sorted.\nThe diffreence is " << diff
                  << ". (Due to float point arithmetic error if it is not 0)\n";

        timer.print();
        timer.reset();

        std::cout << '\n';
    };

    // 3d case
    auto test_3D = [&eval_coord_3d, &eval_coord_3d_sorted](auto order) {
        constexpr size_t len_3d = 1 << (len_power / 3);
        constexpr double dt_3d = 2 * M_PI / len_3d;
        constexpr double dt_3d_aperiodic = 1. / (len_3d - 1);

        auto& timer = Timer::get_timer();
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

        timer.pause_and_start("3D Interpolation");

        InterpolationFunction<double, 3, decltype(order)::value> interp3d(
            {true, true, false}, mesh_3d, std::make_pair(-M_PI, M_PI),
            std::make_pair(-M_PI, M_PI), std::make_pair(-.5, .5));

        timer.pause_and_start("3D Evaluation (Random)");

        double diff{};
        for (auto& x : eval_coord_3d) { diff += interp3d(x); }

        timer.pause_and_start("3D Evaluation (Sequential)");

        for (auto& x : eval_coord_3d_sorted) { diff -= interp3d(x); }

        timer.pause();

        std::cout << "Interpolation on a 3D mesh consisting " << mesh_3d.size()
                  << " points. Then evaluate the function on " << eval_count
                  << " points, unsorted and sorted.\nThe diffreence is " << diff
                  << ". (Due to float point arithmetic error if it is not 0)\n";

        timer.print();
        timer.reset();

        std::cout << '\n';
    };

    ([&]<auto... order>(std::index_sequence<order...>) {
        (test_1D(std::integral_constant<std::size_t, order>{}), ...);
    })(std::index_sequence<3, 4, 5>{});

    ([&]<auto... order>(std::index_sequence<order...>) {
        (test_2D(std::integral_constant<std::size_t, order>{}), ...);
    })(std::index_sequence<3, 4, 5>{});

    ([&]<auto... order>(std::index_sequence<order...>) {
        (test_3D(std::integral_constant<std::size_t, order>{}), ...);
    })(std::index_sequence<3, 4, 5>{});

    return assertion.status();
}
