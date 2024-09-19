#include <Interpolation.hpp>
#include "include/Timer.h"

#include <algorithm>  // sort
#include <iomanip>    // setw
#include <iostream>
#include <random>
#include <vector>

// M_PI is not part of the standard
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

constexpr std::size_t max_len_power = 6 * 4;
constexpr std::size_t eval_count_power = 20;
constexpr std::size_t eval_count = 1 << eval_count_power;

template <std::size_t dim>
auto generate_coordinates() {
    std::mt19937 rand_gen(std::random_device{}());
    std::uniform_real_distribution<> uni(-M_PI, M_PI);

    std::vector<std::array<double, dim>> eval_coords;
    eval_coords.reserve(eval_count);
    for (std::size_t i = 0; i < eval_count; ++i) {
        eval_coords.push_back({});
        for (std::size_t j = 0; j < dim; ++j) {
            eval_coords.back()[j] = uni(rand_gen);
        }
    }
    return eval_coords;
}

auto sort_coord(auto coords) {
    constexpr std::size_t dim =
        std::tuple_size<typename decltype(coords)::value_type>::value;
    constexpr double dt_min = 2. * M_PI / (1 << (max_len_power / dim));

    std::sort(coords.begin(), coords.end(), [](auto c0, auto c1) {
        auto to_linear =
            []<std::size_t... idx>(std::index_sequence<idx...>, auto coord) {
            return (... +
                    (static_cast<std::size_t>((coord[idx] + M_PI) / dt_min)
                     << ((dim - idx - 1) * max_len_power / dim)));
        };
        auto i0 = to_linear(std::make_index_sequence<dim>{}, c0);
        auto i1 = to_linear(std::make_index_sequence<dim>{}, c1);
        return i0 < i1;
    });
    return coords;
}

int main() {
    using namespace intp;
    auto& timer = Timer::get_timer();

    // random sample points

    constexpr std::size_t min_len_power = 6 * 2;
    constexpr std::size_t max_len_power = 6 * 4;

    auto test_nDx = [&]<std::size_t... dim>(std::index_sequence<dim...>,
                                            auto order, auto point_nums) {
        constexpr auto dimension = sizeof...(dim);
        std::array<double, dimension> dt{(2 * M_PI / point_nums[dim])...};

        const static auto eval_coord = generate_coordinates<dimension>();
        const static auto eval_coord_sorted = sort_coord(eval_coord);

        auto& timer = Timer::get_timer();
        timer.reset();

        timer.start("Mesh");

        auto pn = point_nums;
#ifndef INTP_PERIODIC_NO_DUMMY_POINT
        for (std::size_t i = 0; i < mesh.size(); ++i) {
            pn[i] =
                pn[i] + ((i + 1) % 2);  // 0th, 2nd, ... dimension is periodic
        }
#endif
        Mesh<double, dimension> mesh(pn);
        for (std::size_t i = 0; i < mesh.size(); ++i) {
            auto index = mesh.dimension().dimwise_indices(i);
            mesh(index) =
                (... *
                 std::cos(static_cast<double>(index[dim]) * dt[dim] - M_PI));
        }

        constexpr auto interp_name = "Interpolation";
        timer.pause_and_start(interp_name);

        InterpolationFunction<double, dimension, order> interpND(
            {(dim % 2 == 0)...}, mesh, (dim, std::make_pair(-M_PI, M_PI))...);

        constexpr auto eval_rand_name = "Evaluation (Random)";
        timer.pause_and_start(eval_rand_name);

        double diff{};
        for (auto& x : eval_coord) { diff += interpND(x); }

        constexpr auto eval_seq_name = "Evaluation (Sequential)";
        timer.pause_and_start(eval_seq_name);

        for (auto& x : eval_coord_sorted) { diff -= interpND(x); }

        timer.pause();

        std::cout << dimension << "D mesh(" << mesh.size() << "), order "
                  << order.value << " complete\n";
#ifdef EMME_DEBUG
        std::cout << "Evaluate " << eval_coord.size()
                  << " times, unsorted and sorted. The diffreence is " << diff
                  << ". (Due to float point arithmetic error if it is not 0)\n";

        timer.print();
        std::cout << '\n';
#endif

        std::array result{timer.get_duration(interp_name),
                          timer.get_duration(eval_rand_name),
                          timer.get_duration(eval_seq_name)};

        return result;
    };

    using namespace std::chrono;

    using spline_orders =
        std::make_index_sequence<6>;  // interpolation orders vary from 0 to 5

    auto time_consumptions = ([&]<auto... dim>(std::index_sequence<dim...>) {
        return std::array{([&]<auto d>(std::integral_constant<std::size_t, d>) {
            auto func = [&]<auto... idx>(std::index_sequence<idx...> dim_seq,
                                         auto order, auto p) {
                // expand dimension-wise size
                return test_nDx(dim_seq, order,
                                std::array<std::size_t, sizeof...(idx)>{
                                    1u << ((p + idx) / sizeof...(idx))...});
            };
            return ([&]<auto... order>(std::index_sequence<order...>) {
                std::vector<
                    std::array<std::array<high_resolution_clock::duration, 3>,
                               sizeof...(order)>>
                    time_consumption;
                for (std::size_t p = min_len_power; p <= max_len_power;
                     ++p) {  // 2^p = point number
                    // expand order
                    time_consumption.push_back({func(
                        std::make_index_sequence<d>{},
                        std::integral_constant<std::size_t, order>{}, p)...});
                }
                return time_consumption;
            })(spline_orders{});
        })(std::integral_constant<std::size_t, dim>{})...};
    })(std::index_sequence<1, 2, 3>{} /* dimensions */);

    constexpr std::size_t col_width_1 = 8;

    auto print_table_head = [](std::size_t col_width_2) {
        const auto col_width_2s = spline_orders::size() * col_width_2;
        std::cout << std::left << '+';
        for (std::size_t i = 0; i < col_width_1 - 1; ++i) { std::cout << '-'; }
        std::cout << '+';
        for (std::size_t i = 0; i < col_width_2s - 1; ++i) { std::cout << '-'; }
        std::cout << "+\n";
        std::cout << "|mesh   |" << std::setw(col_width_2s - 1) << std::left
                  << "spline order"
                  << "|\n";
        std::cout << "|point  ";
        for (std::size_t i = 0; i < spline_orders::size(); ++i) {
            std::cout << '+';
            for (std::size_t j = 0; j < col_width_2 - 1; ++j) {
                std::cout.put('-');
            }
        }
        std::cout << "+\n|number |";
        for (std::size_t i = 0; i < spline_orders::size(); ++i) {
            std::cout << std::setw(col_width_2 - 1) << i << '|';
        }
        std::cout << "\n+";
        for (std::size_t i = 0; i < col_width_1 - 1; ++i) { std::cout << '-'; }
        for (std::size_t i = 0; i < spline_orders::size(); ++i) {
            std::cout << '+';
            for (std::size_t j = 0; j < col_width_2 - 1; ++j) {
                std::cout.put('-');
            }
        }
        std::cout << "+\n";
    };
    auto print_table_foot = [](std::size_t col_width_2) {
        std::cout << std::left << '+';
        for (std::size_t i = 0; i < col_width_1 - 1; ++i) { std::cout << '-'; }
        for (std::size_t i = 0; i < spline_orders::size(); ++i) {
            std::cout << '+';
            for (std::size_t j = 0; j < col_width_2 - 1; ++j) {
                std::cout.put('-');
            }
        }
        std::cout << "+\n";
    };

    const auto default_precision = std::cout.precision();
    std::cout << std::fixed << std::setprecision(2);
    std::size_t dim = 1;
    for (auto& dimensional_result : time_consumptions) {
        std::cout << "\n" << dim << "D Interpolation time consumption\n";
        // table head
        constexpr std::size_t col_width_interp_time = 11;
        print_table_head(col_width_interp_time);
        auto p = min_len_power;
        for (auto& row : dimensional_result) {
            std::cout << "|2^" << std::setw(col_width_1 - 3) << std::left
                      << p++;
            for (auto& ts : row) {
                std::cout
                    << "|" << std::setw(col_width_interp_time - 1) << std::right
                    << duration<double, milliseconds::period>(ts[0]).count();
            }
            std::cout << "|\n";
        }
        print_table_foot(col_width_interp_time);

        std::cout << "\n"
                  << dim << "D Evaluation(2^" << eval_count_power
                  << ", random and sequential) time consumption\n";
        constexpr std::size_t col_width_eval_time = 14;
        print_table_head(col_width_eval_time);
        p = min_len_power;
        for (auto& row : dimensional_result) {
            std::cout << "|2^" << std::setw(col_width_1 - 3) << std::left
                      << p++;
            for (auto& ts : row) {
                std::cout
                    << "|" << std::setw(6) << std::right
                    << duration<double, milliseconds::period>(ts[1]).count()
                    << ',' << std::setw(6) << std::right
                    << duration<double, milliseconds::period>(ts[2]).count();
            }
            std::cout << "|\n";
        }
        print_table_foot(col_width_eval_time);
        ++dim;
    }
    std::cout << std::setprecision(default_precision);

    return 0;
}
