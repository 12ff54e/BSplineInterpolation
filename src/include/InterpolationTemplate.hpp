#ifndef INTP_TEMPLATE
#define INTP_TEMPLATE

#include "BSpline.hpp"
#include "BandLU.hpp"
#include "DedicatedThreadPool.hpp"
#include "Mesh.hpp"
#include "util.hpp"

#if __cplusplus >= 201703L
#include <variant>
#endif

#ifdef _TRACE
#include <iostream>
#endif

namespace intp {

template <typename T, size_t D, typename U = double>
class InterpolationFunction;  // Forward declaration, since template has
                              // a member of it.

/**
 * @brief Template for interpolation with only coordinates, and generate
 * interpolation function when fed by function values.
 *
 */
template <typename T, size_t D, typename U = double>
class InterpolationFunctionTemplate {
   public:
    using function_type = InterpolationFunction<T, D, U>;
    using size_type = typename function_type::size_type;
    using coord_type = typename function_type::coord_type;
    using val_type = typename function_type::val_type;
    using diff_type = typename function_type::diff_type;

    using ctrl_pt_type =
        typename function_type::spline_type::ControlPointContainer;

    static constexpr size_type dim = D;

    template <typename V>
    using DimArray = std::array<V, dim>;

    using MeshDim = MeshDimension<dim>;

    /**
     * @brief Construct a new Interpolation Function Template object, all other
     * constructors delegate to this one
     *
     * @param order Order of BSpline
     * @param periodicity Periodicity of each dimension
     * @param interp_mesh_dimension The structure of coordinate mesh
     * @param x_ranges Begin and end iterator/value pairs of each dimension
     */
    template <typename... Ts>
    InterpolationFunctionTemplate(size_type order,
                                  DimArray<bool> periodicity,
                                  MeshDim interp_mesh_dimension,
                                  std::pair<Ts, Ts>... x_ranges)
        : input_coords_{},
          mesh_dimension_(interp_mesh_dimension),
          base_(order,
                periodicity,
                input_coords_,
                mesh_dimension_,
                x_ranges...),
          solvers_{} {
        // active union member accordingly
        for (size_type i = 0; i < dim; ++i) {
#if __cplusplus >= 201703L
            if (periodicity[i]) {
                solvers_[i].template emplace<extended_solver_type>();
            } else {
                solvers_[i].template emplace<base_solver_type>();
            }
#endif
        }
        build_solver_();
    }

    /**
     * @brief Construct a new 1D Interpolation Function Template object.
     *
     * @param order order of interpolation, the interpolated function is of
     * $C^{order-1}$
     * @param periodicity whether to construct a periodic spline
     * @param f_length point number of to-be-interpolated data
     * @param x_range a pair of x_min and x_max
     */
    template <typename C1, typename C2>
    InterpolationFunctionTemplate(size_type order,
                                  bool periodicity,
                                  size_type f_length,
                                  std::pair<C1, C2> x_range)
        : InterpolationFunctionTemplate(order,
                                        {periodicity},
                                        MeshDim{f_length},
                                        x_range) {
        static_assert(
            dim == size_type{1},
            "You can only use this overload of constructor in 1D case.");
    }

    template <typename C1, typename C2>
    InterpolationFunctionTemplate(size_type order,
                                  size_type f_length,
                                  std::pair<C1, C2> x_range)
        : InterpolationFunctionTemplate(order, false, f_length, x_range) {}

    /**
     * @brief Construct a new (aperiodic) Interpolation Function Template
     * object
     *
     * @param order Order of BSpline
     * @param interp_mesh_dimension The structure of coordinate mesh
     * @param x_ranges Begin and end iterator/value pairs of each dimension
     */
    template <typename... Ts>
    InterpolationFunctionTemplate(size_type order,
                                  MeshDim interp_mesh_dimension,
                                  std::pair<Ts, Ts>... x_ranges)
        : InterpolationFunctionTemplate(order,
                                        {},
                                        interp_mesh_dimension,
                                        x_ranges...) {}

    template <typename MeshOrIterPair>
    function_type interpolate(MeshOrIterPair&& mesh_or_iter_pair) const& {
        function_type interp{base_};
        interp.spline_.load_ctrlPts(
            solve_for_control_points_(Mesh<val_type, dim>{
                std::forward<MeshOrIterPair>(mesh_or_iter_pair)}));
        return interp;
    }

    template <typename MeshOrIterPair>
    function_type interpolate(MeshOrIterPair&& mesh_or_iter_pair) && {
        base_.spline_.load_ctrlPts(
            solve_for_control_points_(Mesh<val_type, dim>{
                std::forward<MeshOrIterPair>(mesh_or_iter_pair)}));
        return std::move(base_);
    }

    // modify the given interpolation function
    template <typename MeshOrIterPair>
    void interpolate(function_type& interp,
                     MeshOrIterPair&& mesh_or_iter_pair) const& {
        interp = base_;
        interp.spline_.load_ctrlPts(
            solve_for_control_points_(Mesh<val_type, dim>{
                std::forward<MeshOrIterPair>(mesh_or_iter_pair)}));
    }

   private:
    using base_solver_type = BandLU<BandMatrix<coord_type>>;
    using extended_solver_type = BandLU<ExtendedBandMatrix<coord_type>>;

    // input coordinates, needed only in nonuniform case
    DimArray<typename function_type::spline_type::KnotContainer> input_coords_;

    MeshDim mesh_dimension_;

    // the base interpolation function with unspecified weights
    function_type base_;

#if __cplusplus >= 201703L
    using EitherSolver = std::variant<base_solver_type, extended_solver_type>;
#else
    // A union-like class storing BandLU solver for Either band matrix or
    // extended band matrix
    struct EitherSolver {
        // set default active union member, or g++ compiled code will throw
        // err
        EitherSolver() {}
        // Active union member and tag it.
        EitherSolver(bool is_periodic)
            : is_active_(true), is_periodic_(is_periodic) {
            if (is_periodic_) {
                new (&solver_periodic) extended_solver_type;
            } else {
                new (&solver_aperiodic) base_solver_type;
            }
        }
        EitherSolver& operator=(bool is_periodic) {
            if (is_active_) {
                throw std::runtime_error("Can not switch solver type.");
            }

            is_active_ = true;
            is_periodic_ = is_periodic;
            if (is_periodic_) {
                new (&solver_periodic) extended_solver_type;
            } else {
                new (&solver_aperiodic) base_solver_type;
            }

            return *this;
        }

        // Copy constructor is required by aggregate initialization, but never
        // invoked in this code since the content of this union is never
        // switched.
        EitherSolver(const EitherSolver&) {}
        // Destructor needed and it invokes either member's destructor according
        // to the tag "is_periodic";
        ~EitherSolver() {
            if (is_active_) {
                if (is_periodic_) {
                    solver_periodic.~extended_solver_type();
                } else {
                    solver_aperiodic.~base_solver_type();
                }
            }
        }

        union {
            base_solver_type solver_aperiodic;
            extended_solver_type solver_periodic;
        };

        bool is_active_ = false;
        // tag for union member
        bool is_periodic_ = false;
    };
#endif

    // solver for weights
    DimArray<EitherSolver> solvers_;

    void build_solver_() {
        const auto& order = base_.order_;
        // adjust dimension according to periodicity
        {
            DimArray<size_type> dim_size_tmp;
            for (size_type d = 0; d < dim; ++d) {
                dim_size_tmp[d] = mesh_dimension_.dim_size(d) -
                                  (base_.periodicity(d) ? 1 : 0);
            }
            mesh_dimension_.resize(dim_size_tmp);
        }

        DimArray<typename function_type::spline_type::BaseSpline>
            base_spline_vals_per_dim;

        const auto& spline = base_.spline();

        // pre-calculate base spline of periodic dimension, since it never
        // changes due to its even-spaced knots
        for (size_type d = 0; d < dim; ++d) {
            if (base_.periodicity(d) && base_.uniform(d)) {
                base_spline_vals_per_dim[d] = spline.base_spline_value(
                    d, spline.knots_begin(d) + static_cast<diff_type>(order),
                    spline.knots_begin(d)[static_cast<diff_type>(order)] +
                        (1 - order % 2) * base_.dx_[d] * coord_type{.5});
            }
        }

#ifdef _TRACE
        std::cout << "\n[TRACE] Coefficient Matrices\n";
#endif

        // loop through each dimension to construct coefficient matrix
        for (size_type d = 0; d < dim; ++d) {
            bool periodic = base_.periodicity(d);
            bool uniform = base_.uniform(d);
            auto mat_dim = mesh_dimension_.dim_size(d);
            auto band_width = periodic ? order / 2 : order - 1;

#if __cplusplus >= 201703L
            std::variant<typename base_solver_type::matrix_type,
                         typename extended_solver_type::matrix_type>
                coef_mat;
            if (periodic) {
                coef_mat.template emplace<
                    typename extended_solver_type::matrix_type>(
                    mat_dim, band_width, band_width);
            } else {
                coef_mat
                    .template emplace<typename base_solver_type::matrix_type>(
                        mat_dim, band_width, band_width);
            }
#else
            typename extended_solver_type::matrix_type coef_mat(
                mat_dim, band_width, band_width);
#endif

#ifdef _TRACE
            std::cout << "\n[TRACE] Dimension " << d << '\n';
            std::cout << "[TRACE] {0, 0} -> 1\n";
#endif

            for (size_type i = 0; i < mesh_dimension_.dim_size(d); ++i) {
                if (!periodic) {
                    // In aperiodic case, first and last data point can only
                    // be covered by one base spline, and the base spline at
                    // these ending points eval to 1.
                    if (i == 0 || i == mesh_dimension_.dim_size(d) - 1) {
#if __cplusplus >= 201703L
                        std::visit([i](auto& m) { m(i, i) = 1; }, coef_mat);
#else
                        coef_mat.main_bands_val(i, i) = 1;
#endif
                        continue;
                    }
                }

                const auto knot_num = spline.knots_num(d);
                // This is the index of knot point to the left of i-th
                // interpolated value's coordinate, notice that knot points has
                // a larger gap in both ends in non-periodic case.
                size_type knot_ind{};
                // flag for internal points in uniform aperiodic case
                const bool is_internal =
                    i > order / 2 &&
                    i < mesh_dimension_.dim_size(d) - order / 2 - 1;

                if (uniform) {
                    knot_ind =
                        periodic ? i + order
                                 : std::min(knot_num - order - 2,
                                            i > order / 2 ? i + (order + 1) / 2
                                                          : order);
                    if (!periodic) {
                        if (knot_ind <= 2 * order + 1 ||
                            knot_ind >= knot_num - 2 * order - 2) {
                            // out of the zone of even-spaced knots, update base
                            // spline
                            const auto iter = spline.knots_begin(d) +
                                              static_cast<diff_type>(knot_ind);
                            const coord_type x =
                                spline.range(d).first +
                                static_cast<coord_type>(i) * base_.dx_[d];
                            base_spline_vals_per_dim[d] =
                                spline.base_spline_value(d, iter, x);
                        }
                    }
                } else {
                    coord_type x = input_coords_[d][i];
                    // using BSpline::get_knot_iter to find current
                    // knot_ind
                    const auto iter =
                        periodic ? spline.knots_begin(d) +
                                       static_cast<diff_type>(i + order)
                        : i == 0 ? spline.knots_begin(d) +
                                       static_cast<diff_type>(order)
                        : i == input_coords_[d].size() - 1
                            ? spline.knots_end(d) -
                                  static_cast<diff_type>(order + 2)
                            : spline.get_knot_iter(
                                  d, x, i + 1,
                                  std::min(knot_num - order - 1, i + order));
                    knot_ind =
                        static_cast<size_type>(iter - spline.knots_begin(d));
                    base_spline_vals_per_dim[d] =
                        spline.base_spline_value(d, iter, x);
                }

                // number of base spline that covers present data point.
                const size_type s_num = periodic                 ? order | 1
                                        : order == 1             ? 1
                                        : uniform && is_internal ? order | 1
                                                                 : order + 1;
                for (size_type j = 0; j < s_num; ++j) {
                    const size_type row = (i + (periodic ? band_width : 0)) %
                                          mesh_dimension_.dim_size(d);
                    const size_type col =
                        (knot_ind - order + j) % mesh_dimension_.dim_size(d);
#if __cplusplus >= 201703L
                    std::visit(
                        [&](auto& m) {
                            m(row, col) = base_spline_vals_per_dim[d][j];
                        },
                        coef_mat);
#else
                    if (periodic) {
                        coef_mat((i + band_width) % mesh_dimension_.dim_size(d),
                                 (knot_ind - order + j) %
                                     mesh_dimension_.dim_size(d)) =
                            base_spline_vals_per_dim[d][j];
                    } else {
                        coef_mat.main_bands_val(i, knot_ind - order + j) =
                            base_spline_vals_per_dim[d][j];
                    }
#endif
#ifdef _TRACE
                    std::cout
                        << "[TRACE] {"
                        << (periodic
                                ? (i + band_width) % mesh_dimension_.dim_size(d)
                                : i)
                        << ", "
                        << (knot_ind - order + j) % mesh_dimension_.dim_size(d)
                        << "} -> " << base_spline_vals_per_dim[d][j] << '\n';
#endif
                }
            }

#ifdef _TRACE
            std::cout << "[TRACE] {" << mesh_dimension_.dim_size(d) - 1 << ", "
                      << mesh_dimension_.dim_size(d) - 1 << "} -> 1\n";
#endif

#if __cplusplus >= 201703L
            std::visit(
                [&](auto& solver) {
                    using matrix_t = typename util::remove_cvref_t<
                        decltype(solver)>::matrix_type;
                    solver.compute(std::get<matrix_t>(std::move(coef_mat)));
                },
                solvers_[d]);
#else
            solvers_[d] = periodic;
            if (periodic) {
                solvers_[d].solver_periodic.compute(std::move(coef_mat));
            } else {
                solvers_[d].solver_aperiodic.compute(
                    static_cast<typename base_solver_type::matrix_type&&>(
                        coef_mat));
            }
#endif
        }
    }

    ctrl_pt_type solve_for_control_points_(
        const Mesh<val_type, dim>& f_mesh) const {
        ctrl_pt_type weights{mesh_dimension_};
        ctrl_pt_type weights_tmp(1);  // auxilary weight for swapping between
        if CPP17_CONSTEXPR_ (dim > 1) { weights_tmp.resize(mesh_dimension_); }

        auto check_idx =
            [&](typename Mesh<val_type, dim>::index_type& indices) {
                bool keep_flag = true;
                for (size_type d = 0; d < dim; ++d) {
                    if (base_.periodicity(d)) {
                        // Skip last point of periodic dimension
                        keep_flag =
                            keep_flag && indices[d] != weights.dim_size(d);
                        indices[d] = (indices[d] + weights.dim_size(d) +
                                      base_.order_ / 2) %
                                     weights.dim_size(d);
                    }
                }
                return keep_flag;
            };

        // Copy interpolating values into weights mesh as the initial state of
        // the iterative control points solving algorithm
        for (auto it = f_mesh.begin(); it != f_mesh.end(); ++it) {
            auto f_indices = f_mesh.iter_indices(it);
            if (check_idx(f_indices)) { weights(f_indices) = *it; }
        }

        auto array_right_shift = [](DimArray<size_type> arr) {
            DimArray<size_type> new_arr{};
            for (size_type d_ = 0; d_ < dim; ++d_) {
                new_arr[d_] = arr[(d_ + dim - 1) % dim];
            }
            return new_arr;
        };

        // loop through each dimension to solve for control points
        for (size_type d = 0; d < dim; ++d) {
            ctrl_pt_type& old_weight = d % 2 == 0 ? weights : weights_tmp;
            ctrl_pt_type& new_weight = d % 2 != 0 ? weights : weights_tmp;

            if CPP17_CONSTEXPR_ (dim > 1) {
                new_weight.resize(array_right_shift(old_weight.dimension()));
            }

            const auto line_size = old_weight.dim_size(dim - 1);
            // size of hyperplane orthogonal to last dim axis
            const auto hyperplane_size = old_weight.size() / line_size;

            // prepare variables being captured by lambda
            auto& solver_wrapper = solvers_[dim - 1 - d];
            auto solve_and_rearrange_block = [&](size_type begin,
                                                 size_type end) {
                for (size_type j = begin; j < end; ++j) {
                    auto ind_arr =
                        old_weight.dimension().dimwise_indices(j * line_size);

                    auto old_iter_begin = old_weight.begin(dim - 1, ind_arr);
                    auto old_iter_end = old_weight.end(dim - 1, ind_arr);
                    auto new_iter_begin =
                        new_weight.begin(0, array_right_shift(ind_arr));
#if __cplusplus >= 201703L
                    std::visit(
                        [&](auto& solver) { solver.solve(old_iter_begin); },
                        solver_wrapper);
#else
                    if (base_.periodicity(dim - 1 - d)) {
                        solver_wrapper.solver_periodic.solve(old_iter_begin);
                    } else {
                        solver_wrapper.solver_aperiodic.solve(old_iter_begin);
                    }
#endif
                    if CPP17_CONSTEXPR_ (dim > 1) {
                        for (auto old_it = old_iter_begin,
                                  new_it = new_iter_begin;
                             old_it != old_iter_end; ++old_it, ++new_it) {
                            *new_it = *old_it;
                        }
                    }
                }
            };

#ifdef _MULTITHREAD
            // TODO: use a more robust task division strategy
            const size_type block_num = static_cast<size_type>(
                std::sqrt(static_cast<double>(hyperplane_size)));
            const size_type task_per_block =
                block_num == 0 ? hyperplane_size : hyperplane_size / block_num;

            auto& thread_pool = DedicatedThreadPool<void>::get_instance(8);
            std::vector<std::future<void>> res;

            for (size_type i = 0; i < block_num; ++i) {
                // use thread pool here may seem to be a bit overhead
                res.push_back(thread_pool.queue_task([=]() {
                    solve_and_rearrange_block(i * task_per_block,
                                              (i + 1) * task_per_block);
                }));
            }
            // main thread deals with the remaining part in case hyperplane_size
            // not divisible by thread_num
            solve_and_rearrange_block(block_num * task_per_block,
                                      hyperplane_size);
            // wait for all tasks are complete
            for (auto&& f : res) { f.get(); }
#else
            solve_and_rearrange_block(0, hyperplane_size);
#endif  // _MULTITHREAD
        }

        if CPP17_CONSTEXPR_ (dim % 2 == 0 || dim == 1) {
            return weights;
        } else {
            return weights_tmp;
        }
    }
};

template <typename T = double, typename U = double>
class InterpolationFunctionTemplate1D
    : public InterpolationFunctionTemplate<T, size_t{1}, U> {
   private:
    using base = InterpolationFunctionTemplate<T, size_t{1}, U>;

   public:
    InterpolationFunctionTemplate1D(typename base::size_type f_length,
                                    typename base::size_type order = 3,
                                    bool periodicity = false)
        : InterpolationFunctionTemplate1D(
              std::make_pair(
                  typename base::coord_type{},
                  static_cast<typename base::coord_type>(f_length - 1)),
              f_length,
              order,
              periodicity) {}

    template <typename C1, typename C2>
    InterpolationFunctionTemplate1D(std::pair<C1, C2> x_range,
                                    typename base::size_type f_length,
                                    typename base::size_type order = 3,
                                    bool periodicity = false)
        : base(order, periodicity, f_length, x_range) {}
};

}  // namespace intp

#endif
