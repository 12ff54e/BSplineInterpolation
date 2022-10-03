#pragma once

#include "BSpline.hpp"
#include "BandLU.hpp"
#include "Mesh.hpp"

#ifdef _TRACE
#include <iostream>
#endif

namespace intp {

template <typename T, size_t D>
class InterpolationFunction;  // Forward declaration, since template has
                              // a member of it.

/**
 * @brief Template for interpolation with only coordinates, and generate
 * interpolation function when fed by function values.
 *
 */
template <typename T, size_t D>
class InterpolationFunctionTemplate {
   public:
    using function_type = InterpolationFunction<T, D>;
    using size_type = typename function_type::size_type;
    using coord_type = typename function_type::coord_type;
    using val_type = typename function_type::val_type;

    static constexpr size_type dim = D;

    template <typename U>
    using DimArray = std::array<U, dim>;

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
        : input_coords{},
          mesh_dimension(interp_mesh_dimension),
          base(order, periodicity, input_coords, mesh_dimension, x_ranges...),
          solvers{} {
        // active union member accordingly
        for (size_type i = 0; i < dim; ++i) { solvers[i] = periodicity[i]; }
        __build_solver();
    }

    /**
     * @brief Construct a new 1D Interpolation Function Template object.
     *
     * @param order order of interpolation, the interpolated function is of
     * $C^{order-1}$
     * @param periodic whether to construct a periodic spline
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
        function_type interp{base};
        interp.__spline.load_ctrlPts(
            __solve_for_control_points(Mesh<val_type, dim>{
                std::forward<MeshOrIterPair>(mesh_or_iter_pair)}));
        return interp;
    }

    template <typename MeshOrIterPair>
    function_type interpolate(MeshOrIterPair&& mesh_or_iter_pair) && {
        base.__spline.load_ctrlPts(
            __solve_for_control_points(Mesh<val_type, dim>{
                std::forward<MeshOrIterPair>(mesh_or_iter_pair)}));
        return std::move(base);
    }

   private:
    // input coordinates, needed only in nonuniform case
    DimArray<typename function_type::spline_type::KnotContainer> input_coords;

    MeshDim mesh_dimension;

    // the base interpolation function with unspecified weights
    function_type base;

    // A union-like class storing BandLU solver for Either band matrix or
    // extended band matrix
    struct EitherSolver {
        EitherSolver() : is_active_(false) {}
        // Active union member and tag it.
        EitherSolver(bool is_periodic)
            : is_active_(true), is_periodic_(is_periodic) {
            if (is_periodic_) {
                solver_periodic = {};
            } else {
                solver_aperiodic = {};
            }
        }
        EitherSolver& operator=(bool is_periodic) {
            if (is_active_) {
                throw std::runtime_error("Can not switch solver type.");
            }

            is_active_ = true;
            is_periodic_ = is_periodic;
            if (is_periodic_) {
                solver_periodic = {};
            } else {
                solver_aperiodic = {};
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
                    solver_periodic.~BandLU<ExtendedBandMatrix<val_type>>();
                } else {
                    solver_aperiodic.~BandLU<BandMatrix<val_type>>();
                }
            }
        }

        union {
            BandLU<BandMatrix<val_type>> solver_aperiodic;
            BandLU<ExtendedBandMatrix<val_type>> solver_periodic;
        };

        bool is_active_;
        // tag for union member
        bool is_periodic_;
    };

    // solver for weights
    DimArray<EitherSolver> solvers;

    void __build_solver() {
        const auto& order = base.order;
        // adjust dimension according to periodicity
        {
            DimArray<size_type> dim_size_tmp;
            for (size_type d = 0; d < dim; ++d) {
                dim_size_tmp[d] =
                    mesh_dimension.dim_size(d) - (base.periodicity(d) ? 1 : 0);
            }
            mesh_dimension.resize(dim_size_tmp);
        }

        DimArray<typename function_type::spline_type::BaseSpline>
            base_spline_vals_per_dim;

        const auto& spline = base.spline();

        // pre-calculate base spline of periodic dimension, since it never
        // changes due to its even-spaced knots
        for (size_type d = 0; d < dim; ++d) {
            if (base.periodicity(d) && base.uniform(d)) {
                base_spline_vals_per_dim[d] = spline.base_spline_value(
                    d, spline.knots_begin(d) + order,
                    spline.knots_begin(d)[order] +
                        (1 - order % 2) * base.__dx[d] * .5);
            }
        }

#ifdef _TRACE
        std::cout << "\n[TRACE] Coefficient Matrices\n";
#endif

        // loop through each dimension to construct coefficient matrix
        for (size_type d = 0; d < dim; ++d) {
            bool periodic = base.periodicity(d);
            bool uniform = base.uniform(d);
            auto mat_dim = mesh_dimension.dim_size(d);
            auto band_width = periodic ? order / 2 : order - 1;
            ExtendedBandMatrix<val_type> coef_mat{mat_dim, band_width,
                                                  band_width};

#ifdef _TRACE
            std::cout << "\n[TRACE] Dimension " << d << '\n';
            std::cout << "[TRACE] {0, 0} -> 1\n";
#endif

            for (size_type i = 0; i < mesh_dimension.dim_size(d); ++i) {
                if (!periodic) {
                    // In aperiodic case, first and last data point can only
                    // be covered by one base spline, and the base spline at
                    // these ending points eval to 1.
                    if (i == 0 || i == mesh_dimension.dim_size(d) - 1) {
                        coef_mat.main_bands_val(i, i) = 1;
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
                    i < mesh_dimension.dim_size(d) - order / 2 - 1;

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
                            const auto iter = spline.knots_begin(d) + knot_ind;
                            const coord_type x =
                                spline.range(d).first + i * base.__dx[d];
                            base_spline_vals_per_dim[d] =
                                spline.base_spline_value(d, iter, x);
                        }
                    }
                } else {
                    coord_type x = input_coords[d][i];
                    // using BSpline::get_knot_iter to find current
                    // knot_ind
                    const auto iter =
                        periodic ? spline.knots_begin(d) + i + order
                        : i == 0 ? spline.knots_begin(d) + order
                        : i == input_coords[d].size() - 1
                            ? spline.knots_end(d) - order - 2
                            : spline.get_knot_iter(
                                  d, x, i + 1,
                                  std::min(knot_num - order - 1, i + order));
                    knot_ind = iter - spline.knots_begin(d);
                    base_spline_vals_per_dim[d] =
                        spline.base_spline_value(d, iter, x);
                }

                // number of base spline that covers present data point.
                const size_type s_num = periodic                 ? order | 1
                                        : order == 1             ? 1
                                        : uniform && is_internal ? order | 1
                                                                 : order + 1;
                for (size_type j = 0; j < s_num; ++j) {
                    if (periodic) {
                        coef_mat((i + band_width) % mesh_dimension.dim_size(d),
                                 (knot_ind - order + j) %
                                     mesh_dimension.dim_size(d)) =
                            base_spline_vals_per_dim[d][j];
                    } else {
                        coef_mat.main_bands_val(i, knot_ind - order + j) =
                            base_spline_vals_per_dim[d][j];
                    }
#ifdef _TRACE
                    std::cout
                        << "[TRACE] {"
                        << (periodic
                                ? (i + band_width) % mesh_dimension.dim_size(d)
                                : i)
                        << ", "
                        << (knot_ind - order + j) % mesh_dimension.dim_size(d)
                        << "} -> " << base_spline_vals_per_dim[d][j] << '\n';
#endif
                }
            }

#ifdef _TRACE
            std::cout << "[TRACE] {" << mesh_dimension.dim_size(d) - 1 << ", "
                      << mesh_dimension.dim_size(d) - 1 << "} -> 1\n";
#endif

            if (periodic) {
                solvers[d].solver_periodic.compute(coef_mat);
            } else {
                solvers[d].solver_aperiodic.compute(
                    static_cast<BandMatrix<val_type>>(coef_mat));
            }
        }
    }

    Mesh<val_type, dim> __solve_for_control_points(
        const Mesh<val_type, dim>& f_mesh) const {
        Mesh<val_type, dim> weights{mesh_dimension};

        auto check_idx =
            [&](typename Mesh<val_type, dim>::index_type& indices) {
                bool keep_flag = true;
                for (size_type d = 0; d < dim; ++d) {
                    if (base.periodicity(d)) {
                        // Skip last point of periodic dimension
                        keep_flag = indices[d] != weights.dim_size(d);
                        indices[d] = (indices[d] + weights.dim_size(d) +
                                      base.order / 2) %
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

        // loop through each dimension to solve for control points
        for (size_type d = 0; d < dim; ++d) {
            // size of hyperplane when given dimension is fixed
            size_type hyperplane_size = weights.size() / weights.dim_size(d);

            // loop over each point (representing a 1D spline) of hyperplane
            for (size_type i = 0; i < hyperplane_size; ++i) {
                DimArray<size_type> ind_arr{};
                for (size_type d_ = 0, total_ind = i; d_ < dim; ++d_) {
                    if (d_ == d) { continue; }
                    ind_arr[d_] = total_ind % weights.dim_size(d_);
                    total_ind /= weights.dim_size(d_);
                }

                // Loop through one dimension, update interpolating value to
                // control points.
                // In periodic case, rows are shifted to make coefficient matrix
                // diagonal dominate so weights column should be shifted
                // accordingly.
                // auto iter = weights.begin(d, ind_arr);
                if (base.periodicity(d)) {
                    solvers[d].solver_periodic.solve(weights.begin(d, ind_arr));
                } else {
                    solvers[d].solver_aperiodic.solve(
                        weights.begin(d, ind_arr));
                }
            }
        }

        return weights;
    }
};

template <typename T = double>
class InterpolationFunctionTemplate1D
    : public InterpolationFunctionTemplate<T, size_t{1}> {
   private:
    using base = InterpolationFunctionTemplate<T, size_t{1}>;

   public:
    InterpolationFunctionTemplate1D(typename base::size_type f_length,
                                    typename base::size_type order = 3,
                                    bool periodicity = false)
        : InterpolationFunctionTemplate1D(
              std::make_pair(
                  (typename base::coord_type){},
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
