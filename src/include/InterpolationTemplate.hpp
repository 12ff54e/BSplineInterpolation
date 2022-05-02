#pragma once

#include "BSpline.hpp"
#include "Mesh.hpp"

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
     * @brief Construct a new Interpolation Function Template object
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
          base(order, periodicity, input_coords, mesh_dimension, x_ranges...) {
        // adjust dimension according to periodicity
        {
            DimArray<size_type> dim_size_tmp;
            bool p_flag = false;
            for (size_type d = 0; d < dim; ++d) {
                dim_size_tmp[d] =
                    mesh_dimension.dim_size(d) -
                    (base.periodicity(d) ? ((p_flag = true), 1) : 0);
            }
            if (p_flag) { mesh_dimension.resize(dim_size_tmp); }
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

        // loop through each dimension to construct coefficient matrix
        for (size_type d = 0; d < dim; ++d) {
            std::vector<Eigen::Triplet<val_type>> coef_list;
            // rough estimate of upper limit of coefficient number, reserve
            // space to make sure no re-allocation occurs during the filling
            // process
            coef_list.reserve(mesh_dimension.dim_size(d) * (order + 1));

            for (size_type i = 0; i < mesh_dimension.dim_size(d); ++i) {
                const auto knot_num = spline.knots_num(d);
                // This is the index of knot point to the left of i-th
                // interpolated value's coordinate, notice that knot points has
                // a larger gap in both ends in non-periodic case.
                size_type knot_ind{};

                if (base.uniform(d)) {
                    knot_ind =
                        base.periodicity(d)
                            ? i + order
                            : std::min(
                                  knot_num - order - 2,
                                  i > order / 2 ? i + (order + 1) / 2 : order);
                    if (!base.periodicity(d)) {
                        if (knot_ind <= 2 * order + 1 ||
                            knot_ind >= knot_num - 2 * order - 2) {
                            // update base spline
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
                        base.periodicity(d) ? spline.knots_begin(d) + i + order
                        : i == 0            ? spline.knots_begin(d) + order
                        : i == input_coords[d].size() - 1
                            ? spline.knots_end(d) - order - 2
                            : spline.get_knot_iter(
                                  d, x, i + 1,
                                  std::min(knot_num - order - 1, i + order));
                    knot_ind = iter - spline.knots_begin(d);
                    base_spline_vals_per_dim[d] =
                        spline.base_spline_value(d, iter, x);
                }

                for (size_type j = 0; j < order + 1; ++j) {
                    coef_list.emplace_back(
                        i, (knot_ind - order + j) % mesh_dimension.dim_size(d),
                        base_spline_vals_per_dim[d][j]);
                }
            }

            // fill coefficient matrix
            {
                Eigen::SparseMatrix<double> coef(mesh_dimension.dim_size(d),
                                                 mesh_dimension.dim_size(d));
                coef.setFromTriplets(coef_list.begin(), coef_list.end());
                solver[d].compute(coef);
            }

#ifdef _DEBUG
            if (solver[d].info() != Eigen::Success) {
                throw std::runtime_error(
                    std::string{
                        "Coefficient matrix decomposition failed at dim "} +
                    std::to_string(d) + std::string{".\n"});
            }
#endif
        }
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
     * @brief Construct a new (nonperiodic) Interpolation Function Template
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

    // solver for weights
    DimArray<Eigen::SparseLU<Eigen::SparseMatrix<double>,
                             Eigen::COLAMDOrdering<int>>>
        solver;

    Mesh<val_type, dim> __solve_for_control_points(
        const Mesh<val_type, dim>& f_mesh) const {
        Mesh<val_type, dim> weights{mesh_dimension};

        // Copy interpolating values into weights mesh as the initial state of
        // the iterative control points solving algorithm
        for (auto it = f_mesh.begin(); it != f_mesh.end(); ++it) {
            auto f_indices = f_mesh.iter_indices(it);
            bool skip_flag = false;
            for (size_type d = 0; d < dim; ++d) {
                // Skip last point of periodic dimension
                if (f_indices[d] == weights.dim_size(d)) {
                    skip_flag = true;
                    break;
                }
            }
            if (!skip_flag) { weights(f_indices) = *it; }
        }

        // loop through each dimension to solve for control points
        for (size_type d = 0; d < dim; ++d) {
            // size of hyperplane when given dimension is fixed
            size_type hyperplane_size = weights.size() / weights.dim_size(d);

            // loop over each point (representing a 1D spline) of hyperplane
            for (size_type i = 0; i < hyperplane_size; ++i) {
                DimArray<size_type> ind_arr;
                for (size_type d_ = 0, total_ind = i; d_ < dim; ++d_) {
                    if (d_ == d) { continue; }
                    ind_arr[d_] = total_ind % weights.dim_size(d_);
                    total_ind /= weights.dim_size(d_);
                }

                // loop through one dimension, update interpolating value to
                // control points
                Eigen::VectorXd mesh_val_1d(weights.dim_size(d));
                size_type ind_1d{};
                for (auto it = weights.begin(d, ind_arr);
                     it != weights.end(d, ind_arr); ++it, ++ind_1d) {
                    mesh_val_1d(ind_1d) = *it;
                }
                Eigen::VectorXd weights_1d = solver[d].solve(mesh_val_1d);
                ind_1d = 0;
                for (auto it = weights.begin(d, ind_arr);
                     it != weights.end(d, ind_arr); ++it, ++ind_1d) {
                    *it = weights_1d(ind_1d);
                }
            }
        }

        return weights;
    }
};

template <typename T>
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
