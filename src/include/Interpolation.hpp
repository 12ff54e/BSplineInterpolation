#pragma once

#include <cmath>  // ceil
#include <initializer_list>

#include <Eigen/SparseLU>

#include "BSpline.hpp"

namespace intp {

template <typename T, size_t D>
class InterpolationFunctionTemplate;  // Forward declaration, since constructor
                                      // of InterpolationFunction will use it.

template <typename T, size_t D>
class InterpolationFunction {  // TODO: Add integration
   private:
    using val_type = T;
    using spline_type = BSpline<T, D>;
    using size_type = typename spline_type::size_type;
    using coord_type = typename spline_type::knot_type;

    const size_type order;
    const static size_type dim = D;

    spline_type __spline;

    template <typename _T>
    using DimArray = std::array<_T, dim>;

    DimArray<coord_type> __dx;
    DimArray<bool> __periodicity;
    DimArray<bool> __uniform;

    friend class InterpolationFunctionTemplate<T, D>;

    // auxiliary methods

    template <size_type... di>
    inline val_type call_op_helper(util::index_sequence<di...>,
                                   DimArray<coord_type> c) const {
        return __spline(std::make_pair(
            c[di], __uniform[di]
                       ? std::min(__spline.knots_num(di) - order - 2,
                                  (size_type)std::ceil(std::max(
                                      0., (c[di] - range(di).first) / __dx[di] -
                                              (__periodicity[di]
                                                   ? 1.
                                                   : .5 * (order + 1)))) +
                                      order)
                       : order)...);
    }

    template <size_type... di>
    inline val_type derivative_helper(util::index_sequence<di...>,
                                      DimArray<coord_type> c,
                                      DimArray<size_type> d) const {
        return __spline.derivative_at(std::make_tuple(
            (coord_type)c[di], (size_type)d[di],
            __uniform[di]
                ? std::min(
                      __spline.knots_num(di) - order - 2,
                      (size_type)std::ceil(std::max(
                          0.,
                          (c[di] - range(di).first) / __dx[di] -
                              (__periodicity[di] ? 1. : .5 * (order + 1)))) +
                          order)
                : order)...);
    }

    // overload for uniform knots
    template <typename _T>
    typename std::enable_if<std::is_arithmetic<_T>::value>::type
    __create_knot_vector(size_type dim_ind,
                         const MeshDimension<dim>& mesh_dimension,
                         DimArray<typename spline_type::KnotContainer>&,
                         std::pair<_T, _T> x_range) {
        __uniform[dim_ind] = true;
        const size_type n = mesh_dimension.dim_size(dim_ind);
        __dx[dim_ind] = (x_range.second - x_range.first) / (n - 1);

        const size_t extra =
            __periodicity[dim_ind] ? 2 * order + (1 - order % 2) : order + 1;

        std::vector<typename spline_type::knot_type> xs(n + extra,
                                                        x_range.first);

        if (__periodicity[dim_ind]) {
            for (size_type i = 0; i < xs.size(); ++i) {
                xs[i] = x_range.first + (i - .5 * extra) * __dx[dim_ind];
            }
        } else {
            for (size_type i = order + 1; i < xs.size() - order - 1; ++i) {
                xs[i] = x_range.first + (i - .5 * extra) * __dx[dim_ind];
            }
            for (size_type i = xs.size() - order - 1; i < xs.size(); ++i) {
                xs[i] = x_range.second;
            }
        }

        __spline.load_knots(dim_ind, std::move(xs), __periodicity[dim_ind]);
    }

    // overload for nonuniform knots, given by iterator pair
    template <typename _T>
    typename std::enable_if<std::is_convertible<
        typename std::iterator_traits<_T>::iterator_category,
        std::input_iterator_tag>::value>::type
    __create_knot_vector(
        size_type dim_ind,
        const MeshDimension<dim>& mesh_dimension,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        std::pair<_T, _T> x_range) {
        __uniform[dim_ind] = false;
        const size_type n = std::distance(x_range.first, x_range.second);
        if (n != mesh_dimension.dim_size(dim_ind)) {
            throw std::range_error(
                std::string("Inconsistency between knot number and "
                            "interpolated value number at dimension ") +
                std::to_string(dim_ind));
        }
        typename spline_type::KnotContainer xs(
            __periodicity[dim_ind] ? n + 2 * order + (1 - order % 2)
                                   : n + order + 1);

        input_coords[dim_ind].reserve(n);
        if (__periodicity[dim_ind]) {
            auto iter = x_range.first;

            input_coords[dim_ind].push_back(*iter);
            for (size_type i = order + 1; i < order + n; ++i) {
                val_type present = *(++iter);
                xs[i] = order % 2 == 0
                            ? .5 * (input_coords[dim_ind].back() + present)
                            : present;
                input_coords[dim_ind].push_back(present);
            }
            val_type period =
                input_coords[dim_ind].back() - input_coords[dim_ind].front();
            for (size_type i = 0; i < order + 1; ++i) {
                xs[i] = xs[n + i - 1] - period;
                xs[xs.size() - i - 1] = xs[xs.size() - i - n] + period;
            }
        } else {
            auto it = x_range.first;
            // Notice that *it++ is not guarantee to work as what you
            // expected for input iterators.
            auto l_knot = *it;
            // fill lestmost *order+1* identical knots
            for (size_type i = 0; i < order + 1; ++i) { xs[i] = l_knot; }
            // first knot is same as first input coordinate
            input_coords[dim_ind].emplace_back(l_knot);
            // Every knot in middle is average of *order* input
            // coordinates. This var is to track the sum of a moving window with
            // width *order*.
            coord_type window_sum{};
            for (size_type i = 1; i < order; ++i) {
                input_coords[dim_ind].emplace_back(*(++it));
                window_sum += input_coords[dim_ind][i];
            }
            for (size_type i = order + 1; i < n; ++i) {
                input_coords[dim_ind].emplace_back(*(++it));
                window_sum += input_coords[dim_ind][i - 1];
                xs[i] = window_sum / order;
                window_sum -= input_coords[dim_ind][i - order];
            }
            auto r_knot = *(++it);
            // fill rightmost *order+1* identical knots
            for (size_type i = n; i < n + order + 1; ++i) { xs[i] = r_knot; }
            // last knot is same as last input coordinate
            input_coords[dim_ind].emplace_back(r_knot);
        }
#ifdef _DEBUG
        std::cout << "[DEBUG] Nonuniform knots along dimension" << dim_ind
                  << ":\n";
        for (auto& c : xs) { std::cout << "[DEBUG] " << c << '\n'; }
        std::cout << std::endl;
#endif

        __spline.load_knots(dim_ind, std::move(xs));
    }

    template <typename... Ts, size_type... di>
    void __create_knots(
        util::index_sequence<di...>,
        MeshDimension<dim> mesh_dimension,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        std::pair<Ts, Ts>... x_ranges) {
#if __cplusplus >= 201703L
        (__create_knot_vector(di, f_mesh, input_coords, x_ranges), ...);
#else
        // polyfill of C++17 fold expression over comma
        std::array<std::nullptr_t, sizeof...(Ts)>{
            (__create_knot_vector(di, mesh_dimension, input_coords, x_ranges),
             nullptr)...};
#endif
    }

    void __solve_for_control_points(
        const Mesh<val_type, dim>& f_mesh,
        DimArray<typename spline_type::KnotContainer>& input_coords) {
        // initialize mesh storing weights of spline, adjust dimension according
        // to periodicity
        Mesh<val_type, dim> weights;
        {
            DimArray<size_type> dim_size_tmp;
            for (size_type d = 0; d < dim; ++d) {
                dim_size_tmp[d] =
                    f_mesh.dim_size(d) - (__periodicity[d] ? 1 : 0);
            }
            weights.resize(dim_size_tmp);
        }

        DimArray<typename spline_type::BaseSpline> base_spline_vals_per_dim;

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
            if (skip_flag) { continue; }
            weights(f_indices) = *it;
        }

        // pre-calculate base spline of periodic dimension, since it never
        // changes due to its even-spaced knots
        for (size_type d = 0; d < dim; ++d) {
            if (__periodicity[d] && __uniform[d]) {
                base_spline_vals_per_dim[d] = __spline.base_spline_value(
                    d, __spline.knots_begin(d) + order,
                    __spline.knots_begin(d)[order] +
                        (1 - order % 2) * __dx[d] * .5);
            }
        }

        // loop through each dimension to solve for control points
        for (size_type d = 0; d < dim; ++d) {
            std::vector<Eigen::Triplet<val_type>> coef_list;
            // rough estimate of upper limit of coefficient number, reserve
            // space to make sure no re-allocation occurs during the filling
            // process
            coef_list.reserve(weights.dim_size(d) * (order + 1));

            for (size_type i = 0; i < weights.dim_size(d); ++i) {
                const auto knot_num = __spline.knots_num(d);
                // This is the index of knot point to the left of i-th
                // interpolated value's coordinate, notice that knot points has
                // a larger gap in both ends in non-periodic case.
                size_type knot_ind{};
                if (__uniform[d]) {
                    knot_ind =
                        __periodicity[d]
                            ? i + order
                            : std::min(
                                  knot_num - order - 2,
                                  i > order / 2 ? i + (order + 1) / 2 : order);
                    if (!__periodicity[d]) {
                        if (knot_ind <= 2 * order + 1 ||
                            knot_ind >= knot_num - 2 * order - 2) {
                            // update base spline
                            const auto iter =
                                __spline.knots_begin(d) + knot_ind;
                            const coord_type x =
                                __spline.range(d).first + i * __dx[d];
                            base_spline_vals_per_dim[d] =
                                __spline.base_spline_value(d, iter, x);
                        }
                    }
                } else {
                    coord_type x = input_coords[d][i];
                    // using BSpline::get_knot_iter to find current
                    // knot_ind
                    const auto iter =
                        __periodicity[d] ? __spline.knots_begin(d) + i + order
                        : i == 0         ? __spline.knots_begin(d) + order
                        : i == input_coords[d].size() - 1
                            ? __spline.knots_end(d) - order - 2
                            : __spline.get_knot_iter(
                                  d, x, i + 1,
                                  std::min(knot_num - order - 1, i + order));
                    knot_ind = iter - __spline.knots_begin(d);
                    base_spline_vals_per_dim[d] =
                        __spline.base_spline_value(d, iter, x);
                }

                for (size_type j = 0; j < order + 1; ++j) {
                    coef_list.emplace_back(
                        i, (knot_ind - order + j) % weights.dim_size(d),
                        base_spline_vals_per_dim[d][j]);
                }
            }

            Eigen::SparseMatrix<double> coef(weights.dim_size(d),
                                             weights.dim_size(d));
            coef.setFromTriplets(coef_list.begin(), coef_list.end());
            Eigen::SparseLU<Eigen::SparseMatrix<double>,
                            Eigen::COLAMDOrdering<int>>
                solver(coef);

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
                Eigen::VectorXd weights_1d = solver.solve(mesh_val_1d);
                ind_1d = 0;
                for (auto it = weights.begin(d, ind_arr);
                     it != weights.end(d, ind_arr); ++it, ++ind_1d) {
                    *it = weights_1d(ind_1d);
                }
            }
        }

        __spline.load_ctrlPts(std::move(weights));
    }

    inline void __boundary_check(const DimArray<coord_type>& coord) const {
        for (size_type d = 0; d < dim; ++d) {
            if (!__periodicity[d] &&
                (coord[d] < range(d).first || coord[d] > range(d).second)) {
                throw std::domain_error(
                    "Given coordinate out of interpolation function range!");
            }
        }
    }

   public:
    /**
     * @brief Construct a new 1D Interpolation Function object, mimicking
     * Mathematica's `Interpolation` function, with option `Method->"Spline"`.
     *
     * @tparam InputIter
     * @param order order of interpolation, the interpolated function is of
     * $C^{order-1}$
     * @param periodic whether to construct a periodic spline
     * @param f_range a pair of iterators defining to-be-interpolated data
     * @param x_range a pair of x_min and x_max
     */
    template <
        typename InputIter,
        typename C1,
        typename C2,
        typename std::enable_if<
            dim == 1u &&
            std::is_convertible<
                typename std::iterator_traits<InputIter>::iterator_category,
                std::input_iterator_tag>::value>::type* = nullptr>
    InterpolationFunction(size_type order,
                          bool periodic,
                          std::pair<InputIter, InputIter> f_range,
                          std::pair<C1, C2> x_range)
        : InterpolationFunction(
              order,
              {periodic},
              Mesh<val_type, 1u>{std::make_pair(f_range.first, f_range.second)},
              static_cast<std::pair<typename std::common_type<C1, C2>::type,
                                    typename std::common_type<C1, C2>::type>>(
                  x_range)) {}

    template <
        typename InputIter,
        typename C1,
        typename C2,
        typename std::enable_if<
            dim == 1u &&
            std::is_convertible<
                typename std::iterator_traits<InputIter>::iterator_category,
                std::input_iterator_tag>::value>::type* = nullptr>
    InterpolationFunction(size_type order,
                          std::pair<InputIter, InputIter> f_range,
                          std::pair<C1, C2> x_range)
        : InterpolationFunction(order, false, f_range, x_range) {}

    /**
     * @brief Construct a new nD Interpolation Function object, mimicking
     * Mathematica's `Interpolation` function, with option `Method->"Spline"`.
     * Notice: last value of periodic dimension will be discarded since it is
     * considered same of the first value. Thus inconsistency input data for
     * periodic interpolation will be accepted.
     *
     * @param order order of interpolation, the interpolated function is of
     * $C^{order-1}$
     * @param periodicity an array describing periodicity of each dimension
     * @param f_mesh a mesh containing data to be interpolated
     * @param x_ranges pairs of x_min and x_max or begin and end iterator
     */
    template <typename... Ts>
    InterpolationFunction(size_type order,
                          DimArray<bool> periodicity,
                          const Mesh<val_type, dim>& f_mesh,
                          std::pair<Ts, Ts>... x_ranges)
        : InterpolationFunction(InterpolationFunctionTemplate<val_type, dim>{
              order, periodicity, f_mesh.dimension(), x_ranges...}
                                    .interpolate(f_mesh)) {}

    // Non-periodic for all dimension
    template <typename... Ts>
    InterpolationFunction(size_type order,
                          const Mesh<val_type, dim>& f_mesh,
                          std::pair<Ts, Ts>... x_ranges)
        : InterpolationFunction(order, {}, f_mesh, x_ranges...) {}

    // constructor for partial construction, that is, without interpolated
    // values
    template <typename... Ts>
    InterpolationFunction(
        size_type order,
        DimArray<bool> periodicity,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        MeshDimension<dim> mesh_dimension,
        std::pair<Ts, Ts>... x_ranges)
        : order(order),
          __spline(periodicity, order),
          __periodicity(periodicity) {
        // load knots into spline
        __create_knots(util::make_index_sequence_for<Ts...>{},
                       std::move(mesh_dimension), input_coords, x_ranges...);
    }

    /**
     * @brief Get spline value.
     *
     * @param x coordinates
     */
    template <typename... Coords,
              typename Indices = util::make_index_sequence_for<Coords...>,
              typename = typename std::enable_if<std::is_arithmetic<
                  typename std::common_type<Coords...>::type>::value>::type>
    val_type operator()(Coords... x) const {
        return call_op_helper(Indices{}, DimArray<coord_type>{x...});
    }

    /**
     * @brief Get spline value.
     *
     * @param coord coordinate array
     */
    val_type operator()(DimArray<coord_type> coord) const {
        return call_op_helper(util::make_index_sequence<dim>{}, coord);
    }

    /**
     * @brief Get spline value, but with out of boundary check.
     *
     * @param coord coordinate array
     */
    val_type at(DimArray<coord_type> coord) const {
        __boundary_check(coord);
        return call_op_helper(util::make_index_sequence<dim>{}, coord);
    }

    /**
     * @brief Get spline value, but with out of boundary check.
     *
     * @param x coordinates
     */
    template <typename... Coords,
              typename Indices = util::make_index_sequence_for<Coords...>,
              typename = typename std::enable_if<std::is_arithmetic<
                  typename std::common_type<Coords...>::type>::value>::type>
    val_type at(Coords... x) const {
        return at(DimArray<coord_type>{static_cast<coord_type>(x)...});
    }

    /**
     * @brief Get spline derivative value.
     *
     * @param coord coordinate array
     * @param derivatives derivative order array
     */
    val_type derivative(DimArray<coord_type> coord,
                        DimArray<size_type> derivatives) const {
        return derivative_helper(util::make_index_sequence<dim>{}, coord,
                                 derivatives);
    }

    /**
     * @brief Get spline derivative value.
     *
     * @param coord coordinate array
     * @param deriOrder derivative orders
     */
    template <typename... Args>
    val_type derivative(DimArray<coord_type> coord, Args... deriOrder) const {
        return derivative(coord, DimArray<size_type>{(size_type)deriOrder...});
    }

    /**
     * @brief Get spline derivative value.
     *
     * @param coord_deriOrder_pair pairs of coordinate and derivative order
     */
    template <typename... CoordDeriOrderPair>
    val_type derivative(CoordDeriOrderPair... coord_deriOrder_pair) const {
        return derivative(
            DimArray<coord_type>{coord_deriOrder_pair.first...},
            DimArray<size_type>{(size_type)coord_deriOrder_pair.second...});
    }

    /**
     * @brief Get spline derivative value, but with out of boundary check.
     *
     * @param coord coordinate array
     * @param derivatives derivative order array
     */
    val_type derivative_at(DimArray<coord_type> coord,
                           DimArray<size_type> derivatives) const {
        __boundary_check(coord);
        return derivative_helper(util::make_index_sequence<dim>{}, coord,
                                 derivatives);
    }

    /**
     * @brief Get spline derivative value, but with out of boundary check.
     *
     * @param coord coordinate array
     * @param deriOrder derivative orders
     */
    template <typename... Args>
    val_type derivative_at(DimArray<coord_type> coord,
                           Args... deriOrder) const {
        return derivative_at(coord,
                             DimArray<size_type>{(size_type)deriOrder...});
    }

    /**
     * @brief Get spline derivative value, but with out of boundary check.
     *
     * @param coord_deriOrder_pair pairs of coordinate and derivative order
     */
    template <typename... CoordDeriOrderPair>
    val_type derivative_at(CoordDeriOrderPair... coord_deriOrder_pair) const {
        return derivative_at(
            DimArray<coord_type>{coord_deriOrder_pair.first...},
            DimArray<size_type>{(size_type)coord_deriOrder_pair.second...});
    }

    // properties

    bool periodicity(size_type dim_ind) const { return __periodicity[dim_ind]; }

    bool uniform(size_type dim_ind) const { return __uniform[dim_ind]; }

    const std::pair<typename spline_type::knot_type,
                    typename spline_type::knot_type>&
    range(size_type dim_ind) const {
        return __spline.range(dim_ind);
    }

    /**
     * @brief Get a ref of underlying spline object
     *
     * @return spline_type&
     */
    const spline_type& spline() const { return __spline; }
};

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

}  // namespace intp
