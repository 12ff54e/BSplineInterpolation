#include <cmath>  // ceil
#include <initializer_list>

#include <Eigen/SparseLU>

#include "BSpline.hpp"

namespace intp {

template <typename T, unsigned D>
class InterpolationFunction {
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

    // auxiliary methods

    template <typename>
    using __cast_to_bool = bool;

    template <unsigned... di>
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

    template <unsigned... di>
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
                         const Mesh<val_type, dim>& f_mesh,
                         DimArray<typename spline_type::KnotContainer>&,
                         std::pair<_T, _T> x_range) {
        __uniform[dim_ind] = true;
        const size_type n = f_mesh.dim_size(dim_ind);
        __dx[dim_ind] = (x_range.second - x_range.first) / (n - 1);

        std::vector<typename spline_type::knot_type> xs(
            __periodicity[dim_ind] ? n + 2 * order : n + order + 1,
            x_range.first);

        if (__periodicity[dim_ind]) {
            for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
                xs[i] = x_range.first +
                        (i - static_cast<int>(order)) * __dx[dim_ind];
            }
        } else {
            for (size_type i = order + 1; i < xs.size() - order - 1; ++i) {
                xs[i] =
                    x_range.first + .5 * (2 * i - order - 1) * __dx[dim_ind];
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
        const Mesh<val_type, dim>& f_mesh,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        std::pair<_T, _T> x_range) {
        __uniform[dim_ind] = false;
        const unsigned n = std::distance(x_range.first, x_range.second);
        if (n != f_mesh.dim_size(dim_ind)) {
            throw std::range_error(
                std::string("Inconsistency between knot number and "
                            "interpolated value number at dimension ") +
                std::to_string(dim_ind));
        }
        typename spline_type::KnotContainer xs(
            __periodicity[dim_ind] ? n + 2 * order : n + order + 1);

        if (__periodicity[dim_ind]) {
            std::copy(x_range.first, x_range.second, xs.begin() + order);
            for (unsigned i = 0; i < order; ++i) {
                xs[i] = xs[order] - xs[n + order - 1] + xs[n + i - 1];
                xs[n + order + i] =
                    xs[n + order - 1] + xs[order + i + 1] - xs[order];
            }
        } else {
            input_coords[dim_ind].reserve(n);
            auto it = x_range.first;
            // Notice that *it++ is not guarantee to work as what you
            // expected for input iterators.
            auto l_knot = *it;
            // fill lestmost *order+1* identical knots
            for (unsigned i = 0; i < order + 1; ++i) { xs[i] = l_knot; }
            // first knot is same as first input coordinate
            input_coords[dim_ind].emplace_back(l_knot);
            // Every knot in middle is average of *order* input
            // coordinates. This var is to track the sum of a moving window with
            // width *order*.
            coord_type window_sum{};
            for (unsigned i = 1; i < order; ++i) {
                input_coords[dim_ind].emplace_back(*(++it));
                window_sum += input_coords[dim_ind][i];
            }
            for (unsigned i = order + 1; i < n; ++i) {
                input_coords[dim_ind].emplace_back(*(++it));
                window_sum += input_coords[dim_ind][i - 1];
                xs[i] = window_sum / order;
                window_sum -= input_coords[dim_ind][i - order];
            }
            auto r_knot = *(++it);
            // fill rightmost *order+1* identical knots
            for (unsigned i = n; i < n + order + 1; ++i) { xs[i] = r_knot; }
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
        const Mesh<val_type, dim>& f_mesh,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        std::pair<Ts, Ts>... x_ranges) {
        // polyfill of C++17 fold expression over comma
        std::array<std::nullptr_t, sizeof...(Ts)>{
            (__create_knot_vector(di, f_mesh, input_coords, x_ranges),
             nullptr)...};
    }

    void __solve_for_control_points(
        const Mesh<val_type, dim>& f_mesh,
        DimArray<typename spline_type::KnotContainer>& input_coords) {
        // initialize mesh storing weights of spline, adjust dimension according
        // to periodicity
        Mesh<val_type, dim> weights{f_mesh};
        {
            DimArray<size_type> dim_size_tmp;
            bool p_flag = false;
            for (size_type d = 0; d < dim; ++d) {
                dim_size_tmp[d] = weights.dim_size(d) -
                                  (__periodicity[d] ? ((p_flag = true), 1) : 0);
            }
            if (p_flag) { weights.resize(dim_size_tmp); }
        }

        // estimate numbers of base spline covering, accurate for uniform case
        // and nonuniform case with even order
        size_type entry_count{};
        for (size_type total_ind = 0; total_ind < f_mesh.size(); ++total_ind) {
            auto indices = f_mesh.dimwise_indices(total_ind);
            size_type c = 1;
            for (size_type d = 0; d < dim; ++d) {
                auto ind = indices[d];
                c *= __periodicity[d]                            ? order
                     : ind == 0 || ind == f_mesh.dim_size(d) - 1 ? 1
                     : ind <= order / 2 ||
                             ind >= f_mesh.dim_size(d) - 1 - order / 2
                         ? order + 1
                         : order | 1;
            }
            entry_count += c;
        }
        Eigen::VectorXd mesh_val(weights.size());

        DimArray<typename spline_type::BaseSpline> base_spline_vals_per_dim;

        // pre-calculate base spline of periodic dimension, since it never
        // changes due to its even-spaced knots
        for (size_type d = 0; d < dim; ++d) {
            if (__periodicity[d] && __uniform[d]) {
                base_spline_vals_per_dim[d] = __spline.base_spline_value(
                    d, __spline.knots_begin(d) + order, 0.);
            }
        }

        std::vector<Eigen::Triplet<val_type>> coef_list;
        coef_list.reserve(entry_count);

        size_type actual_ind = 0;
        // loop over dimension to calculate 1D base spline for each dim
        for (auto it = f_mesh.begin(); it != f_mesh.end(); ++it) {
            auto f_indices = f_mesh.iter_indices(it);

            // first base spline index (also the index of weights) in each
            // dimension
            decltype(f_indices) base_spline_anchor;

            // This flag indicates the current point is redundent due to
            // periodicity
            bool skip_flag = false;
            for (size_type d = 0; d < dim; ++d) {
                if (__periodicity[d] &&
                    f_indices[d] == f_mesh.dim_size(d) - 1) {
                    skip_flag = true;
                    break;
                }

                const auto knot_num = __spline.knots_num(d);
                // This is the index of i-th dimension knot vector to the left
                // of current f_indices[i] position, notice that knot points has
                // a larger gap in both ends in non-periodic case.
                size_type knot_ind{};
                if (__uniform[d]) {
                    knot_ind =
                        __periodicity[d]
                            ? f_indices[d] + order
                            : std::min(knot_num - order - 2,
                                       f_indices[d] > order / 2
                                           ? f_indices[d] + (order + 1) / 2
                                           : order);
                    if (!__periodicity[d]) {
                        if (knot_ind <= 2 * order + 1 ||
                            knot_ind >= knot_num - 2 * order - 2) {
                            // update base spline
                            const auto iter =
                                __spline.knots_begin(d) + knot_ind;
                            const coord_type x = __spline.range(d).first +
                                                 f_indices[d] * __dx[d];
                            base_spline_vals_per_dim[d] =
                                __spline.base_spline_value(d, iter, x);
                        }
                    }
                } else {
                    const coord_type x =
                        __periodicity[d]
                            ? __spline.knots_begin(d)[f_indices[d] + order]
                            : input_coords[d][f_indices[d]];
                    // using BSpline::get_knot_iter to find current
                    // knot_ind
                    const auto iter =
                        __periodicity[d]
                            ? __spline.knots_begin(d) + f_indices[d] + order
                        : f_indices[d] == 0 ? __spline.knots_begin(d) + order
                        : f_indices[d] == input_coords[d].size() - 1
                            ? __spline.knots_end(d) - order - 2
                            : __spline.get_knot_iter(
                                  d, x, f_indices[d] + 1,
                                  std::min(knot_num - order - 1,
                                           f_indices[d] + order));
                    knot_ind = iter - __spline.knots_begin(d);
                    base_spline_vals_per_dim[d] =
                        __spline.base_spline_value(d, iter, x);
                }

                base_spline_anchor[d] = knot_ind - order;
            }
            if (skip_flag) { continue; }

            // loop over nD base splines that contributes to current f_mesh
            // point, fill matrix
            for (size_type i = 0; i < util::pow(order + 1, dim); ++i) {
                DimArray<size_type> ind_arr;
                val_type spline_val = 1;
                for (int d = dim - 1, local_ind = i; d >= 0; --d) {
                    ind_arr[d] = local_ind % (order + 1);
                    local_ind /= (order + 1);

                    spline_val *= base_spline_vals_per_dim[d][ind_arr[d]];
                    ind_arr[d] += base_spline_anchor[d];

                    if (__periodicity[d]) {
                        ind_arr[d] %= (__spline.knots_num(d) - 2 * order - 1);
                    }
                }
                if (spline_val != val_type{0}) {
#ifdef _TRACE
                    std::cout << "[TRACE] {" << actual_ind << ','
                              << weights.indexing(ind_arr) << "} -> "
                              << spline_val << '\n';
#endif

                    coef_list.emplace_back(
                        actual_ind, weights.indexing(ind_arr), spline_val);
                }
            }

            // fill rhs vector
            mesh_val(actual_ind++) = *it;
        }

        // fill coefficient matrix
        Eigen::SparseMatrix<double> coef(weights.size(), weights.size());
        coef.setFromTriplets(coef_list.begin(), coef_list.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
            solver;
        solver.compute(coef);

        Eigen::VectorXd weights_eigen_vec;
        weights_eigen_vec = solver.solve(mesh_val);

#if EIGEN_VERSION_AT_LEAST(3, 4, 0)
        weights.storage = std::vector<val_type>(weights_eigen_vec.begin(),
                                                weights_eigen_vec.end());
#else
        std::vector<val_type> tmp;
        tmp.reserve(weights_eigen_vec.size());
        for (unsigned i = 0; i < weights_eigen_vec.size(); ++i) {
            tmp.emplace_back(weights_eigen_vec[i]);
        }
        weights.storage = std::move(tmp);
#endif

        __spline.load_ctrlPts(std::move(weights));
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
              Mesh<val_type, 1u>{f_range.first, f_range.second},
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
        : order(order),
          __spline(periodicity, order),
          __periodicity(periodicity) {
        // a temperory storage for input coordinates, only needed in nonuniform
        // and non-periodic case
        DimArray<typename spline_type::KnotContainer> input_coords;
        // load knots into spline
        __create_knots(util::make_index_sequence_for<Ts...>{}, f_mesh,
                       input_coords, x_ranges...);
        // solve control points
        __solve_for_control_points(f_mesh, input_coords);
    }

    // Non-periodic for all dimension
    template <typename... Ts>
    InterpolationFunction(size_type order,
                          const Mesh<val_type, dim>& f_mesh,
                          std::pair<Ts, Ts>... x_ranges)
        : InterpolationFunction(order, {}, f_mesh, x_ranges...) {}

    const std::pair<typename spline_type::knot_type,
                    typename spline_type::knot_type>&
    range(size_type dim_ind) const {
        return __spline.range(dim_ind);
    }

    template <typename... Coords,
              typename Indices = util::make_index_sequence_for<Coords...>,
              typename = typename std::enable_if<std::is_arithmetic<
                  typename std::common_type<Coords...>::type>::value>::type>
    val_type operator()(Coords... x) const {
        return call_op_helper(Indices{}, DimArray<coord_type>{x...});
    }

    val_type operator()(DimArray<coord_type> coord) const {
        return call_op_helper(util::make_index_sequence<dim>{}, coord);
    }

    template <typename... Args>
    val_type derivative_at(DimArray<coord_type> coord,
                           Args... deriOrder) const {
        return derivative_helper(util::make_index_sequence<dim>{}, coord,
                                 DimArray<size_type>{(size_type)deriOrder...});
    }

    val_type derivative_at(DimArray<coord_type> coord,
                           DimArray<size_type> derivatives) const {
        return derivative_helper(util::make_index_sequence<dim>{}, coord,
                                 derivatives);
    }

    template <typename... CoordDeriOrderPair>
    val_type derivative_at(CoordDeriOrderPair... coord_deriOrder_pair) const {
        return derivative_helper(
            util::make_index_sequence<dim>{},
            DimArray<coord_type>{coord_deriOrder_pair.first...},
            DimArray<size_type>{(size_type)coord_deriOrder_pair.second...});
    }

    bool periodicity(size_type dim_ind) const { return __periodicity[dim_ind]; }

    bool uniform(size_type dim_ind) const { return __uniform[dim_ind]; }

    /**
     * @brief Get a ref of underlying spline object
     *
     * @return spline_type&
     */
    const spline_type& spline() const { return __spline; }
};

}  // namespace intp
