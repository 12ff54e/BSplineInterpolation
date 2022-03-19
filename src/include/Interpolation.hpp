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

    bool __uniform;
    spline_type __spline;

    template <typename _T>
    using DimArray = std::array<_T, dim>;

    DimArray<coord_type> __dx;
    DimArray<bool> __periodicity;

    // auxiliary methods

    template <unsigned... di>
    inline val_type call_op_helper(util::index_sequence<di...>,
                                   DimArray<coord_type> c) const {
        return __uniform
                   ? __spline(std::make_pair(
                         c[di],
                         std::min(__spline.knots_num(di) - order - 2,
                                  (size_type)std::ceil(std::max(
                                      0., (c[di] - range(di).first) / __dx[di] -
                                              (__periodicity[di]
                                                   ? 1.
                                                   : .5 * (order + 1)))) +
                                      order))...)
                   : __spline((coord_type)c[di]...);
    }

    template <unsigned... di>
    inline val_type derivative_helper(util::index_sequence<di...>,
                                      DimArray<coord_type> c,
                                      DimArray<size_type> d) const {
        return __uniform
                   ? __spline.derivative_at(std::make_tuple(
                         (coord_type)c[di], (size_type)d[di],
                         std::min(__spline.knots_num(di) - order - 2,
                                  (size_type)std::ceil(std::max(
                                      0., (c[di] - range(di).first) / __dx[di] -
                                              (__periodicity[di]
                                                   ? 1.
                                                   : .5 * (order + 1)))) +
                                      order))...)
                   : __spline.derivative_at(std::make_pair(c[di], d[di])...);
    }

   public:
    /**
     * @brief Construct a new 1D Interpolation Function object, mimicking
     * Mathematica's `Interpolation` function, with option `Method->"Spline"`.
     *
     * @tparam InputIter
     * @param order order of interpolation, the interpolated function is of
     * $C^{order-1}$
     * @param f_range a pair of iterators defining to-be-interpolated data
     * @param x_range a pair of x_min and x_max
     */
    template <
        typename InputIter,
        typename = typename std::enable_if<
            dim == 1u && std::is_convertible<typename std::iterator_traits<
                                                 InputIter>::iterator_category,
                                             std::input_iterator_tag>::value,
            int>::type>
    InterpolationFunction(size_type order,
                          bool periodic,
                          std::pair<InputIter, InputIter> f_range,
                          std::pair<double, double> x_range)
        : InterpolationFunction(
              order,
              {periodic},
              Mesh<val_type, 1u>{f_range.first, f_range.second},
              x_range) {}

    template <
        typename InputIter,
        typename = typename std::enable_if<
            dim == 1u && std::is_convertible<typename std::iterator_traits<
                                                 InputIter>::iterator_category,
                                             std::input_iterator_tag>::value,
            int>::type>
    InterpolationFunction(size_type order,
                          std::pair<InputIter, InputIter> f_range,
                          std::pair<double, double> x_range)
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
     * @param x_ranges pairs of x_min and x_max
     */
    template <typename... Ts>
    InterpolationFunction(size_type order,
                          DimArray<bool> periodicity,
                          Mesh<val_type, dim> f_mesh,
                          std::pair<Ts, Ts>... x_ranges)
        : order(order),
          __uniform(true),
          __spline(periodicity, order),
          __periodicity(periodicity) {
        // load knots into spline
        util::dispatch_indexed(
            [&](size_type dim_ind,
                std::pair<typename spline_type::knot_type,
                          typename spline_type::knot_type> x_range) {
                const size_type n = f_mesh.dim_size(dim_ind);
                __dx[dim_ind] = (x_range.second - x_range.first) / (n - 1);

                std::vector<typename spline_type::knot_type> xs(
                    __periodicity[dim_ind] ? n + 2 * order : n + order + 1,
                    x_range.first);

                if (__periodicity[dim_ind]) {
                    for (int i = 0; i < xs.size(); ++i) {
                        xs[i] =
                            x_range.first + (i - (int)order) * __dx[dim_ind];
                    }
                } else {
                    for (int i = order + 1; i < xs.size() - order - 1; ++i) {
                        xs[i] = x_range.first +
                                .5 * (2 * i - order - 1) * __dx[dim_ind];
                    }
                    for (int i = xs.size() - order - 1; i < xs.size(); ++i) {
                        xs[i] = x_range.second;
                    }
                }

                __spline.load_knots(dim_ind, std::move(xs));
            },
            x_ranges...);

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

        // numbers of base spline covering
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
            if (__periodicity[d]) {
                base_spline_vals_per_dim[d] = __spline.base_spline_value(
                    d, __spline.knots_begin(d) + order, 0.);
            }
        }

        std::vector<Eigen::Triplet<val_type>> coef_list;
        coef_list.reserve(entry_count);

        size_type actual_ind = 0;
        // loop over dimension to calculate 1D base spline for each dim
        for (auto it = f_mesh.begin(); it != f_mesh.end(); ++it) {
            // size_type f_total_ind = std::distance(f_mesh.begin(), it);
            auto f_indices = f_mesh.iter_indices(it);

            // first base spline index (also the index of weights) in each
            // dimension
            decltype(f_indices) base_spline_anchor;

            // This flag indicates the current point is redundent due to
            // periodicity
            bool skip_flag = false;
            for (int d = 0; d < dim; ++d) {
                const auto knot_num = __spline.knots_num(d);
                // This is the index of i-th dimension knot vector to the left
                // of current f_indices[i] position, notice that knot points has
                // a larger gap in both ends in non-periodic case.
                const auto knot_ind =
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
                        const auto iter = __spline.knots_begin(d) + knot_ind;
                        const double x =
                            __spline.range(d).first + f_indices[d] * __dx[d];
                        base_spline_vals_per_dim[d] =
                            __spline.base_spline_value(d, iter, x);
                    }
                } else {
                    if (f_indices[d] == f_mesh.dim_size(d) - 1) {
                        skip_flag = true;
                        break;
                    }
                }
                base_spline_anchor[d] = knot_ind - order;
            }
            if (skip_flag) { continue; }

            // loop over nD base splines that contributes to current f_mesh
            // point, fill matrix
            for (int i = 0; i < util::pow(order + 1, dim); ++i) {
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
#ifdef _DEBUG
                    std::cout << "[DEBUG] {" << actual_ind << ','
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

        weights.storage = std::vector<val_type>(weights_eigen_vec.begin(),
                                                weights_eigen_vec.end());

        __spline.load_ctrlPts(std::move(weights));
    }

    template <typename... Ts>
    InterpolationFunction(size_type order,
                          Mesh<val_type, dim> f_mesh,
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
        return call_op_helper(util::make_index_sequence<dim>{},
                              std::move(coord));
    }

    template <typename... Args>
    val_type derivative_at(DimArray<coord_type> coord,
                           Args... deriOrder) const {
        return derivative_helper(util::make_index_sequence<dim>{},
                                 std::move(coord),
                                 DimArray<size_type>{deriOrder...});
    }

    val_type derivative_at(DimArray<coord_type> coord,
                           DimArray<size_type> derivatives) const {
        return derivative_helper(util::make_index_sequence<dim>{},
                                 std::move(coord), std::move(derivatives));
    }

    template <typename... CoordDeriOrderPair>
    val_type derivative_at(CoordDeriOrderPair... coord_deriOrder_pair) const {
        return derivative_helper(
            util::make_index_sequence<dim>{},
            DimArray<coord_type>{coord_deriOrder_pair.first...},
            DimArray<size_type>{coord_deriOrder_pair.second...});
    }

    bool periodicity(size_type dim_ind) const { return __periodicity[dim_ind]; }

    /**
     * @brief Get a ref of underlying spline object
     *
     * @return spline_type&
     */
    const spline_type& spline() const { return __spline; }
};

}  // namespace intp
