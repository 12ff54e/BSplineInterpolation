#include <cmath>  // ceil

#include <Eigen/SparseLU>

#include "BSpline.hpp"

template <typename T, unsigned D>
class InterpolationFunction {
   private:
    using val_type = T;
    using spline_type = BSpline<T, D>;
    using size_type = typename spline_type::size_type;

    const size_type order;
    const static size_type dim = D;

    spline_type _spline;
    bool _uniform;
    std::array<val_type, dim> _dx;
    std::array<bool, dim> _periodicity;

    // auxiliary methods

    template <typename... Coords, unsigned... indices>
    val_type call_op_helper(util::index_sequence<indices...>,
                            Coords... coords) const {
        return _uniform ? _spline(std::make_pair(
                              coords,
                              std::min(_spline.knots_num(indices) - order - 2,
                                       (size_type)std::ceil(std::max(
                                           0., (coords - range(indices).first) /
                                                       _dx[indices] -
                                                   .5 * (order + 1))) +
                                           order))...)
                        : _spline(coords...);
    }

    template <typename Arr, size_type... indices>
    val_type call_op_arr_helper(util::index_sequence<indices...> i,
                                Arr& x) const {
        return call_op_helper(i, x[indices]...);
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
                          std::array<bool, dim> periodicity,
                          Mesh<val_type, dim> f_mesh,
                          std::pair<Ts, Ts>... x_ranges)
        : order(order),
          _uniform(true),
          _spline(periodicity, order),
          _periodicity(periodicity) {
        // load knots into spline
        util::dispatch_indexed(
            [&](size_type dim_ind,
                std::pair<typename spline_type::knot_type,
                          typename spline_type::knot_type> x_range) {
                const size_type n = f_mesh.dim_size(dim_ind);
                _dx[dim_ind] = (x_range.second - x_range.first) / (n - 1);

                std::vector<typename spline_type::knot_type> xs(
                    _periodicity[dim_ind] ? n + 2 * order : n + order + 1,
                    x_range.first);

                if (_periodicity[dim_ind]) {
                    for (int i = 0; i < xs.size(); ++i) {
                        xs[i] = x_range.first + (i - (int)order) * _dx[dim_ind];
                    }
                } else {
                    for (int i = order + 1; i < xs.size() - order - 1; ++i) {
                        xs[i] = x_range.first +
                                .5 * (2 * i - order - 1) * _dx[dim_ind];
                    }
                    for (int i = xs.size() - order - 1; i < xs.size(); ++i) {
                        xs[i] = x_range.second;
                    }
                }

                _spline.load_knots(dim_ind, std::move(xs));
            },
            x_ranges...);

        // initialize mesh storing weights of spline, adjust dimension according
        // to periodicity
        Mesh<val_type, dim> weights{f_mesh};
        std::array<size_type, dim> dim_size_tmp;
        bool p_flag = false;
        for (size_type i = 0; i < dim; ++i) {
            dim_size_tmp[i] = weights.dim_size(i) -
                              (_periodicity[i] ? ((p_flag = true), 1) : 0);
        }
        if (p_flag) { weights.resize(dim_size_tmp); }

        Eigen::SparseMatrix<double, Eigen::RowMajor> coef(weights.size(),
                                                          weights.size());
        // numbers of base spline covering
        std::vector<size_type> entry_per_row(f_mesh.size());
        for (size_type total_ind = 0; total_ind < f_mesh.size(); ++total_ind) {
            auto indices = f_mesh.dimwise_indices(total_ind);
            size_type c = 1;
            for (size_type d = 0; d < dim; ++d) {
                auto ind = indices[d];
                c *= _periodicity[d]                             ? order
                     : ind == 0 || ind == f_mesh.dim_size(d) - 1 ? 1
                     : ind <= order / 2 ||
                             ind >= f_mesh.dim_size(d) - 1 - order / 2
                         ? order + 1
                         : order | 1;
            }
            entry_per_row[total_ind] = c;
        }
        coef.reserve(entry_per_row);
        Eigen::VectorXd mesh_val(weights.size());

        std::array<typename spline_type::BaseSpline, dim>
            base_spline_vals_per_dim;

        // pre-calculate base spline of periodic dimension, since it never
        // changes due to its even-spaced knots
        for (size_type d = 0; d < dim; ++d) {
            if (_periodicity[d]) {
                base_spline_vals_per_dim[d] = _spline.base_spline_value(
                    d, _spline.knots_begin(d) + order, 0.);
            }
        }

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
            for (int i = 0; i < dim; ++i) {
                const auto knot_num = _spline.knots_num(i);
                // This is the index of i-th dimension knot vector to the left
                // of current f_indices[i] position, notice that knot points has
                // a larger gap in both ends in non-periodic case.
                const auto knot_ind =
                    _periodicity[i]
                        ? f_indices[i] + order
                        : std::min(knot_num - order - 2,
                                   f_indices[i] > order / 2
                                       ? f_indices[i] + (order + 1) / 2
                                       : order);
                if (!_periodicity[i]) {
                    if (knot_ind <= 2 * order + 1 ||
                        knot_ind >= knot_num - 2 * order - 2) {
                        // update base spline
                        const auto iter = _spline.knots_begin(i) + knot_ind;
                        const double x =
                            _spline.range(i).first + f_indices[i] * _dx[i];
                        base_spline_vals_per_dim[i] =
                            _spline.base_spline_value(i, iter, x);
                    }
                } else {
                    if (f_indices[i] == f_mesh.dim_size(i) - 1) {
                        skip_flag = true;
                        break;
                    }
                }
                base_spline_anchor[i] = knot_ind - order;
            }
            if (skip_flag) { continue; }

            // loop over nD base splines that contributes to current f_mesh
            // point, fill matrix
            for (int i = 0; i < util::pow(order + 1, dim); ++i) {
                std::array<size_type, dim> ind_arr;
                val_type spline_val = 1;
                for (int j = dim - 1, local_ind = i; j >= 0; --j) {
                    ind_arr[j] = local_ind % (order + 1);
                    local_ind /= (order + 1);

                    spline_val *= base_spline_vals_per_dim[j][ind_arr[j]];
                    ind_arr[j] += base_spline_anchor[j];

                    if (_periodicity[j]) {
                        ind_arr[j] %= (_spline.knots_num(j) - 2 * order - 1);
                    }
                }
                if (spline_val != val_type{0}) {
#ifdef _DEBUG
                    std::cout << "[DEBUG] {" << actual_ind << ','
                              << weights.indexing(ind_arr) << "} -> "
                              << spline_val << '\n';
#endif

                    coef.insert(actual_ind, weights.indexing(ind_arr)) =
                        spline_val;
                }
            }

            // fill rhs vector
            mesh_val(actual_ind++) = *it;
        }

        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
            solver;
        solver.compute(coef);

        Eigen::VectorXd weights_vec;
        weights_vec = solver.solve(mesh_val);

        weights.storage =
            std::vector<val_type>(weights_vec.begin(), weights_vec.end());

        _spline.load_ctrlPts(std::move(weights));
    }

    template <typename... Ts>
    InterpolationFunction(size_type order,
                          Mesh<val_type, dim> f_mesh,
                          std::pair<Ts, Ts>... x_ranges)
        : InterpolationFunction(order, {}, f_mesh, x_ranges...) {}

    const std::pair<typename spline_type::knot_type,
                    typename spline_type::knot_type>&
    range(size_type dim_ind) const {
        return _spline.range(dim_ind);
    }

    template <typename... Coords,
              typename Indices = util::make_index_sequence_for<Coords...>,
              typename = typename std::enable_if<std::is_scalar<
                  typename std::common_type<Coords...>::type>::value>::type>
    val_type operator()(Coords... x) const {
        return call_op_helper(Indices{}, x...);
    }

    template <
        typename Arr,
        typename = typename std::enable_if<!std::is_scalar<Arr>::value>::type>
    val_type operator()(Arr& x) const {
        using Indices = util::make_index_sequence<dim>;
        return call_op_arr_helper(Indices{}, x);
    }

    bool periodicity(size_type dim_ind) const { return _periodicity[dim_ind]; }

    /**
     * @brief Get a ref of underlying spline object
     *
     * @return spline_type&
     */
    const spline_type& spline() const { return _spline; }
};
