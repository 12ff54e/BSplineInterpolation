#include <cmath>  // ceil

#include "BSpline.hpp"
#include "BandMatrix.hpp"

template <typename T, unsigned D>
class InterpolationFunction {
   private:
    using val_type = T;
    using spline_type = BSpline<T, D>;
    using size_type = typename spline_type::size_type;

    const size_type order;
    const static size_type dim = D;

    spline_type spline;
    bool _uniform;
    val_type _dx;

    // auxiliary methods

    template <typename... Coords, unsigned... indices>
    val_type call_op_helper(util::index_sequence<indices...>,
                            Coords... coords) {
        return _uniform
                   ? spline(std::make_pair(
                         coords,
                         std::min(
                             spline.knots_num(indices) - order - 2,
                             (size_type)std::ceil(std::max(
                                 0., (coords - range(indices).first) / _dx -
                                         .5 * (order + 1))) +
                                 order))...)
                   : spline(coords...);
    }

   public:
    /**
     * @brief Construct a new Interpolation Function object, mimicking
     * Mathematica's `Interpolation` function, with option `Method->"Spline"`.
     *
     * @tparam InputIter
     * @param order order of interpolation, the interpolated function is of
     * $C^{order-1}$
     * @param f_range a pair of iterators defining to-be-interpolated data
     * @param x_range a pair of x_min and x_max
     */
    template <typename InputIter>
    InterpolationFunction(size_type order,
                          std::pair<InputIter, InputIter> f_range,
                          std::pair<double, double> x_range)
        : order(order), _uniform(true), spline(order) {
        static_assert(dim == 1u,
                      "This constructor can only be used in 1D interpolation");

        const size_type n = std::distance(f_range.first, f_range.second);
        if (n < order + 1) { throw "Requested order is too high!"; }

        _dx = (x_range.second - x_range.first) / (n - 1);
        std::vector<typename spline_type::knot_type> xs(n + order + 1,
                                                        x_range.first);

        for (int i = order + 1; i < xs.size() - order - 1; ++i) {
            xs[i] = x_range.first + .5 * _dx * (2 * i - order - 1);
        }
        for (int i = xs.size() - order - 1; i < xs.size(); ++i) {
            xs[i] = x_range.second;
        }

        spline.load_knots(0, std::move(xs));
        BandMatrix<typename spline_type::knot_type> coef_mat{n, order - 1,
                                                             order - 1};
        auto knots_iter = spline.knots_begin(0) + order;
        auto& weights = spline.base_spline_value();
        for (int i = 0; i < n; ++i) {
            if (i == 0 || i == n - 1) {
                // base spline at boundary spanning one segment and their value
                // on endpoints are always 1
                coef_mat(i, i) = 1;
            } else if (i <= order / 2 || i >= n - order / 2 - 1) {
                // there are floor(order/2) point(s) in first (last) segment
                if (order % 2 == 0 && i == n - order / 2 - 1) { knots_iter++; }
                spline.base_spline_value(
                    0, knots_iter,
                    *knots_iter + (i <= order
                                       ? (double)i
                                       : .5 * (order + 3 - 2 * (n - i))) *
                                      _dx);
                for (int j = 0; j < order + 1; ++j) {
                    coef_mat(i, i <= order / 2 ? j : j + n - order - 1) =
                        weights[j];
                }
            } else {
                //
                ++knots_iter;
                if (i <= 1 + 3 * order / 2 || i >= n - (2 + 3 * order / 2)) {
                    // base spline function near boundary has different shape
                    spline.base_spline_value(
                        0, knots_iter,
                        *knots_iter + (1 - order % 2) * 0.5 * _dx);
                }
                for (int j = i - order / 2; j < i + order / 2 + 1; ++j) {
                    coef_mat(i, j) = weights[j - i + order / 2];
                }
            }
        }

        spline.load_ctrlPts(Mesh<val_type, 1>{std::move(coef_mat).linear_solve(
            std::vector<val_type>{f_range.first, f_range.second})});
    }

    const std::pair<typename BSpline<T, dim>::knot_type,
                    typename BSpline<T, dim>::knot_type>&
    range(size_type dim_ind) const {
        return spline.range(dim_ind);
    }

    template <typename... Coords,
              typename Indices = util::make_index_sequence_for<Coords...>>
    val_type operator()(Coords... x) {
        return call_op_helper(Indices{}, x...);
    }
};
