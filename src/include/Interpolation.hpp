#include <cmath>

#include "BSpline.hpp"
#include "BandMatrix.hpp"

template <typename T>
class InterpolationFunction {
   private:
    using val_type = T;
    using spline_type = BSpline<T>;
    using size_type = typename spline_type::size_type;
    const size_type order;

    spline_type spline;
    bool _uniform;
    val_type _dx;

   public:
    /**
     * @brief Construct a new Interpolation Function object, mimicking
     * Mathematica's Interpolation function.
     *
     * @tparam InputIter
     * @param order order of interpolation, the interpolated function is
     * $C^{order-1}$
     * @param x_range a pair of x_min and x_max
     * @param f_begin iterator points to first value of to-be-interpolated data
     * @param f_end iterator points to one-past-end of to-be-interpolated data
     */
    template <typename InputIter>
    InterpolationFunction(size_type order,
                          std::pair<double, double> x_range,
                          InputIter f_begin,
                          InputIter f_end)
        : order(order), _uniform(true), spline(order) {
        const size_type n = std::distance(f_begin, f_end);
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

        spline.load_knots(std::move(xs));
        BandMatrix<typename spline_type::knot_type> coef_mat{n, order - 1,
                                                             order - 1};
        auto knots_iter = spline.knots_begin() + order;
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
                    knots_iter,
                    (i <= order ? (double)i : .5 * (order + 3 - 2 * (n - i))) *
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
                    spline.base_spline_value(knots_iter,
                                             (1 - order % 2) * 0.5 * _dx);
                }
                for (int j = i - order / 2; j < i + order / 2 + 1; ++j) {
                    coef_mat(i, j) = weights[j - i + order / 2];
                }
            }
        }

        spline.load_ctrlPts(std::move(coef_mat).linear_solve(
            std::vector<val_type>{f_begin, f_end}));
    }
    template <typename InputIterX, typename InputIterF>
    InterpolationFunction(InputIterX x_begin,
                          InputIterX x_end,
                          InputIterF f_begin,
                          InputIterF f_end)
        : _uniform(false) {}

    const std::pair<typename BSpline<T>::knot_type,
                    typename BSpline<T>::knot_type>&
    range() const {
        return spline.range();
    }

    val_type operator()(double x) {
        return _uniform ? spline(x, std::min(spline.knots_num() - order - 2,
                                             (size_type)std::ceil(std::max(
                                                 0., (x - range().first) / _dx -
                                                         .5 * (order + 1))) +
                                                 order))
                        : spline(x);
    }
};
