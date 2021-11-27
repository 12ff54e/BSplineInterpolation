#pragma once

#include <algorithm>
#include <array>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

// Debug
#include <iostream>

/**
 * @brief Calculate cardinal B-spline value on knots recursively (for constexpr
 * sake). Can be modified to iterative method in C++17, where std::array is
 * fully constexpr
 *
 * @param n order of spline
 * @param i index of knot
 * @return constexpr double
 */
constexpr double __spline_knot_val(unsigned n, unsigned i) {
    return n == 1 ? 1.
                  : (i == 0 || i == n - 1
                         ? __spline_knot_val(n - 1, 0) / n
                         : ((i + 1) * __spline_knot_val(n - 1, i) +
                            (n - i) * __spline_knot_val(n - 1, i - 1)) /
                               n);
}

/**
 * @brief Calculate cardinal B-spline value on knots during runtime,
 * iteratively.
 *
 * @param n order of spline
 * @param b number of points outside left boundary
 * @return std::vector<double>
 */
inline std::vector<double> __spline_knot_val_arr_rt(unsigned n, unsigned b) {
    std::vector<double> _spline_knot_vals(n, 0);
    if (b <= n) {
        _spline_knot_vals[n - 1] = 1.;
        for (int j = 2; j <= n; ++j) {
            const int _left_idx = n - j;
            const int _l_offset = _left_idx < b ? b - _left_idx : 0;
            const int _r_offset = _left_idx < b - 1 ? b - _left_idx - 1 : 0;
            for (int k = _left_idx + _r_offset; k < n; ++k) {
                _spline_knot_vals[k] =
                    (_spline_knot_vals[k] * (n - k) / (j - _r_offset) +
                     (k == n - 1 ? 0.
                                 : _spline_knot_vals[k + 1] *
                                       (k - _left_idx + 1 - _l_offset) /
                                       (j - _l_offset)));
            }
        }
    }
    return _spline_knot_vals;
}

// recursive ends, returns all parameters as an array
template <unsigned Length, typename... Vals>
constexpr typename std::enable_if<sizeof...(Vals) == (Length - 1),
                                  std::array<double, Length>>::type
__spline_knot_val_arr_cmpl(Vals... values) {
    return {values..., __spline_knot_val(Length, sizeof...(Vals))};
}
// construct array by recursively append values at the end of parameter list
template <unsigned Length, typename... Vals>
constexpr typename std::enable_if<sizeof...(Vals) < (Length - 1),
                                  std::array<double, Length>>::type
__spline_knot_val_arr_cmpl(Vals... values) {
    return __spline_knot_val_arr_cmpl<Length, Vals..., double>(
        values..., __spline_knot_val(Length, sizeof...(Vals)));
}

template <typename T>
class BSpline {
   public:
    using size_type = unsigned int;
    using val_type = T;
    using knot_type = double;

    using KnotContainer = std::vector<knot_type>;
    using ControlPointContainer = std::vector<val_type>;

    const size_type order;

   private:
    KnotContainer knots;
    ControlPointContainer control_points;

    mutable std::vector<double> buf;

   public:
    BSpline(size_type order = 3) : order(order), buf(order + 1, 0){};
    template <typename KnotInputIter, typename ControlPointInputIter>
    BSpline(size_type order,
            KnotInputIter knot_begin,
            KnotInputIter knot_end,
            ControlPointInputIter ctrl_pt_begin,
            ControlPointInputIter ctrl_pt_end)
        : order(order),
          buf(order + 1, 0),
          knots(knot_begin, knot_end),
          control_points(ctrl_pt_begin, ctrl_pt_end) {
        if (knots.size() - control_points.size() != order + 1) {
            throw std::range_error(
                "Inconsistency between knots number and control point number.");
        }
    }

    template <typename C>
    typename std::enable_if<
        std::is_same<typename std::remove_reference<C>::type,
                     KnotContainer>::value,
        void>::type
    load_knots(C&& _knots) {
        knots = std::forward<C>(_knots);
    }

    template <typename C>
    typename std::enable_if<
        std::is_same<typename std::remove_reference<C>::type,
                     ControlPointContainer>::value,
        void>::type
    load_ctrlPts(C&& _control_points) {
        control_points = std::forward<C>(_control_points);
    }

    inline const std::vector<knot_type>& base_spline_value(
        KnotContainer::const_iterator seg_idx_iter,
        knot_type x_offset) const {
        if (*seg_idx_iter == knots.back()) {
            // Special case, only used when calculating control point in
            // constructing interpolation function.
            std::fill(buf.begin(), buf.end(), 0);
            buf[order - 1] = 1.;
        } else {
            buf[order] = 1;
            for (size_type i = 1; i <= order; ++i) {
                // Each iteration will expand buffer zone by one, from back to
                // front.
                const size_type idx_begin = order - i;
                for (size_type j = 0; j <= i; ++j) {
                    const auto left_iter = seg_idx_iter - (i - j);
                    const auto right_iter = seg_idx_iter + (j + 1);
                    buf[idx_begin + j] =
                        (j == 0 ? 0
                                : buf[idx_begin + j] *
                                      (*seg_idx_iter - *left_iter + x_offset) /
                                      (*(right_iter - 1) - *left_iter)) +
                        (idx_begin + j == order
                             ? 0
                             : buf[idx_begin + j + 1] *
                                   (*right_iter - *seg_idx_iter - x_offset) /
                                   (*right_iter - *(left_iter + 1)));
                }
            }
        }
        return buf;
    }

    inline const std::vector<knot_type>& base_spline_value() const {
        return buf;
    }

    /**
     * @brief Calculate spline value at x, use a position hint indicating a
     * knot being left of the segment where x might located at.
     *
     * @param x
     * @param pos_hint
     * @return val_type
     */
    val_type operator()(val_type x, size_type pos_hint) const {
        // TODO: Change to higher order extrapolation
        if (x < knots.front()) {
            return control_points.front();
        } else if (x >= knots.back()) {
            return control_points.back();
        }

        const auto temp_iter =
            std::upper_bound(knots.begin() + pos_hint, knots.end(), x);
        const auto seg_idx_iter =
            prev(temp_iter == knots.end()
                     ? std::upper_bound(knots.begin(), knots.end(), x)
                     : temp_iter);
        const size_type seg_idx = distance(knots.begin(), seg_idx_iter);

        base_spline_value(seg_idx_iter,
                          x - *seg_idx_iter);  // This method modifies buf

        val_type v{};
        for (size_type i = 0; i <= order; ++i) {
            v += control_points[seg_idx - order + i] * buf[i];
        }
        return v;
    }

    val_type operator()(val_type x) const { return operator()(x, 0); }

    // iterators

    KnotContainer::const_iterator knots_begin() const { return knots.cbegin(); }
    KnotContainer::const_iterator knots_end() const { return knots.cend(); }

    void __debug_output() {
        auto end = knots.end() - order;

        std::cout << order << " base spline value on knot point:\n";
        for (auto iter = knots.begin() + order; iter != end; ++iter) {
            std::cout << "Index@" << *iter << " knot point: ";
            base_spline_value(iter, 0);
            std::copy(buf.begin(), buf.end(),
                      std::ostream_iterator<knot_type>(std::cout, " "));
            std::cout << '\n';
        }
    }
};
