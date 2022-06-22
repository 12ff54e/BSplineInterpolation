#pragma once

#include <algorithm>  // upper_bound
#include <array>
#include <cmath>       // fmod
#include <functional>  // ref
#include <iterator>    // distance
#include <stdexcept>    // range_error
#include <type_traits>  // is_same, is_arithmatic
#include <vector>

#ifdef _DEBUG
#include <iostream>
#endif

#include "Mesh.hpp"
#include "util.hpp"

namespace intp {

/**
 * @brief B-Spline function
 *
 * @tparam T Type of control point
 * @tparam D Dimension
 */
template <typename T, size_t D>
class BSpline {
   public:
    using size_type = size_t;
    using val_type = T;
    using knot_type = double;

    using KnotContainer = std::vector<knot_type>;
    using ControlPointContainer = Mesh<val_type, D>;

    using BaseSpline = std::vector<knot_type>;

    const static size_type dim = D;
    const size_type order;

    // Container for dimension-wise storage
    template <typename _T>
    using DimArray = std::array<_T, dim>;

    /**
     * @brief Calculate values on base spline function. This is the core of
     * B-Spline. Note: when the given order is smaller than order of spline
     * (used in calculating derivative), spline value is aligned at right in
     * result vector.
     *
     * @param dim_ind which dimension
     * @param seg_idx_iter the iterator points to left knot point of a segment
     * @param x
     * @param spline_order order of base spline, defaulted to be spline function
     * order
     * @return a reference to local buffer
     */
    inline const BaseSpline& base_spline_value(
        size_type dim_ind,
        KnotContainer::const_iterator seg_idx_iter,
        knot_type x,
        size_type spline_order) const {
        std::fill(base_spline_buf.begin(), base_spline_buf.end(), 0);
        base_spline_buf[order] = 1;
        for (size_type i = 1; i <= spline_order; ++i) {
            // Each iteration will expand buffer zone by one, from back
            // to front.
            const size_type idx_begin = order - i;
            for (size_type j = 0; j <= i; ++j) {
                const auto left_iter = seg_idx_iter - (i - j);
                const auto right_iter = seg_idx_iter + (j + 1);
                base_spline_buf[idx_begin + j] =
                    (j == 0
                         ? 0
                         : base_spline_buf[idx_begin + j] * (x - *left_iter) /
                               (*(right_iter - 1) - *left_iter)) +
                    (idx_begin + j == order
                         ? 0
                         : base_spline_buf[idx_begin + j + 1] *
                               (*right_iter - x) /
                               (*right_iter - *(left_iter + 1)));
            }
        }
        // }
        return base_spline_buf;
    }

    inline const BaseSpline& base_spline_value(
        size_type dim_ind,
        KnotContainer::const_iterator seg_idx_iter,
        knot_type x) const {
        return base_spline_value(dim_ind, seg_idx_iter, x, order);
    }

    /**
     * @brief Return the result of last call for this method
     *
     * @return a reference to local buffer
     */
    inline const std::vector<knot_type>& base_spline_value() const {
        return base_spline_buf;
    }

    /**
     * @brief Get the lower knot iter points to the segment where given x
     * locates. If x is out of range of knot vector, the iterator is rather
     * begin or end of knot vector.
     *
     * @param dim_ind specify the dimension
     * @param x coordinate
     * @param first a hint for iter offset
     * @param last an upper bound for iter offset, this function will not search
     * knots beyond it.
     * @return KnotContainer::const_iterator
     */
    inline KnotContainer::const_iterator get_knot_iter(size_type dim_ind,
                                                       knot_type& x,
                                                       size_type hint,
                                                       size_type last) const {
        const auto iter = knots_begin(dim_ind) + hint;
        if (__periodicity[dim_ind]) {
            const knot_type period =
                range(dim_ind).second - range(dim_ind).first;
            x = range(dim_ind).first +
                (x < range(dim_ind).first
                     ? std::fmod(x - range(dim_ind).first, period) + period
                     : std::fmod(x - range(dim_ind).first, period));
        }
#ifdef _DEBUG
        if ((*iter > x || *(iter + 1) < x) && x >= range(dim_ind).first &&
            x <= range(dim_ind).second) {
            std::cout << "[DEBUG] knot hint miss at dim = " << dim_ind
                      << ", hint = " << hint << ", x = " << x << '\n';
        }
#endif
        return *iter <= x && *(iter + 1) > x
                   // If the hint is accurate, use that iter
                   ? iter
                   // else, use binary search in the range of distinct knots
                   // (excluding beginning and ending knots that have same
                   // value)
                   : --(std::upper_bound(knots_begin(dim_ind) + order + 1,
                                         knots_begin(dim_ind) + (last + 1), x));
    }

    inline KnotContainer::const_iterator get_knot_iter(size_type dim_ind,
                                                       knot_type& x,
                                                       size_type hint) const {
        return get_knot_iter(dim_ind, x, hint, knots_num(dim_ind) - order - 2);
    }

    template <typename... CoordWithHints, size_type... indices>
    inline DimArray<KnotContainer::const_iterator> get_knot_iters(
        util::index_sequence<indices...>,
        CoordWithHints&&... coords) const {
        return {get_knot_iter(indices, coords.first, coords.second)...};
    }

   private:
    DimArray<bool> __periodicity;
    DimArray<bool> __uniform;

    DimArray<KnotContainer> knots;
    ControlPointContainer control_points;

    mutable BaseSpline base_spline_buf;
    DimArray<std::pair<knot_type, knot_type>> __range;

    const size_type buf_size;

    // maximum stack buffer size
    // This buffer is for storing weights when calculating spline derivative
    // value.
    constexpr static size_type MAX_BUF_SIZE = 1000;

    // auxiliary methods

    /**
     * @brief Calculate base spline value of each dimension
     *
     * @tparam Coords knot_type ...
     * @tparam indices 0, 1, ...
     * @param knot_iters an array of knot iters
     * @param spline_order an array of spline order
     * @param coords a bunch of coordinates
     * @return std::array<decltype(base_spline_buf), dim>
     */
    template <typename... Coords, size_type... indices>
    inline DimArray<decltype(base_spline_buf)> calc_base_spline_vals(
        util::index_sequence<indices...>,
        const DimArray<KnotContainer::const_iterator>& knot_iters,
        const DimArray<size_type>& spline_order,
        Coords... coords) const {
        return {base_spline_value(indices, knot_iters[indices], coords,
                                  spline_order[indices])...};
    }

   public:
    /**
     * @brief Construct a new BSpline object, with periodicity of each dimension
     * specified.
     *
     * @param periodicity
     * @param order
     */
    explicit BSpline(DimArray<bool> periodicity, size_type order = 3)
        : order(order),
          __periodicity(periodicity),
          control_points(size_type{}),
          base_spline_buf(order + 1, 0),
          buf_size(util::pow(order + 1, dim)){};

    /**
     * @brief Basically the default constructor, initialize an empty, non-closed
     * B-Spline with order defaulted to be 3.
     *
     * @param order
     */
    BSpline(size_type order = 3) : BSpline(DimArray<bool>{}, order){};

    template <typename... InputIters>
    BSpline(size_type order,
            DimArray<bool> periodicity,
            ControlPointContainer ctrl_points,
            std::pair<InputIters, InputIters>... knot_iter_pairs)
        : order(order),
          __periodicity(periodicity),
          knots{
              KnotContainer(knot_iter_pairs.first, knot_iter_pairs.second)...},
          control_points(std::move(ctrl_points)),
          base_spline_buf(order + 1, 0),
          __range{std::make_pair(
              (knot_iter_pairs.first)[order],
              (knot_iter_pairs.second)[-static_cast<int>(order) - 1])...},
          buf_size(util::pow(order + 1, dim)) {
        for (size_type d = 0; d < dim; ++d) {
            if (knots[d].size() - control_points.dim_size(d) !=
                (__periodicity[d] ? 2 * order + 1 : order + 1)) {
                throw std::range_error(
                    "Inconsistency between knot number and control point "
                    "number.");
            }
        }
        __uniform.fill(true);
    }

    template <typename... InputIters>
    BSpline(size_type order,
            ControlPointContainer ctrl_points,
            std::pair<InputIters, InputIters>... knot_iter_pairs)
        : BSpline(order, DimArray<bool>{}, ctrl_points, knot_iter_pairs...) {}

    template <typename C>
    typename std::enable_if<
        std::is_same<typename std::remove_reference<C>::type,
                     KnotContainer>::value,
        void>::type
    load_knots(size_type dim_ind, C&& _knots, bool is_uniform = false) {
        knots[dim_ind] = std::forward<C>(_knots);
        __range[dim_ind].first = knots[dim_ind][order];
        __range[dim_ind].second =
            knots[dim_ind][knots[dim_ind].size() - order - (2 - order % 2)];
        __uniform[dim_ind] = is_uniform;
    }

    template <typename C>
    typename std::enable_if<
        std::is_same<typename std::remove_reference<C>::type,
                     ControlPointContainer>::value,
        void>::type
    load_ctrlPts(C&& _control_points) {
        control_points = std::forward<C>(_control_points);
    }

    /**
     * @brief Get spline value at given pairs of coordinate and position hint
     * (hopefully lower knot point index of the segment where coordinate
     * locates, dimension wise).
     *
     * @tparam CoordWithHints std::pair<knot_type, size_type>, ...
     * @param coord_with_hints a bunch of (coordinate, position hint) pairs
     * @return val_type
     */
    template <
        typename... CoordWithHints,
        typename Indices = util::make_index_sequence_for<CoordWithHints...>>
    typename std::enable_if<
        std::is_arithmetic<typename std::common_type<
            typename CoordWithHints::first_type...>::type>::value &&
            std::is_integral<typename std::common_type<
                typename CoordWithHints::second_type...>::type>::value,
        val_type>::type
    operator()(CoordWithHints... coord_with_hints) const {
        // get knot point iter, it will modifies coordinate value into
        // interpolation range of periodic dimension.
        const auto knot_iters = get_knot_iters(Indices{}, coord_with_hints...);

        DimArray<size_type> spline_order;
        spline_order.fill(order);
        // calculate basic spline (out of boundary check also conducted here)
        const auto base_spline_values_1d = calc_base_spline_vals(
            Indices{}, knot_iters, spline_order, coord_with_hints.first...);

        // combine control points and basic spline values to get spline value
        val_type v{};
        for (size_type i = 0; i < buf_size; ++i) {
            DimArray<size_type> ind_arr;
            for (size_type d = 0, combined_ind = i; d < dim; ++d) {
                ind_arr[d] = combined_ind % (order + 1);
                combined_ind /= (order + 1);
            }

            val_type coef = 1;
            for (size_type d = 0; d < dim; ++d) {
                coef *= base_spline_values_1d[d][ind_arr[d]];

                // Shift index array according to knot iter of each dimension.
                // When the coordinate is out of range in some dimensions, the
                // corresponding iterator was set to be begin or end iterator of
                // knot vector in `get_knot_iters` method and it will be treated
                // separately.
                ind_arr[d] +=
                    knot_iters[d] == knots_begin(d) ? 0
                    : knot_iters[d] == knots_end(d)
                        ? control_points.dim_size(d) - order - 1
                        : distance(knots_begin(d), knot_iters[d]) - order;

                // check periodicity, put out-of-right-boundary index to left
                if (__periodicity[d]) {
                    ind_arr[d] %= control_points.dim_size(d);
                }
            }

            v += coef * control_points(ind_arr);
        }

        return v;
    }

    /**
     * @brief Get spline value at given coordinates
     *
     * @tparam Coords some arithmetic type, ...
     * @param coords a bunch of cartesian coordinates
     * @return val_type
     */
    template <typename... Coords>
    typename std::enable_if<
        std::is_arithmetic<typename std::common_type<Coords...>::type>::value,
        val_type>::type
    operator()(Coords... coords) const {
        return operator()(std::make_pair((knot_type)coords, order)...);
    }

    /**
     * @brief Get derivative value at given pairs of coordinate and position
     * hint (possibly lower knot point index of the segment where coordinate
     * locates, dimension wise).
     *
     * @tparam CoordDeriOrderHintTuple std::tuple<knot_type, size_type,
     * size_type>, ...
     * @param coord_deriOrder_hint_tuple a bunch of (coordinate, derivative
     * order, position hint) tuple
     * @return val_type
     */
    template <typename... CoordDeriOrderHintTuple,
              typename Indices =
                  util::make_index_sequence_for<CoordDeriOrderHintTuple...>>
    typename std::enable_if<std::tuple_size<typename std::common_type<
                                CoordDeriOrderHintTuple...>::type>::value == 3,
                            val_type>::type
    derivative_at(CoordDeriOrderHintTuple... coord_deriOrder_hint_tuple) const {
        // get spline order
        DimArray<size_type> spline_order{
            (order >= std::get<1>(coord_deriOrder_hint_tuple)
                 ? order - std::get<1>(coord_deriOrder_hint_tuple)
                 : order + 1)...};
        // if derivative order is larger than spline order, derivative is 0.
        for (auto o : spline_order) {
            if (o > order) { return val_type{}; }
        }

        // get knot point iter
        const auto knot_iters = get_knot_iters(
            Indices{},
            std::make_pair(std::ref(std::get<0>(coord_deriOrder_hint_tuple)),
                           std::get<2>(coord_deriOrder_hint_tuple))...);

        // calculate basic spline (out of boundary check also conducted here)
        const auto base_spline_values_1d =
            calc_base_spline_vals(Indices{}, knot_iters, spline_order,
                                  std::get<0>(coord_deriOrder_hint_tuple)...);

#ifdef STACK_ALLOCATOR
        // create local buffer
        val_type buffer[MAX_BUF_SIZE];
        util::stack_allocator<val_type, MAX_BUF_SIZE> alloc(buffer);

        Mesh<val_type, dim, util::stack_allocator<val_type, MAX_BUF_SIZE>>
            local_control_points(order + 1, alloc);
        auto local_spline_val = local_control_points;
#else
        Mesh<val_type, dim> local_control_points(order + 1);
        auto local_spline_val = local_control_points;
#endif

        // get local control points and basic spline values
        for (size_type i = 0; i < buf_size; ++i) {
            DimArray<size_type> local_ind_arr{};
            for (size_type d = 0, combined_ind = i; d < dim; ++d) {
                local_ind_arr[d] = combined_ind % (order + 1);
                combined_ind /= (order + 1);
            }

            val_type coef = 1;
            DimArray<size_type> ind_arr{};
            for (size_type d = 0; d < dim; ++d) {
                coef *= base_spline_values_1d[d][local_ind_arr[d]];

                ind_arr[d] =
                    local_ind_arr[d] +
                    (knot_iters[d] == knots_begin(d) ? 0
                     : knot_iters[d] == knots_end(d)
                         ? control_points.dim_size(d) - order - 1
                         : distance(knots_begin(d), knot_iters[d]) - order);

                // check periodicity, put out-of-right-boundary index to left
                if (__periodicity[d]) {
                    ind_arr[d] %= (knots_num(d) - 2 * order - 1);
                }
            }

            local_spline_val(local_ind_arr) = coef;
            local_control_points(local_ind_arr) = control_points(ind_arr);
        }

        for (size_type d = 0; d < dim; ++d) {
            if (spline_order[d] == order) { continue; }
            // calculate control points for derivative along this dimension

            const size_type hyper_surface_size =
                local_control_points.size() / local_control_points.dim_size(d);
            // transverse the hyper surface of fixing dimension d
            for (size_type i = 0; i < hyper_surface_size; ++i) {
                DimArray<size_type> local_ind_arr;
                for (size_type dd = 0, combined_ind = i; dd < dim; ++dd) {
                    if (dd == d) { continue; }
                    local_ind_arr[dd] = combined_ind % (order + 1);
                    combined_ind /= (order + 1);
                }

                auto iter = local_control_points.begin(d, local_ind_arr);
                // Taking derivative is effectively computing new control
                // points. Number of iteration is order of derivative.
                for (int k = order; k > static_cast<int>(spline_order[d]);
                     --k) {
                    // Each reduction reduce control points number by one.
                    // Reduce backward to match pattern of local_spline_val.
                    for (int j = k; j > 0; --j) {
                        iter[order + j - k] =
                            k *
                            (iter[order + j - k] - iter[order + j - k - 1]) /
                            (knot_iters[d][j] - knot_iters[d][j - k]);
                    }
                }
            }
        }

        // combine spline value and control points to get spline derivative
        // value
        val_type v{};
        for (auto s_it = local_spline_val.begin(),
                  c_it = local_control_points.begin();
             s_it != local_spline_val.end(); ++s_it, ++c_it) {
            v += (*s_it) * (*c_it);
        }

        return v;
    }

    /**
     * @brief Get derivative value at given coordinates
     *
     * @tparam CoordDeriOrderPair std::pair<knot_type, size_type>, ...
     * @param coords a bunch of (coordinate, derivative order) tuple
     * @return val_type
     */
    template <typename... CoordDeriOrderPair>
    typename std::enable_if<std::tuple_size<typename std::common_type<
                                CoordDeriOrderPair...>::type>::value == 2,
                            val_type>::type
    derivative_at(CoordDeriOrderPair... coords) const {
        return derivative_at(std::make_tuple(
            (knot_type)coords.first, (size_type)coords.second, order)...);
    }

    // iterators

    /**
     * @brief Returns a read-only (constant) iterator that points to the
     *  first element in the knot vector of one dimension.
     *
     * @param dim_ind dimension index
     */
    inline KnotContainer::const_iterator knots_begin(size_type dim_ind) const {
        return knots[dim_ind].cbegin();
    }
    /**
     * @brief Returns a read-only (constant) iterator that points to the
     *  first element in the knot vector of one dimension.
     *
     * @param dim_ind dimension index
     */
    inline KnotContainer::const_iterator knots_end(size_type dim_ind) const {
        return knots[dim_ind].cend();
    }

    // properties

    /**
     * @brief Get range of one dimension
     *
     * @param dim_ind
     * @return const std::pair<knot_type, knot_type>&
     */
    inline const std::pair<knot_type, knot_type>& range(
        size_type dim_ind) const {
        return __range[dim_ind];
    }

    /**
     * @brief Get knot number of one dimension
     *
     * @param dim_ind
     * @return size_type
     */
    inline size_type knots_num(size_type dim_ind) const {
        return knots[dim_ind].size();
    }

    /**
     * @brief Get periodicity of one dimension
     *
     * @param dim_ind dimension index
     * @return a bool
     */
    bool periodicity(size_type dim_ind) const {
        return __periodicity[dim_ind];
    }

    /**
     * @brief Get uniformity of one dimension
     *
     * @param dim_ind dimension index
     * @return a bool
     */
    bool uniform(size_type dim_ind) const {
        return __uniform[dim_ind];
    }

#ifdef _DEBUG
    void __debug_output() const {
        std::cout << "\n[DEBUG] Control Points (raw data):\n";

        // 17 digits for double precision
        std::cout.precision(17);
        int idx = 1;
        for (auto v : control_points) {
            if (idx % control_points.dim_size(dim - 1) == 1) {
                std::cout << "[DEBUG] ";
            }
            std::cout << v << ' ';
            if (idx++ % control_points.dim_size(dim - 1) == 0) {
                std::cout << '\n';
            }
        }
        std::cout << '\n';
    }
#endif
};

}  // namespace intp
