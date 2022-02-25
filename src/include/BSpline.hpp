#pragma once

#include <algorithm>  // upper_bound
#include <array>
#include <iterator>     // advance, distance, etc.
#include <stdexcept>    // range_error
#include <type_traits>  // is_same, is_arithmatic
#include <vector>

#ifdef DEBUG
#include <iostream>
#endif

#include "util.hpp"

/**
 * @brief A multi dimension mesh storing data on each mesh point
 *
 * @tparam T Type of data stored
 * @tparam D Dimension
 */
template <typename T, unsigned D>
class Mesh {
   public:
    using size_type = unsigned int;
    using val_type = T;
    const static size_type dim = D;
    using Indices = std::array<size_type, dim>;

   private:
    using iterator = typename std::vector<val_type>::const_iterator;

    /**
     * @brief Stores the mesh content in row-major format.
     */
    std::vector<val_type> storage;
    /**
     * @brief The i-th element stores the sub-mesh size when the first (dim-i)
     * coordinates are specified.
     */
    std::array<size_type, dim + 1> __dim_acc_size;
    std::array<size_type, dim> __dim_size;

    // auxilary methods

    /**
     * @brief Set the __dim_acc_size object. Check description of __dim_acc_size
     * for details.
     */
    template <typename... DimSizes>
    void set_dim_acc_size(size_type dimSize, DimSizes&... dimSizes) {
        for (size_type i = sizeof...(dimSizes) + 1; i <= dim; ++i) {
            __dim_acc_size[i] *= dimSize;
        }
        set_dim_acc_size(dimSizes...);
        return;
    }

    constexpr void set_dim_acc_size() { return; }

    /**
     * @brief Convert multi-dimesion index to one dimension index in storage
     * vector.
     *
     * @param ind
     * @param indices
     * @return size_type
     */
    template <typename... _Indices>
    size_type indexing(size_type ind, _Indices... indices) const {
        return ind * __dim_acc_size[sizeof...(indices)] + indexing(indices...);
    }

    size_type indexing(Indices ind_arr) const {
        size_type ind{};
        for (size_type i = 0; i < dim; ++i) {
            ind += ind_arr[i] * __dim_acc_size[dim - i - 1];
        }
        return ind;
    }

    /**
     * @brief Convert one dimension index in storage vector to multi-dimension
     * indices
     *
     * @param total_ind
     * @return std::array<size_type, dim>
     */
    Indices dimwise_indices(size_type total_ind) {
        Indices indices;

        for (unsigned i = 0; i < dim; ++i) {
            indices[i] = total_ind / __dim_acc_size[dim - i - 1];
            total_ind %= __dim_acc_size[dim - i - 1];
        }

        return indices;
    }

    constexpr size_type indexing() const { return 0; }

    template <typename U, unsigned DD>
    friend class InterpolationFunction;  // friend class forward declaration,
                                         // InterpolationFunction need access
                                         // to indexing function during
                                         // interpolation procedure

   public:
    Mesh() = delete;

    template <typename... DimSizes>
    explicit Mesh(DimSizes... dimSizes) : __dim_size{(size_type)dimSizes...} {
        __dim_acc_size.fill(1u);
        set_dim_acc_size(dimSizes...);

        storage.resize(__dim_acc_size.back(), val_type{});
    }

    template <typename InputIter,
              typename = typename std::enable_if<
                  dim == 1u &&
                  std::is_convertible<typename std::iterator_traits<
                                          InputIter>::iterator_category,
                                      std::input_iterator_tag>::value>::type>
    explicit Mesh(InputIter iter_begin, InputIter iter_end)
        : storage(iter_begin, iter_end),
          __dim_size{(size_type)storage.size()},
          __dim_acc_size{1u, (size_type)storage.size()} {}

    template <typename Array,
              typename = typename std::enable_if<
                  dim == 1u && !std::is_scalar<Array>::value,
                  int>::type>
    explicit Mesh(const Array& array) : Mesh(array.begin(), array.end()) {}

   public:
    // properties

    size_type size() const { return storage.size(); }

    size_type dim_size(size_type dim_ind) const { return __dim_size[dim_ind]; }

    // modifiers

    void resize(Indices sizes) {
        __dim_size = sizes;
        for (size_type i = 0; i <= dim; ++i) {
            __dim_acc_size[i] = 1;
            for (size_type j = 0; j < i; ++j) {
                __dim_acc_size[i] *= sizes[dim - j - 1];
            }
        }
        storage.resize(__dim_acc_size.back());
    }

    // element access

    template <typename... _Indices>
    val_type& operator()(_Indices... indices) {
        return storage[indexing(indices...)];
    }

    template <typename... _Indices>
    val_type operator()(_Indices... indices) const {
        return storage[indexing(indices...)];
    }

    const val_type* data() const { return storage.data(); }

    // iterator

    iterator begin() const { return storage.cbegin(); }
    iterator end() const { return storage.cend(); }

    Indices iter_indices(iterator iter) {
        return dimwise_indices(std::distance(begin(), iter));
    }

    // debug

    const std::array<unsigned, dim + 1>& dim_acc_size() const {
        return __dim_acc_size;
    }
};

/**
 * @brief B-Spline function
 *
 * @tparam T Type of control point
 * @tparam D Dimension
 */
template <typename T, unsigned D>
class BSpline {
   public:
    using size_type = unsigned int;
    using val_type = T;
    using knot_type = double;

    using KnotContainer = std::vector<knot_type>;
    using ControlPointContainer = Mesh<val_type, D>;

    using BaseSpline = std::vector<knot_type>;

    const static size_type dim = D;
    const size_type order;

    /**
     * @brief Calculate values on base spline function. This is the core of
     * B-Spline.
     *
     * @param seg_idx_iter the iterator points to left knot point of a segment
     * @param x
     * @return a reference to local buffer
     */
    inline const BaseSpline& base_spline_value(
        size_type dim_ind,
        KnotContainer::const_iterator seg_idx_iter,
        knot_type x) const {
        std::fill(base_spline_buf.begin(), base_spline_buf.end(), 0);
        // out of boundary check
        if (seg_idx_iter == knots_begin(dim_ind)) {
            base_spline_buf[0] = 1.;
        } else if (seg_idx_iter == knots_end(dim_ind)) {
            base_spline_buf[order] = 1.;
        } else {
            base_spline_buf[order] = 1;
            for (size_type i = 1; i <= order; ++i) {
                // Each iteration will expand buffer zone by one, from back
                // to front.
                const size_type idx_begin = order - i;
                for (size_type j = 0; j <= i; ++j) {
                    const auto left_iter = seg_idx_iter - (i - j);
                    const auto right_iter = seg_idx_iter + (j + 1);
                    base_spline_buf[idx_begin + j] =
                        (j == 0 ? 0
                                : base_spline_buf[idx_begin + j] *
                                      (x - *left_iter) /
                                      (*(right_iter - 1) - *left_iter)) +
                        (idx_begin + j == order
                             ? 0
                             : base_spline_buf[idx_begin + j + 1] *
                                   (*right_iter - x) /
                                   (*right_iter - *(left_iter + 1)));
                }
            }
        }
        return base_spline_buf;
    }

    /**
     * @brief Return the result of last call for this method
     *
     * @return a reference to local buffer
     */
    inline const std::vector<knot_type>& base_spline_value() const {
        return base_spline_buf;
    }

    std::vector<knot_type> composite_spline_value() const {}

   private:
    template <typename _T>
    using DimArray = std::array<_T, dim>;

    DimArray<KnotContainer> knots;
    ControlPointContainer control_points;

    mutable BaseSpline base_spline_buf;
    DimArray<std::pair<knot_type, knot_type>> _range;

    DimArray<bool> _periodicity;
    DimArray<bool> _uniform;

    const size_type buf_size;
    // auxiliary methods

    /**
     * @brief Get the lower knot iter points to the segment where given x
     * locates. If x is out of range of knot vector, the iterator is rather
     * begin or end of knot vector.
     *
     * @param dim_ind specify the dimension
     * @param x coordinate
     * @param knot_ind a hint for iter offset
     * @return KnotContainer::const_iterator
     */
    inline KnotContainer::const_iterator
    get_knot_iter(size_type dim_ind, knot_type x, size_type knot_ind) const {
        if (x < knots[dim_ind].front()) {
            return knots_begin(dim_ind);
        } else if (x >= knots[dim_ind].back()) {
            return knots_end(dim_ind);
        }
        const auto iter = knots_begin(dim_ind) + knot_ind;
        return *iter <= x && *(iter + 1) >= x
                   // If the hint is accurate, use that iter
                   ? iter
                   // else, use binary search in the range of distint knots
                   // (excluding beginning and ending knots that have same
                   // value)
                   : prev(std::upper_bound(knots_begin(dim_ind) + order,
                                           knots_end(dim_ind), x));
    }

    template <typename... CoordWithHints, unsigned... indices>
    inline DimArray<KnotContainer::const_iterator> get_knot_iters(
        util::index_sequence<indices...>,
        const CoordWithHints&... coords) const {
        return {get_knot_iter(indices, coords.first, coords.second)...};
    }

    /**
     * @brief Calculate base spline value of each dimension
     *
     * @tparam Coords std::pair<knot_type, size_type>, ...
     * @tparam indices 0, 1, ...
     * @param knot_iters an array of knot iters
     * @param coords a bunch of coordinates
     * @return std::array<decltype(base_spline_buf), dim>
     */
    template <typename... Coords, unsigned... indices>
    inline DimArray<decltype(base_spline_buf)> calc_base_spline_vals(
        util::index_sequence<indices...>,
        const DimArray<KnotContainer::const_iterator>& knot_iters,
        Coords... coords) const {
        return {base_spline_value(indices, knot_iters[indices], coords)...};
    }

   public:
    BSpline(DimArray<bool> periodicity, size_type order = 3)
        : order(order),
          base_spline_buf(order + 1, 0),
          control_points(size_type{}),
          buf_size(util::pow(order + 1, dim)),
          _periodicity(periodicity){};

    BSpline(size_type order = 3) : BSpline(DimArray<bool>{}, order){};

    template <typename... InputIters>
    BSpline(size_type order,
            DimArray<bool> periodicity,
            ControlPointContainer ctrl_points,
            std::pair<InputIters, InputIters>... knot_iter_pairs)
        : order(order),
          base_spline_buf(order + 1, 0),
          control_points(std::move(ctrl_points)),
          knots{
              KnotContainer(knot_iter_pairs.first, knot_iter_pairs.second)...},
          buf_size(util::pow(order + 1, dim)),
          _periodicity(periodicity),
          _range{std::make_pair(*(knot_iter_pairs.first),
                                *(knot_iter_pairs.second - 1))...} {
        for (int i = 0; i < dim; ++i) {
            if (knots[i].size() - control_points.dim_size(i) !=
                (_periodicity[i] ? 2 * order + 1 : order + 1)) {
                throw std::range_error(
                    "Inconsistency between knot number and control point "
                    "number.");
            }
        }
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
    load_knots(size_type dim_ind, C&& _knots) {
        knots[dim_ind] = std::forward<C>(_knots);
        _range[dim_ind].first = knots[dim_ind].front();
        _range[dim_ind].second = knots[dim_ind].back();
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
     * (possiblily lower knot point index of the segment where coordinate
     * locates, dimension wise).
     *
     * @tparam CoordWithHints std::pair<knot_type, size_type>, ...
     * @param coords a bunch of (coordinate, position hint) pairs
     * @return val_type
     */
    template <
        typename... CoordWithHints,
        typename Indices = util::make_index_sequence_for<CoordWithHints...>>
    typename std::enable_if<!std::is_arithmetic<typename std::common_type<
                                CoordWithHints...>::type>::value,
                            val_type>::type
    operator()(CoordWithHints... coord_with_hints) const {
        // get knot point iter
        const auto knot_iters = get_knot_iters(Indices{}, coord_with_hints...);

        // calculate basic spline (out of boundary check also conducted here)
        const auto base_spline_values_1d = calc_base_spline_vals(
            Indices{}, knot_iters, coord_with_hints.first...);

        // combine control points and basic spline values to get spline value
        val_type v{};
        for (size_type i = 0; i < buf_size; ++i) {
            DimArray<unsigned> ind_arr;
            for (size_type j = 0, combined_ind = i; j < dim; ++j) {
                ind_arr[j] = combined_ind % (order + 1);
                combined_ind /= (order + 1);
            }

            val_type coef = 1;
            for (size_type j = 0; j < dim; ++j) {
                coef *= base_spline_values_1d[j][ind_arr[j]];

                // Shift index array according to knot iter of each dimension.
                // When the coordinate is out of range in some dimensions, the
                // corresponding iterator was set to be begin or end iterator of
                // knot vector in `get_knot_iters` method and it will treated
                // seperately.
                ind_arr[j] +=
                    knot_iters[j] == knots_begin(j) ? 0
                    : knot_iters[j] == knots_end(j)
                        ? control_points.dim_size(j) - order - 1
                        : distance(knots_begin(j), knot_iters[j]) - order;

                // check periodicity, wrap the boundary
                if (_periodicity[j]) {
                    ind_arr[j] %= (knots_num(j) - 2 * order - 1);
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

    // iterators

    inline KnotContainer::const_iterator knots_begin(size_type dim_ind) const {
        return knots[dim_ind].cbegin();
    }
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
    const std::pair<knot_type, knot_type>& range(size_type dim_ind) const {
        return _range[dim_ind];
    }

    /**
     * @brief Get knot number of one dimension
     *
     * @param dim_ind
     * @return size_type
     */
    size_type knots_num(size_type dim_ind) const {
        return knots[dim_ind].size();
    }

    bool periodicity(size_type dim_ind) const { return _periodicity[dim_ind]; }

    bool uniform(size_type dim_ind) const { return _uniform[dim_ind]; }

#ifdef DEBUG
    void __debug_output() const {
        std::cout << "\n[DEBUG] Control points (raw data):\n";

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
