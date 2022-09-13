#pragma once

#include <cmath>  // ceil
#include <initializer_list>

#include "InterpolationTemplate.hpp"

namespace intp {

template <typename T, size_t D>
class InterpolationFunction {  // TODO: Add integration
   public:
    using val_type = T;
    using spline_type = BSpline<T, D>;
    using size_type = typename spline_type::size_type;
    using coord_type = typename spline_type::knot_type;

    const size_type order;
    const static size_type dim = D;

   private:
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
        // The x_range may be given by input iterators, which can not be
        // multi-passed.
        if (__periodicity[dim_ind]) {
            // In periodic case, the knots are data points, shifted by half of
            // local grid size if spline order is odd.

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
            // In aperiodic case, the internal knots are moving average of data
            // points with windows size equal to spline order.

            auto it = x_range.first;
            auto l_knot = *it;
            // fill leftmost *order+1* identical knots
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
#ifdef _TRACE
        std::cout << "[TRACE] Nonuniform knots along dimension" << dim_ind
                  << ":\n";
        for (auto& c : xs) { std::cout << "[TRACE] " << c << '\n'; }
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
        (__create_knot_vector(di, mesh_dimension, input_coords, x_ranges), ...);
#else
        // polyfill of C++17 fold expression over comma
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
        std::array<std::nullptr_t, sizeof...(Ts)>{
            (__create_knot_vector(di, mesh_dimension, input_coords, x_ranges),
             nullptr)...};
#pragma GCC diagnostic pop
#endif
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
        return call_op_helper(
            Indices{}, DimArray<coord_type>{static_cast<coord_type>(x)...});
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

template <typename T = double>
class InterpolationFunction1D : public InterpolationFunction<T, size_t{1}> {
   private:
    using base = InterpolationFunction<T, size_t{1}>;

   public:
    template <typename InputIter>
    InterpolationFunction1D(std::pair<InputIter, InputIter> f_range,
                            typename base::size_type order = 3,
                            bool periodicity = false)
        : InterpolationFunction1D(
              std::make_pair(typename base::coord_type{},
                             static_cast<typename base::coord_type>(
                                 (f_range.second - f_range.first) - 1)),
              f_range,
              order,
              periodicity){};

    template <typename C1, typename C2, typename InputIter>
    InterpolationFunction1D(std::pair<C1, C2> x_range,
                            std::pair<InputIter, InputIter> f_range,
                            typename base::size_type order = 3,
                            bool periodicity = false)
        : base(order, periodicity, f_range, x_range){};
};

}  // namespace intp
