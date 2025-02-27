#ifndef INTP_INTERPOLATION
#define INTP_INTERPOLATION

#include <cmath>  // ceil
#include <initializer_list>

#ifdef INTP_TRACE
#ifndef INTP_DEBUG
#define INTP_DEBUG
#endif
#endif

#include "InterpolationTemplate.hpp"

namespace intp {

template <typename T, std::size_t D, std::size_t O, typename U>
class InterpolationFunction {
   public:
    using spline_type = BSpline<T, D, O, U>;

    using val_type = typename spline_type::val_type;
    using size_type = typename spline_type::size_type;
    using coord_type = typename spline_type::knot_type;
    using diff_type = typename spline_type::diff_type;

    constexpr static size_type dim = D;
    constexpr static size_type order = O;

    using function_type =
        InterpolationFunction<val_type, dim, order, coord_type>;

    template <typename T_>
    using DimArray = std::array<T_, dim>;

    friend class InterpolationFunctionTemplate<val_type,
                                               dim,
                                               order,
                                               coord_type>;

    /**
     * @brief Construct a new 1D Interpolation Function object, mimicking
     * Mathematica's `Interpolation` function, with option `Method->"Spline"`.
     *
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
    InterpolationFunction(bool periodic,
                          std::pair<InputIter, InputIter> f_range,
                          std::pair<C1, C2> x_range)
        : InterpolationFunction(
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
    InterpolationFunction(std::pair<InputIter, InputIter> f_range,
                          std::pair<C1, C2> x_range)
        : InterpolationFunction(false, f_range, x_range) {}

    /**
     * @brief Construct a new nD Interpolation Function object, mimicking
     * Mathematica's `Interpolation` function, with option `Method->"Spline"`.
     * Notice: last value of periodic dimension will be discarded since it is
     * considered same of the first value. Thus inconsistency input data for
     * periodic interpolation will be accepted.
     *
     * @param periodicity an array describing periodicity of each dimension
     * @param f_mesh a mesh containing data to be interpolated
     * @param x_ranges pairs of x_min and x_max or begin and end iterator
     */
    template <typename... Ts>
    InterpolationFunction(DimArray<bool> periodicity,
                          const Mesh<val_type, dim>& f_mesh,
                          std::pair<Ts, Ts>... x_ranges)
        : InterpolationFunction(
              InterpolationFunctionTemplate<val_type, dim, order, coord_type>(
                  periodicity,
                  f_mesh.dimension(),
                  x_ranges...)
                  .interpolate(f_mesh)) {}

    // Non-periodic for all dimension
    template <typename... Ts>
    InterpolationFunction(const Mesh<val_type, dim>& f_mesh,
                          std::pair<Ts, Ts>... x_ranges)
        : InterpolationFunction({}, f_mesh, x_ranges...) {}

    // constructor for partial construction, that is, without interpolated
    // values
    template <typename... Ts>
    InterpolationFunction(
        DimArray<bool> periodicity,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        MeshDimension<dim> mesh_dimension,
        std::pair<Ts, Ts>... x_ranges)
        : spline_(periodicity) {
        // load knots into spline
        create_knots_(util::make_index_sequence_for<Ts...>{},
                      std::move(mesh_dimension), input_coords, x_ranges...);
    }

    // An empty interpolation function, can only be used after populated by
    // InterpolationTemplate::interpolate().
    InterpolationFunction() = default;

    /**
     * @brief Get spline value.
     *
     * @param x coordinates
     */
    template <typename... Coords,
              typename = typename std::enable_if<std::is_arithmetic<
                  typename std::common_type<Coords...>::type>::value>::type>
    val_type operator()(Coords... x) const {
        return call_op_helper({static_cast<coord_type>(x)...});
    }

    /**
     * @brief Get spline value.
     *
     * @param coord coordinate array
     */
    val_type operator()(DimArray<coord_type> coord) const {
        return call_op_helper(coord);
    }

    /**
     * @brief Get spline value, but with out of boundary check.
     *
     * @param coord coordinate array
     */
    val_type at(DimArray<coord_type> coord) const {
        boundary_check_(coord);
        return call_op_helper(coord);
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
        return derivative(
            coord, DimArray<size_type>{static_cast<size_type>(deriOrder)...});
    }

    /**
     * @brief Get spline derivative value.
     *
     * @param coord_deriOrder_pair pairs of coordinate and derivative order
     */
    template <typename... CoordDeriOrderPair>
    val_type derivative(CoordDeriOrderPair... coord_deriOrder_pair) const {
        return derivative(DimArray<coord_type>{coord_deriOrder_pair.first...},
                          DimArray<size_type>{static_cast<size_type>(
                              coord_deriOrder_pair.second)...});
    }

    /**
     * @brief Get spline derivative value, but with out of boundary check.
     *
     * @param coord coordinate array
     * @param derivatives derivative order array
     */
    val_type derivative_at(DimArray<coord_type> coord,
                           DimArray<size_type> derivatives) const {
        boundary_check_(coord);
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
        return derivative_at(
            coord, DimArray<size_type>{static_cast<size_type>(deriOrder)...});
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
            DimArray<size_type>{
                static_cast<size_type>(coord_deriOrder_pair.second)...});
    }

    // properties

    inline bool periodicity(size_type dim_ind) const {
        return spline_.periodicity(dim_ind);
    }

    inline bool uniform(size_type dim_ind) const { return uniform_[dim_ind]; }

    const std::pair<typename spline_type::knot_type,
                    typename spline_type::knot_type>&
    range(size_type dim_ind) const {
        return spline_.range(dim_ind);
    }

    /**
     * @brief Get a ref of underlying spline object
     *
     * @return spline_type&
     */
    const spline_type& spline() const { return spline_; }

    static constexpr size_type get_order() { return order; }

   private:
    spline_type spline_;

    DimArray<coord_type> dx_;
    DimArray<bool> uniform_;

    // auxiliary methods

    template <size_type... di>
    inline DimArray<std::pair<coord_type, size_type>> add_hint_for_spline(
        util::index_sequence<di...>,
        DimArray<coord_type> c) const {
        return {std::make_pair(
            c[di],
            uniform_[di]
                ? std::min(
                      spline_.knots_num(di) - order - 2,
                      static_cast<size_type>(std::ceil(std::max(
                          coord_type{0.},
                          (c[di] - range(di).first) / dx_[di] -
                              (periodicity(di)
                                   ? coord_type{1.}
                                   : coord_type{.5} * static_cast<coord_type>(
                                                          order + 1))))) +
                          order)
                : order)...};
    }

    inline val_type call_op_helper(DimArray<coord_type> c) const {
        return spline_(
            add_hint_for_spline(util::make_index_sequence<dim>{}, c));
    }

    template <size_type... di>
    inline val_type derivative_helper(util::index_sequence<di...>,
                                      DimArray<coord_type> c,
                                      DimArray<size_type> d) const {
        return spline_.derivative_at({std::make_tuple(
            static_cast<coord_type>(c[di]),
            uniform_[di]
                ? std::min(spline_.knots_num(di) - order - 2,
                           static_cast<size_type>(std::ceil(std::max(
                               0., (c[di] - range(di).first) / dx_[di] -
                                       (periodicity(di)
                                            ? 1.
                                            : .5 * static_cast<coord_type>(
                                                       order + 1))))) +
                               order)
                : order,
            static_cast<size_type>(d[di]))...});
    }

    // overload for uniform knots
    template <typename T_>
    typename std::enable_if<std::is_arithmetic<T_>::value>::type
    create_knot_vector_(size_type dim_ind,
                        const MeshDimension<dim>& mesh_dimension,
                        DimArray<typename spline_type::KnotContainer>&,
                        std::pair<T_, T_> x_range) {
        uniform_[dim_ind] = true;
        const auto periodic = periodicity(dim_ind);
        const size_type n = mesh_dimension.dim_size(dim_ind)
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
                            + (periodic ? 1 : 0)
#endif
            ;
        dx_[dim_ind] =
            (x_range.second - x_range.first) / static_cast<coord_type>(n - 1);

        const size_t extra = periodic ? 2 * order + (1 - order % 2) : order + 1;

        std::vector<coord_type> xs(n + extra, x_range.first);

        if (periodic) {
            for (size_type i = 0; i < xs.size(); ++i) {
                xs[i] = x_range.first +
                        (static_cast<coord_type>(i) -
                         coord_type{.5} * static_cast<coord_type>(extra)) *
                            dx_[dim_ind];
            }
        } else {
            for (size_type i = order + 1; i < xs.size() - order - 1; ++i) {
                xs[i] = x_range.first +
                        (static_cast<coord_type>(i) -
                         coord_type{.5} * static_cast<coord_type>(extra)) *
                            dx_[dim_ind];
            }
            for (size_type i = xs.size() - order - 1; i < xs.size(); ++i) {
                xs[i] = x_range.second;
            }
        }

        spline_.load_knots(dim_ind, std::move(xs));
    }

    // overload for nonuniform knots, given by iterator pair
    template <typename T_>
    typename std::enable_if<std::is_convertible<
        typename std::iterator_traits<T_>::iterator_category,
        std::input_iterator_tag>::value>::type
    create_knot_vector_(
        size_type dim_ind,
        const MeshDimension<dim>& mesh_dimension,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        std::pair<T_, T_> x_range) {
        uniform_[dim_ind] = false;
        const auto periodic = periodicity(dim_ind);

        const size_type n{static_cast<size_type>(
            std::distance(x_range.first, x_range.second))};
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
        INTP_ASSERT(n == mesh_dimension.dim_size(dim_ind) + (periodic ? 1 : 0),
                    std::string("Inconsistency between knot number and "
                                "interpolated value number at dimension ") +
                        std::to_string(dim_ind));
#else
        INTP_ASSERT(n == mesh_dimension.dim_size(dim_ind),
                    std::string("Inconsistency between knot number and "
                                "interpolated value number at dimension ") +
                        std::to_string(dim_ind));
#endif
#ifndef INTP_ENABLE_ASSERTION
        // suppress unused parameter warning
        (void)mesh_dimension;
#endif
        typename spline_type::KnotContainer xs(
            periodic ? n + 2 * order + (1 - order % 2) : n + order + 1);

        auto& input_coord = input_coords[dim_ind];
        input_coord.reserve(n);
        // The x_range may be given by input iterators, which can not be
        // multi-passed.
        if (periodic) {
            // In periodic case, the knots are data points, shifted by half
            // of local grid size if spline order is odd.

            auto iter = x_range.first;
            input_coord.push_back(*iter);
            for (size_type i = order + 1; i < order + n; ++i) {
                coord_type present = *(++iter);
                xs[i] = order % 2 == 0 ? .5 * (input_coord.back() + present)
                                       : present;
                input_coord.push_back(present);
            }
            coord_type period = input_coord.back() - input_coord.front();
            for (size_type i = 0; i < order + 1; ++i) {
                xs[i] = xs[n + i - 1] - period;
                xs[xs.size() - i - 1] = xs[xs.size() - i - n] + period;
            }
        } else {
            // In aperiodic case, the internal knots are moving average of
            // data points with windows size equal to spline order.

            auto it = x_range.first;
            auto l_knot = *it;
            // fill leftmost *order+1* identical knots
            for (size_type i = 0; i < order + 1; ++i) { xs[i] = l_knot; }
            // first knot is same as first input coordinate
            input_coord.emplace_back(l_knot);
            // Every knot in middle is average of *order* input
            // coordinates. This var is to track the sum of a moving window
            // with width *order*.
            coord_type window_sum{};
            for (size_type i = 1; i < order; ++i) {
                input_coord.emplace_back(*(++it));
                window_sum += input_coord[i];
            }
            for (size_type i = order + 1; i < n; ++i) {
                input_coord.emplace_back(*(++it));
                window_sum += input_coord[i - 1];
                xs[i] = window_sum / static_cast<coord_type>(order);
                window_sum -= input_coord[i - order];
            }
            auto r_knot = *(++it);
            // fill rightmost *order+1* identical knots
            for (size_type i = n; i < n + order + 1; ++i) { xs[i] = r_knot; }
            // last knot is same as last input coordinate
            input_coord.emplace_back(r_knot);
        }
        // check whether input coordinates is monotonic
        for (std::size_t i = 0; i < input_coord.size() - 1; ++i) {
            INTP_ASSERT(input_coord[i + 1] > input_coord[i],
                        std::string("Given coordinate is not monotonically "
                                    "increasing at dimension ") +
                            std::to_string(dim_ind));
        }

#ifdef INTP_TRACE
        std::cout << "[TRACE] Nonuniform knots along dimension" << dim_ind
                  << ":\n";
        for (auto& c : xs) { std::cout << "[TRACE] " << c << '\n'; }
        std::cout << std::endl;
#endif

        spline_.load_knots(dim_ind, std::move(xs));
    }

    template <typename... Ts, size_type... di>
    void create_knots_(
        util::index_sequence<di...>,
        MeshDimension<dim> mesh_dimension,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        std::pair<Ts, Ts>... x_ranges) {
#if __cplusplus >= 201703L
        (create_knot_vector_(di, mesh_dimension, input_coords, x_ranges), ...);
#else
        // polyfill of C++17 fold expression over comma
        static_cast<void>(std::array<std::nullptr_t, sizeof...(Ts)>{
            (create_knot_vector_(di, mesh_dimension, input_coords, x_ranges),
             nullptr)...});
#endif
    }

    inline void boundary_check_(const DimArray<coord_type>& coord) const {
        for (size_type d = 0; d < dim; ++d) {
            if (!periodicity(d) &&
                (coord[d] < range(d).first || coord[d] > range(d).second)) {
                throw std::domain_error(
                    "Given coordinate out of interpolation function "
                    "range!");
            }
        }
    }

#ifdef INTP_CELL_LAYOUT
#if __cplusplus >= 201402L
    auto
#else
    std::function<val_type(const function_type&)>
#endif
    eval_proxy(DimArray<coord_type> coords) const {
        auto spline_proxy = spline().pre_calc_coef(
            add_hint_for_spline(util::make_index_sequence<dim>{}, coords));
        return [spline_proxy](const function_type& interp) {
            return spline_proxy(interp.spline());
        };
    }
#endif  // INTP_CELL_LAYOUT
};

template <std::size_t O = std::size_t{3},
          typename T = double,
          typename U = double>
class InterpolationFunction1D
    : public InterpolationFunction<T, std::size_t{1}, O, U> {
   private:
    using base = InterpolationFunction<T, size_t{1}, O, U>;

   public:
    template <typename InputIter>
    InterpolationFunction1D(std::pair<InputIter, InputIter> f_range,
                            bool periodicity = false)
        : InterpolationFunction1D(
              std::make_pair(typename base::coord_type{},
                             static_cast<typename base::coord_type>(
                                 f_range.second - f_range.first
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
                                 - (periodicity ? 0 : 1)
#else
                                 - 1
#endif
                                     )),
              f_range,
              periodicity) {
    }

    template <typename C1, typename C2, typename InputIter>
    InterpolationFunction1D(std::pair<C1, C2> x_range,
                            std::pair<InputIter, InputIter> f_range,
                            bool periodicity = false)
        : base(periodicity, f_range, x_range) {}
};

}  // namespace intp

#endif
