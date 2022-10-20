#pragma once

#include <cmath>  // ceil
#include <initializer_list>

#include "InterpolationTemplate.hpp"

namespace intp {

template <typename T, size_t D>
class InterpolationFunction {  // TODO: Add integration

    /**
     * The instantiation of this class is prior to that of
     * InterpolationTemplate, thus the use of InterpolationTemplate in the
     * following context must accept incomplete type (e.g. friend declaration,
     * in a template member function, etc.)
     *
     */

   public:
    using val_type = T;
    using spline_type = BSpline<T, D>;
    using size_type = typename spline_type::size_type;

    const static size_type dim = D;

    using coord_type = typename spline_type::knot_type;
    using diff_type = typename spline_type::diff_type;
    using template_type = InterpolationFunctionTemplate<val_type, dim>;

    template <typename U>
    using DimArray = std::array<U, dim>;

    /**
     * @brief Parameters for an interpolation
     *
     */
    struct InputParameters {
        size_type order = 3;              // default interpolation order: 3
        DimArray<bool> periodicity = {};  // default periodicity: aperiodic
    };

    /**
     * @brief Parameters for an interpolation
     *
     */
    struct Parameters : public InputParameters {
        DimArray<bool> uniform = {};
        Parameters(InputParameters p) : InputParameters(std::move(p)) {}
    };

   public:
    /**
     * @brief Construct a new nD Interpolation Function object, mimicking
     * Mathematica's `Interpolation` function, with option `Method->"Spline"`.
     * Notice: last value of periodic dimension will be discarded since it is
     * considered same of the first value. Thus inconsistency input data for
     * periodic interpolation will be accepted.
     *
     * @param spline_order order of interpolation, the interpolated function is
     * of $C^{order-1}$
     * @param periodicity an array describing periodicity of each dimension
     * @param f_mesh a mesh containing data to be interpolated
     * @param x_ranges pairs of x_min and x_max or begin and end iterator
     */
    template <typename... Ts>
    InterpolationFunction(size_type spline_order,
                          DimArray<bool> periodicity,
                          const Mesh<val_type, dim>& f_mesh,
                          std::pair<Ts, Ts>... x_ranges)
        : InterpolationFunction(template_type{spline_order, periodicity,
                                              f_mesh.dimension(), x_ranges...}
                                    .interpolate(f_mesh)) {}

    // Non-periodic for all dimension
    template <typename... Ts>
    InterpolationFunction(size_type spline_order,
                          const Mesh<val_type, dim>& f_mesh,
                          std::pair<Ts, Ts>... x_ranges)
        : InterpolationFunction(spline_order, {}, f_mesh, x_ranges...) {}

    // constructor for partial construction, that is, without interpolated
    // values
    template <typename... Ts>
    InterpolationFunction(
        size_type spline_order,
        DimArray<bool> periodicity,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        MeshDimension<dim> mesh_dimension,
        std::pair<Ts, Ts>... x_ranges)
        : parameters_{{spline_order, periodicity}},
          spline_(periodicity, spline_order) {
        // load knots into spline
        create_knots_(util::index_sequence_for<Ts...>{},
                      std::move(mesh_dimension), input_coords, x_ranges...);
    }

    // New constructors

    // Constructor for partial construction, that is, without interpolated
    // values. This overload can only be used in uniform case.
    InterpolationFunction(
        MeshDimension<dim> mesh_dimension,
        util::n_pairs_t<coord_type, dim> xs_ranges,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        InputParameters parameters)
        : parameters_{parameters},
          spline_(parameters.periodicity, parameters.order) {
        // load knots into spline
        create_knots_(util::make_index_sequence<dim>{},
                      std::move(mesh_dimension), input_coords,
                      std::move(xs_ranges));
    }

    // Constructor for partial construction, that is, without interpolated
    // values.
    template <typename... Ts>
    InterpolationFunction(
        MeshDimension<dim> mesh_dimensions,
        std::tuple<std::pair<Ts, Ts>...> xs_ranges,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        InputParameters parameters)
        : parameters_{parameters},
          spline_(parameters.periodicity, parameters.order) {
        // load knots into spline
        create_knots_(util::index_sequence_for<Ts...>{},
                      std::move(mesh_dimensions), input_coords,
                      std::move(xs_ranges));
    }

    /**
     * @brief Construct a new nD Interpolation Function object, mimicking
     * Mathematica's `Interpolation` function, with option `Method->"Spline"`
     * and default interpolation order is 3. Notice: last value of periodic
     * dimension will be discarded since it is considered same of the first
     * value. Thus inconsistency input data for periodic interpolation will be
     * accepted.
     *
     * @param f_mesh a mesh containing data to be interpolated
     * @param xs_ranges a tuple of pairs of x_min and x_max or begin and end
     * iterator
     * @param parameters specifying order and periodicity
     */
    InterpolationFunction(const Mesh<val_type, dim>& f_mesh,
                          util::n_pairs_t<coord_type, dim> xs_ranges,
                          InputParameters parameters = {})
        : InterpolationFunction(
              template_type(f_mesh.dimension(), xs_ranges, parameters)
                  .interpolate(f_mesh)) {}

    template <typename... Ts>
    InterpolationFunction(const Mesh<val_type, dim>& f_mesh,
                          std::tuple<std::pair<Ts, Ts>...> xs_ranges,
                          InputParameters parameters = {})
        : InterpolationFunction(
              template_type(f_mesh.dimension(), xs_ranges, parameters)
                  .interpolate(f_mesh)) {}

    /**
     * @brief Get spline value.
     *
     * @param x coordinates
     */
    template <typename... Coords,
              typename Indices = util::index_sequence_for<Coords...>,
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
        boundary_check_(coord);
        return call_op_helper(util::make_index_sequence<dim>{}, coord);
    }

    /**
     * @brief Get spline value, but with out of boundary check.
     *
     * @param x coordinates
     */
    template <typename... Coords,
              typename Indices = util::index_sequence_for<Coords...>,
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
        return parameters_.periodicity[dim_ind];
    }

    inline bool uniform(size_type dim_ind) const {
        return parameters_.uniform[dim_ind];
    }

    inline size_type order() const { return parameters_.order; }

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

   private:
    Parameters parameters_;
    spline_type spline_;
    DimArray<coord_type> dx_;  // use in uniform case

    friend template_type;

    // auxiliary methods

    template <size_type... di>
    inline val_type call_op_helper(util::index_sequence<di...>,
                                   DimArray<coord_type> c) const {
        return spline_(std::make_pair(
            c[di],
            uniform(di)
                ? std::min(spline_.knots_num(di) - order() - 2,
                           static_cast<size_type>(std::ceil(std::max(
                               0., (c[di] - range(di).first) / dx_[di] -
                                       (periodicity(di)
                                            ? 1.
                                            : .5 * static_cast<coord_type>(
                                                       order() + 1))))) +
                               order())
                : order())...);
    }

    template <size_type... di>
    inline val_type derivative_helper(util::index_sequence<di...>,
                                      DimArray<coord_type> c,
                                      DimArray<size_type> d) const {
        return spline_.derivative_at(std::make_tuple(
            static_cast<coord_type>(c[di]), static_cast<size_type>(d[di]),
            uniform(di)
                ? std::min(spline_.knots_num(di) - order() - 2,
                           static_cast<size_type>(std::ceil(std::max(
                               0., (c[di] - range(di).first) / dx_[di] -
                                       (periodicity(di)
                                            ? 1.
                                            : .5 * static_cast<coord_type>(
                                                       order() + 1))))) +
                               order())
                : order())...);
    }

    // overload for uniform knots
    template <typename T_>
    typename std::enable_if<std::is_arithmetic<T_>::value>::type
    create_knot_vector_(size_type dim_ind,
                        const MeshDimension<dim>& mesh_dimension,
                        DimArray<typename spline_type::KnotContainer>&,
                        std::pair<T_, T_> x_range) {
        parameters_.uniform[dim_ind] = true;
        const size_type n = mesh_dimension.dim_size(dim_ind);
        dx_[dim_ind] =
            (x_range.second - x_range.first) / static_cast<coord_type>(n - 1);

        const size_t extra = periodicity(dim_ind)
                                 ? 2 * order() + (1 - order() % 2)
                                 : order() + 1;

        std::vector<typename spline_type::knot_type> xs(n + extra,
                                                        x_range.first);

        if (periodicity(dim_ind)) {
            for (size_type i = 0; i < xs.size(); ++i) {
                xs[i] = x_range.first + (static_cast<coord_type>(i) -
                                         .5 * static_cast<coord_type>(extra)) *
                                            dx_[dim_ind];
            }
        } else {
            for (size_type i = order() + 1; i < xs.size() - order() - 1; ++i) {
                xs[i] = x_range.first + (static_cast<coord_type>(i) -
                                         .5 * static_cast<coord_type>(extra)) *
                                            dx_[dim_ind];
            }
            for (size_type i = xs.size() - order() - 1; i < xs.size(); ++i) {
                xs[i] = x_range.second;
            }
        }

        spline_.load_knots(dim_ind, std::move(xs), periodicity(dim_ind));
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
        parameters_.uniform[dim_ind] = false;
        const size_type n{static_cast<size_type>(
            std::distance(x_range.first, x_range.second))};
        if (n != mesh_dimension.dim_size(dim_ind)) {
            throw std::range_error(
                std::string("Inconsistency between knot number and "
                            "interpolated value number at dimension ") +
                std::to_string(dim_ind));
        }
        typename spline_type::KnotContainer xs(
            periodicity(dim_ind) ? n + 2 * order() + (1 - order() % 2)
                                 : n + order() + 1);

        auto& input_coord = input_coords[dim_ind];
        input_coord.reserve(n);
        // The x_range may be given by input iterators, which can not be
        // multi-passed.
        if (periodicity(dim_ind)) {
            // In periodic case, the knots are data points, shifted by half of
            // local grid size if spline order is odd.

            auto iter = x_range.first;
            input_coord.push_back(*iter);
            for (size_type i = order() + 1; i < order() + n; ++i) {
                val_type present = *(++iter);
                xs[i] = order() % 2 == 0 ? .5 * (input_coord.back() + present)
                                         : present;
                input_coord.push_back(present);
            }
            val_type period = input_coord.back() - input_coord.front();
            for (size_type i = 0; i < order() + 1; ++i) {
                xs[i] = xs[n + i - 1] - period;
                xs[xs.size() - i - 1] = xs[xs.size() - i - n] + period;
            }
        } else {
            // In aperiodic case, the internal knots are moving average of data
            // points with windows size equal to spline order.

            auto it = x_range.first;
            auto l_knot = *it;
            // fill leftmost *order+1* identical knots
            for (size_type i = 0; i < order() + 1; ++i) { xs[i] = l_knot; }
            // first knot is same as first input coordinate
            input_coord.emplace_back(l_knot);
            // Every knot in middle is average of *order* input
            // coordinates. This var is to track the sum of a moving window with
            // width *order*.
            coord_type window_sum{};
            for (size_type i = 1; i < order(); ++i) {
                input_coord.emplace_back(*(++it));
                window_sum += input_coord[i];
            }
            for (size_type i = order() + 1; i < n; ++i) {
                input_coord.emplace_back(*(++it));
                window_sum += input_coord[i - 1];
                xs[i] = window_sum / static_cast<coord_type>(order());
                window_sum -= input_coord[i - order()];
            }
            auto r_knot = *(++it);
            // fill rightmost *order+1* identical knots
            for (size_type i = n; i < n + order() + 1; ++i) { xs[i] = r_knot; }
            // last knot is same as last input coordinate
            input_coord.emplace_back(r_knot);
        }
#ifdef _DEBUG
        // check whether input coordinates is monotonic
        for (std::size_t i = 0; i < input_coord.size() - 1; ++i) {
            if (input_coord[i + 1] <= input_coord[i]) {
                throw std::range_error(
                    std::string("Given coordinate is not monotonically "
                                "increasing at dimension ") +
                    std::to_string(dim_ind));
            }
        }
#endif
#ifdef _TRACE
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

    template <typename... Ts, size_type... di>
    void create_knots_(
        util::index_sequence<di...> indices,
        MeshDimension<dim> mesh_dimension,
        DimArray<typename spline_type::KnotContainer>& input_coords,
        std::tuple<std::pair<Ts, Ts>...> xs_ranges) {
        create_knots_(indices, std::move(mesh_dimension), input_coords,
                      std::get<di>(xs_ranges)...);
    }

    inline void boundary_check_(const DimArray<coord_type>& coord) const {
        for (size_type d = 0; d < dim; ++d) {
            if (!periodicity(d) &&
                (coord[d] < range(d).first || coord[d] > range(d).second)) {
                throw std::domain_error(
                    "Given coordinate out of interpolation function range!");
            }
        }
    }
};

template <typename T = double>
class InterpolationFunction1D : public InterpolationFunction<T, size_t{1}> {
   private:
    using base = InterpolationFunction<T, size_t{1}>;

   public:
    template <typename C1, typename C2, typename InputIter>
    InterpolationFunction1D(std::pair<InputIter, InputIter> f_range,
                            std::pair<C1, C2> x_range,
                            typename base::InputParameters parameter = {})
        : base(Mesh<typename base::val_type, 1u>{std::move(f_range)},
               std::make_tuple(x_range),
               parameter) {}

    template <typename InputIter>
    InterpolationFunction1D(
        std::pair<InputIter, InputIter> f_range,
        std::pair<typename base::coord_type, typename base::coord_type> x_range,
        typename base::InputParameters parameter = {})
        : base(Mesh<typename base::val_type, 1u>{std::move(f_range)},
               std::make_tuple(x_range),
               parameter) {}
};

}  // namespace intp
