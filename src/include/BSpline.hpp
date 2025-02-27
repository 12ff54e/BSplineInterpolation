#ifndef INTP_BSPLINE
#define INTP_BSPLINE

#include <algorithm>  // upper_bound
#include <array>
#include <cmath>        // fmod
#include <functional>   // ref
#include <iterator>     // distance
#include <type_traits>  // is_same, is_arithmatic
#include <vector>

#ifdef INTP_MULTITHREAD
#include "DedicatedThreadPool.hpp"
#endif

#ifdef INTP_DEBUG
#include <iostream>
#endif

#include "Mesh.hpp"
#include "aligned-allocator.hpp"
#include "util.hpp"

namespace intp {

/**
 * @brief B-Spline function
 *
 * @tparam T Type of control point
 * @tparam D Dimension
 */
template <typename T, std::size_t D, std::size_t O, typename U = double>
class BSpline {
   public:
    using spline_type = BSpline<T, D, O, U>;

    using size_type = std::size_t;
    using val_type = T;
    using knot_type = U;
    constexpr static size_type dim = D;
    constexpr static size_type order = O;

    using KnotContainer = std::vector<knot_type>;

    using ControlPointContainer =
        Mesh<val_type,
             dim,
             util::default_init_allocator<
                 val_type,
                 AlignedAllocator<val_type, Alignment::AVX>>>;
#ifdef INTP_CELL_LAYOUT
    using ControlPointCellContainer =
        Mesh<val_type,
             dim + 1,
             util::default_init_allocator<
                 val_type,
                 AlignedAllocator<val_type, Alignment::AVX>>>;
    using control_point_type = ControlPointCellContainer;
#else
    using control_point_type = ControlPointContainer;
#endif

    using BaseSpline = std::array<knot_type, order + 1>;
    using diff_type = typename KnotContainer::iterator::difference_type;
    using knot_const_iterator = typename KnotContainer::const_iterator;

    // Container for dimension-wise storage
    template <typename T_>
    using DimArray = std::array<T_, dim>;

    /**
     * @brief Calculate values on base spline function. This is the core of
     * B-Spline. Note: when the given order is smaller than order of spline
     * (used in calculating derivative), spline value is aligned at right in
     * result vector.
     *
     * @param seg_idx_iter the iterator points to left knot point of a segment
     * @param x coordinate
     * @param spline_order order of base spline, defaulted to be spline function
     * order
     * @return a reference to local buffer
     */
    inline const BaseSpline base_spline_value(
        size_type,
        knot_const_iterator seg_idx_iter,
        knot_type x,
        size_type spline_order = order) const {
        BaseSpline base_spline{};
        base_spline[order] = 1;

        for (size_type i = 1; i <= spline_order; ++i) {
            // Each iteration will expand buffer zone by one, from back
            // to front.
            const size_type idx_begin = order - i;
            for (size_type j = 0; j <= i; ++j) {
                const auto left_iter =
                    seg_idx_iter - static_cast<diff_type>(i - j);
                const auto right_iter =
                    seg_idx_iter + static_cast<diff_type>(j + 1);
                base_spline[idx_begin + j] =
                    (j == 0 ? 0
                            : base_spline[idx_begin + j] * (x - *left_iter) /
                                  (*(right_iter - 1) - *left_iter)) +
                    (idx_begin + j == order
                         ? 0
                         : base_spline[idx_begin + j + 1] * (*right_iter - x) /
                               (*right_iter - *(left_iter + 1)));
            }
        }
        return base_spline;
    }

    /**
     * @brief Get the lower knot iter points to the segment where given x
     * locates. If x is out of range of knot vector, the iterator is rather
     * begin or end of knot vector.
     *
     * @param dim_ind specify the dimension
     * @param x coordinate
     * @param hint a hint for iter offset
     * @param last an upper bound for iter offset, this function will not search
     * knots beyond it.
     * @return knot_const_iterator
     */
    inline knot_const_iterator get_knot_iter(size_type dim_ind,
                                             knot_type& x,
                                             size_type hint,
                                             size_type last) const {
        const auto iter = knots_begin(dim_ind) + static_cast<diff_type>(hint);
        if (periodicity_[dim_ind]) {
            const knot_type period =
                range(dim_ind).second - range(dim_ind).first;
            x = range(dim_ind).first +
                std::fmod(x - range(dim_ind).first, period) +
                (x < range(dim_ind).first ? period : knot_type{});
        }
#ifdef INTP_TRACE
        if ((*iter > x || *(iter + 1) < x) && x >= range(dim_ind).first &&
            x <= range(dim_ind).second) {
            std::cout << "[TRACE] knot hint miss at dim = " << dim_ind
                      << ", hint = " << hint << ", x = " << x << '\n';
        }
#endif
        // I tried return the iter without checking, but the speed has no
        // significant improves.
        return *iter <= x && *(iter + 1) > x
                   // If the hint is accurate, use that iter
                   ? iter
                   // else, use binary search in the range of distinct knots
                   // (excluding beginning and ending knots that have same
                   // value)
                   : --(std::upper_bound(knots_begin(dim_ind) +
                                             static_cast<diff_type>(order + 1),
                                         knots_begin(dim_ind) +
                                             static_cast<diff_type>(last + 1),
                                         x));
    }

    inline knot_const_iterator get_knot_iter(size_type dim_ind,
                                             knot_type& x,
                                             size_type hint) const {
        return get_knot_iter(dim_ind, x, hint, knots_num(dim_ind) - order - 2);
    }

    template <typename C, size_type... indices>
    inline DimArray<knot_const_iterator> get_knot_iters(
        util::index_sequence<indices...>,
        C&& coords) const {
        return {get_knot_iter(indices, std::get<0>(coords[indices]),
                              std::get<1>(coords[indices]))...};
    }

    /**
     * @brief Construct a new BSpline object, with periodicity of each dimension
     * specified.
     *
     */
    explicit BSpline(DimArray<bool> periodicity)
        : periodicity_(periodicity), control_points_(size_type{}) {}

    /**
     * @brief Basically the default constructor, initialize an empty, non-closed
     * B-Spline
     *
     */
    explicit BSpline() : BSpline(DimArray<bool>{}) {}

    template <typename... InputIters>
    BSpline(DimArray<bool> periodicity,
            ControlPointContainer ctrl_pts,
            std::pair<InputIters, InputIters>... knot_iter_pairs)
        : periodicity_(periodicity),
          knots_{
              KnotContainer(knot_iter_pairs.first, knot_iter_pairs.second)...},
#ifdef INTP_CELL_LAYOUT
          control_points_(generate_cell_layout(ctrl_pts)),
#else
          control_points_(std::move(ctrl_pts)),
#endif
          range_{std::make_pair(
              (knot_iter_pairs.first)[order],
              (knot_iter_pairs.second)[-static_cast<int>(order) - 1])...} {
        for (size_type d = 0; d < dim; ++d) {
            INTP_ASSERT(knots_[d].size() - ctrl_pts.dim_size(d) ==
                            (periodicity_[d] ? 2 * order + 1 : order + 1),
                        std::string("Inconsistency between knot number and "
                                    "control point number at dimension ") +
                            std::to_string(d));
        }
    }

    template <typename... InputIters>
    BSpline(ControlPointContainer ctrl_points,
            std::pair<InputIters, InputIters>... knot_iter_pairs)
        : BSpline(DimArray<bool>{}, ctrl_points, knot_iter_pairs...) {}

    template <typename C>
    typename std::enable_if<
        std::is_same<typename std::remove_reference<C>::type,
                     KnotContainer>::value,
        void>::type
    load_knots(size_type dim_ind, C&& _knots) {
        knots_[dim_ind] = std::forward<C>(_knots);
        range_[dim_ind].first = knots_[dim_ind][order];
        range_[dim_ind].second =
            knots_[dim_ind][knots_[dim_ind].size() - order - (2 - order % 2)];
    }

#ifdef INTP_CELL_LAYOUT
    void load_ctrlPts(const ControlPointContainer& control_points) {
        control_points_ = generate_cell_layout(control_points);
    }
#else
    template <typename C>
    typename std::enable_if<
        std::is_same<typename std::remove_reference<C>::type,
                     ControlPointContainer>::value,
        void>::type
    load_ctrlPts(C&& control_points) {
        control_points_ = std::forward<C>(control_points);
    }
#endif

#ifdef INTP_CELL_LAYOUT
#if __cplusplus >= 201402L
    auto
#else
    std::function<val_type(const spline_type&)>
#endif
    pre_calc_coef(
        DimArray<std::pair<knot_type, size_type>> coord_with_hints) const {
        using Indices = util::make_index_sequence<dim>;
        // get knot point iter, it will modifies coordinate value into
        // interpolation range of periodic dimension.
        const auto knot_iters = get_knot_iters(Indices{}, coord_with_hints);

        DimArray<size_type> spline_order;
        spline_order.fill(order);
        // calculate basic spline (out of boundary check also conducted here)
        const auto base_spline_values_1d = calc_base_spline_vals(
            Indices{}, knot_iters, spline_order, coord_with_hints);

        std::array<size_type, dim + 1> ind_arr{};
        for (size_type d = 0; d < dim; ++d) {
            ind_arr[d] = static_cast<size_type>(
                             distance(knots_begin(d), knot_iters[d])) -
                         order;
        }

        auto total_offset = calculate_cell_dim_from_knots().indexing(ind_arr);

        return
            [base_spline_values_1d, total_offset](const spline_type& spline) {
                val_type v{};

                const auto& control_points = spline.control_points();
                auto cell_iter = control_points.begin() +
                                 static_cast<std::ptrdiff_t>(total_offset);
                for (size_type i = 0;
                     i < control_points.dim_size(dim) * (order + 1); ++i) {
                    knot_type coef = 1;
                    if CPP17_CONSTEXPR_ (dim == 1) {
                        // helps with vectorization in 1D case
                        coef = base_spline_values_1d[0][i];
                    } else {
                        for (size_type d = 0, combined_ind = i; d < dim; ++d) {
                            coef *= base_spline_values_1d[d][combined_ind %
                                                             (order + 1)];
                            combined_ind /= (order + 1);
                        }
                    }
                    v += coef * (*cell_iter++);
                }
                return v;
            };
    }
#endif

    /**
     * @brief Get spline value at given pairs of coordinate and position hint
     * (hopefully lower knot point index of the segment where coordinate
     * locates, dimension wise).
     *
     */
    val_type operator()(
        DimArray<std::pair<knot_type, size_type>> coord_with_hints) const {
        using Indices = util::make_index_sequence<dim>;
        // get knot point iter, it will modifies coordinate value into
        // interpolation range of periodic dimension.
        const auto knot_iters = get_knot_iters(Indices{}, coord_with_hints);

        DimArray<size_type> spline_order;
        spline_order.fill(order);
        // calculate basic spline (out of boundary check also conducted here)
        const auto base_spline_values_1d = calc_base_spline_vals(
            Indices{}, knot_iters, spline_order, coord_with_hints);

        // combine control points and basic spline values to get spline value
        val_type v{};
#ifdef INTP_CELL_LAYOUT
        std::array<size_type, dim + 1> ind_arr{};
        for (size_type d = 0; d < dim; ++d) {
            ind_arr[d] = static_cast<size_type>(
                             distance(knots_begin(d), knot_iters[d])) -
                         order;
        }
        auto cell_iter = control_points_.begin(dim, ind_arr);
        for (size_type i = 0; i < buf_size_; ++i) {
            auto coef = *cell_iter++;
            for (size_type d = 0, combined_ind = i; d < dim; ++d) {
                coef *= base_spline_values_1d[d][combined_ind % (order + 1)];
                combined_ind /= (order + 1);
            }
            v += coef;
        }
#else
        MeshDimension<dim> local_mesh_dim(order + 1);
        for (size_type i = 0; i < buf_size_; ++i) {
            DimArray<size_type> ind_arr = local_mesh_dim.dimwise_indices(i);

            knot_type coef = 1;
            for (size_type d = 0; d < dim; ++d) {
                coef *= base_spline_values_1d[d][ind_arr[d]];

                // Shift index array according to knot iter of each dimension.
                // When the coordinate is out of range in some dimensions, the
                // corresponding iterator was set to be begin or end iterator of
                // knot vector in `get_knot_iters` method and it will be treated
                // separately.
                ind_arr[d] += knot_iters[d] == knots_begin(d) ? 0
                              : knot_iters[d] == knots_end(d)
                                  ? control_points_.dim_size(d) - order - 1
                                  : static_cast<size_type>(distance(
                                        knots_begin(d), knot_iters[d])) -
                                        order;

                // check periodicity, put out-of-right-boundary index to left
                if (periodicity_[d]) {
                    ind_arr[d] %= control_points_.dim_size(d);
                }
            }

            v += coef * control_points_(ind_arr);
        }
#endif

        return v;
    }

    /**
     * @brief Get spline value at given coordinates
     *
     * @param coords a bunch of cartesian coordinates
     * @return val_type
     */
    val_type operator()(DimArray<double> coords) const {
        DimArray<std::pair<knot_type, size_type>> coord_with_hints;
        for (std::size_t d = 0; d < dim; ++d) {
            coord_with_hints[d] = {coords[d], order};
        }
        return operator()(coord_with_hints);
    }

    /**
     * @brief Get derivative value at given pairs of coordinate and position
     * hint (possibly lower knot point index of the segment where coordinate
     * locates, dimension wise).
     *
     * @param coord_deriOrder_hint_tuple a bunch of (coordinate, position hint,
     * derivative order) tuple
     * @return val_type
     */
    val_type derivative_at(DimArray<std::tuple<knot_type, size_type, size_type>>
                               coord_deriOrder_hint_tuple) const {
        // get spline order
        DimArray<size_type> spline_order;
        for (size_type d = 0; d < dim; ++d) {
            spline_order[d] =
                order >= std::get<2>(coord_deriOrder_hint_tuple[d])
                    ? order - std::get<2>(coord_deriOrder_hint_tuple[d])
                    : order + 1;
        }

        // if derivative order is larger than spline order, derivative is 0.
        for (auto o : spline_order) {
            if (o > order) { return val_type{}; }
        }

        using Indices = util::make_index_sequence<dim>;
        // get knot point iter (out of boundary check also conducted here)
        const auto knot_iters =
            get_knot_iters(Indices{}, coord_deriOrder_hint_tuple);

        // calculate basic spline
        const auto base_spline_values_1d = calc_base_spline_vals(
            Indices{}, knot_iters, spline_order, coord_deriOrder_hint_tuple);

#ifdef STACK_ALLOCATOR
        // create local buffer
        val_type buffer[MAX_BUF_SIZE_];
        util::stack_allocator<val_type, MAX_BUF_SIZE_> alloc(buffer);

        Mesh<val_type, dim, util::stack_allocator<val_type, MAX_BUF_SIZE_>>
            local_control_points(order + 1, alloc);
        auto local_spline_val = local_control_points;
#else
        Mesh<val_type, dim> local_control_points(order + 1);
        auto local_spline_val = local_control_points;
#endif

        // get local control points and basic spline values

#ifdef INTP_CELL_LAYOUT
        std::array<size_type, dim + 1> ind_arr{};
        for (size_type d = 0; d < dim; ++d) {
            ind_arr[d] = static_cast<size_type>(
                             distance(knots_begin(d), knot_iters[d])) -
                         order;
        }
        auto cell_iter = control_points_.begin(dim, ind_arr);
        for (size_type i = 0; i < buf_size_; ++i) {
            DimArray<size_type> local_ind_arr{};
            for (size_type d = 0, combined_ind = i; d < dim; ++d) {
                local_ind_arr[d] = combined_ind % (order + 1);
                combined_ind /= (order + 1);
            }

            knot_type coef = 1;
            for (size_type d = 0; d < dim; ++d) {
                coef *= base_spline_values_1d[d][local_ind_arr[d]];
            }

            local_spline_val(local_ind_arr) = coef;
            local_control_points(local_ind_arr) = *cell_iter++;
        }
#else
        for (size_type i = 0; i < buf_size_; ++i) {
            DimArray<size_type> local_ind_arr{};
            for (size_type d = 0, combined_ind = i; d < dim; ++d) {
                local_ind_arr[d] = combined_ind % (order + 1);
                combined_ind /= (order + 1);
            }

            knot_type coef = 1;
            DimArray<size_type> ind_arr{};
            for (size_type d = 0; d < dim; ++d) {
                coef *= base_spline_values_1d[d][local_ind_arr[d]];

                ind_arr[d] = local_ind_arr[d] +
                             (knot_iters[d] == knots_begin(d) ? 0
                              : knot_iters[d] == knots_end(d)
                                  ? control_points_.dim_size(d) - order - 1
                                  : static_cast<size_t>(distance(
                                        knots_begin(d), knot_iters[d])) -
                                        order);

                // check periodicity, put out-of-right-boundary index to
                // left
                if (periodicity_[d]) {
                    ind_arr[d] %= control_points_.dim_size(d);
                }
            }

            local_spline_val(local_ind_arr) = coef;
            local_control_points(local_ind_arr) = control_points_(ind_arr);
        }
#endif

        for (size_type d = 0; d < dim; ++d) {
            if (spline_order[d] == order) { continue; }
            // calculate control points for derivative along this dimension

            const size_type hyper_surface_size =
                local_control_points.size() / local_control_points.dim_size(d);
            // transverse the hyper surface of fixing dimension d
            for (size_type i = 0; i < hyper_surface_size; ++i) {
                DimArray<size_type> local_ind_arr{};
                for (size_type dd = 0, combined_ind = i; dd < dim; ++dd) {
                    if (dd == d) { continue; }
                    local_ind_arr[dd] = combined_ind % (order + 1);
                    combined_ind /= (order + 1);
                }

                auto iter = local_control_points.begin(d, local_ind_arr);
                // Taking derivative is effectively computing new control
                // points. Number of iteration is order of derivative.
                for (diff_type k = static_cast<diff_type>(order);
                     k > static_cast<diff_type>(spline_order[d]); --k) {
                    // Each reduction reduce control points number by one.
                    // Reduce backward to match pattern of local_spline_val.
                    for (diff_type j = k; j > 0; --j) {
                        iter[static_cast<diff_type>(order) + j - k] =
                            static_cast<val_type>(k) *
                            (iter[static_cast<diff_type>(order) + j - k] -
                             iter[static_cast<diff_type>(order) + j - k - 1]) /
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
     * @param coord_deriOrders a bunch of (coordinate, derivative order) tuple
     * @return val_type
     */
    val_type derivative_at(
        DimArray<std::pair<knot_type, size_type>> coord_deriOrders) const {
        DimArray<std::tuple<knot_type, size_type, size_type>>
            coord_deriOrder_hint_tuple;
        for (std::size_t d = 0; d < dim; ++d) {
            coord_deriOrder_hint_tuple[d] = {std::get<0>(coord_deriOrders[d]),
                                             order,
                                             std::get<1>(coord_deriOrders[d])};
        }
        return derivative_at(coord_deriOrder_hint_tuple);
    }

    // iterators

    /**
     * @brief Returns a read-only (constant) iterator that points to the
     *  first element in the knot vector of one dimension.
     *
     * @param dim_ind dimension index
     */
    inline knot_const_iterator knots_begin(size_type dim_ind) const {
        return knots_[dim_ind].cbegin();
    }
    /**
     * @brief Returns a read-only (constant) iterator that points to the
     *  first element in the knot vector of one dimension.
     *
     * @param dim_ind dimension index
     */
    inline knot_const_iterator knots_end(size_type dim_ind) const {
        return knots_[dim_ind].cend();
    }

    // properties

    inline const control_point_type& control_points() const {
        return control_points_;
    }

    /**
     * @brief Get range of one dimension
     *
     */
    inline const std::pair<knot_type, knot_type>& range(
        size_type dim_ind) const {
        return range_[dim_ind];
    }

    /**
     * @brief Get knot number of one dimension
     *
     */
    inline size_type knots_num(size_type dim_ind) const {
        return knots_[dim_ind].size();
    }

    /**
     * @brief Get periodicity of one dimension
     *
     * @param dim_ind dimension index
     * @return a bool
     */
    inline bool periodicity(size_type dim_ind) const {
        return periodicity_[dim_ind];
    }

    inline constexpr size_type get_order() const {
        return order;
    }

#ifdef INTP_DEBUG
    void debug_output() const {
        std::cout << "\n[DEBUG] Control Points (raw data):\n";

        // 17 digits for double precision
        std::cout.precision(17);
        size_type idx = 1;
        for (auto v : control_points_) {
            if (idx % control_points_.dim_size(dim - 1) == 1) {
                std::cout << "[DEBUG] ";
            }
            std::cout << v << ' ';
            if (idx++ % control_points_.dim_size(dim - 1) == 0) {
                std::cout << '\n';
            }
        }
        std::cout << '\n';
    }
#endif

   private:
    DimArray<bool> periodicity_;

    DimArray<KnotContainer> knots_;
#ifdef INTP_CELL_LAYOUT
    ControlPointCellContainer control_points_;
#else
    ControlPointContainer control_points_;
#endif

    DimArray<std::pair<knot_type, knot_type>> range_;

    constexpr static size_type buf_size_ = util::pow(order + 1, dim);

    // maximum stack buffer size
    // This buffer is for storing weights when calculating spline derivative
    // value.
    constexpr static size_type MAX_BUF_SIZE_ = 1000;

    // auxiliary methods

    /**
     * @brief Calculate base spline value of each dimension
     *
     */
    template <typename C, size_type... indices>
    inline DimArray<BaseSpline> calc_base_spline_vals(
        util::index_sequence<indices...>,
        const DimArray<knot_const_iterator>& knot_iters,
        const DimArray<size_type>& spline_order,
        const C& coords) const {
        return {base_spline_value(indices, knot_iters[indices],
                                  std::get<0>(coords[indices]),
                                  spline_order[indices])...};
    }

#ifdef INTP_CELL_LAYOUT
    MeshDimension<dim + 1> calculate_cell_dim_from_knots() const {
        MeshDimension<dim + 1> cell_dim(util::pow(order + 1, dim - 1));
        for (size_type d = 0; d < dim; ++d) {
            cell_dim.dim_size(d) = knots_[d].size() -
                                   (periodicity(d) ? order + 1 : order + 1) -
                                   (d == dim - 1 ? 0 : order);
        }
        return cell_dim;
    }

    ControlPointCellContainer generate_cell_layout(
        const ControlPointContainer& ctrl_pts) const {
        MeshDimension<dim + 1> cell_container_dim(
            util::pow(order + 1, dim - 1));
        for (size_type d = 0; d < dim; ++d) {
            cell_container_dim.dim_size(d) = ctrl_pts.dim_size(d) -
                                             (d == dim - 1 ? 0 : order) +
                                             (periodicity(d) ? order : 0);
        }

        // Size of the last two dimension of control point cell container. The
        // (dim-1)th dimention of cell container has the same length (order
        // points more in periodic case) as the last dimension (which is also
        // the (dim-1)th dimension) of the origin container, while other
        // dimensions are #order shorter, except the last dimension with length
        // (order+1)^(dim-1).
        const auto line_size = cell_container_dim.dim_size(dim - 1) *
                               cell_container_dim.dim_size(dim);

        ControlPointCellContainer control_point_cell(cell_container_dim);
        // size of hyperplane orthogonal to last 2 dimension
        const size_type hyper_surface_size =
            control_point_cell.size() / line_size;

        auto fill_cell = [&](size_type begin, size_type end) {
            // iterate over hyperplane
            for (size_type h_ind = begin; h_ind < end; ++h_ind) {
                const auto ind_arr_on_hyper_surface =
                    control_point_cell.dimension().dimwise_indices(h_ind *
                                                                   line_size);
                const auto line_begin =
                    control_point_cell.begin(dim - 1, ind_arr_on_hyper_surface);
                const auto line_end =
                    control_point_cell.end(dim - 1, ind_arr_on_hyper_surface);
                // iterate along the (dim-1)th dimension
                for (auto iter = line_begin; iter != line_end; ++iter) {
                    auto cell_ind_arr = control_point_cell.iter_indices(iter);
                    MeshDimension<dim - 1> cell_dim(order + 1);
                    // iterate the (dim)th dimension
                    for (size_type i = 0; i < cell_dim.size(); ++i) {
                        auto local_shift_ind = cell_dim.dimwise_indices(i);
                        cell_ind_arr[dim] = i;
                        auto cp_ind_arr =
                            typename ControlPointContainer::index_type{};
                        for (size_type d = 0; d < dim; ++d) {
                            cp_ind_arr[d] =
                                (cell_ind_arr[d] +
                                 (d == dim - 1
                                      ? 0
                                      : local_shift_ind[dim - 2 - d])) %
                                ctrl_pts.dim_size(d);
                        }
                        control_point_cell(cell_ind_arr) = ctrl_pts(cp_ind_arr);
                    }
                }
            }
        };
#ifdef INTP_MULTITHREAD
        const size_type block_num = static_cast<size_type>(
            std::sqrt(static_cast<double>(hyper_surface_size)));
        const size_type task_per_block = block_num == 0
                                             ? hyper_surface_size
                                             : hyper_surface_size / block_num;

        auto& thread_pool = DedicatedThreadPool<void>::get_instance(8);
        std::vector<std::future<void>> res;

        for (size_type i = 0; i < block_num; ++i) {
            res.push_back(thread_pool.queue_task([=]() {
                fill_cell(i * task_per_block, (i + 1) * task_per_block);
            }));
        }
        // main thread deals with the remaining part in case hyper_surface_size
        // not divisible by thread_num
        fill_cell(block_num * task_per_block, hyper_surface_size);
        // wait for all tasks are complete
        for (auto&& f : res) { f.get(); }
#else
        fill_cell(0, hyper_surface_size);
#endif  // INTP_MULTITHREAD
        return control_point_cell;
    }
#endif  // INTP_CELL_LAYOUT
};

}  // namespace intp

#endif
