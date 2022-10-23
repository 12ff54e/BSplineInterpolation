#pragma once

#include <type_traits>

#include "BandMatrix.hpp"
#include "util.hpp"

namespace intp {

/**
 * @brief Base class of band matrix LU solver. The actual solver should inherits
 * from it and implement corresponding impl methods, per the CRTP.
 */
template <template <typename> class Solver, typename Matrix>
class BandLUBase : public util::CRTP<Solver<Matrix>> {
   public:
    using matrix_type = Matrix;

    BandLUBase() noexcept : is_computed_(false) {}

    template <typename Mat_>
    BandLUBase(Mat_&& mat) : lu_store_(std::forward<Mat_>(mat)) {
        this->cast().compute_impl();
        is_computed_ = true;
    }

    template <typename Mat_>
    void compute(Mat_&& mat) {
        static_assert(
            std::is_same<util::remove_cvref_t<Mat_>, matrix_type>::value,
            "Matrix type mismatch");
        if (!is_computed_) {
            lu_store_ = std::forward<matrix_type>(mat);
            this->cast().compute_impl();
            is_computed_ = true;
        }
    }

    template <typename Vec>
    util::remove_cvref_t<Vec> solve(Vec&& vec) const {
        util::remove_cvref_t<Vec> vec_tmp(std::forward<Vec>(vec));
        solve_in_place(vec_tmp);
        return vec_tmp;  // Thanks to NRVO, no copy/move will be performed here.
        // But `return std::move(solve_in_place(vec_tmp));` will do extra moves.
    }

    template <typename Iter>  // TODO: Accept raw point
    void solve_in_place(Iter&& iter) const {
        this->cast().solve_in_place_impl(iter);
    }

   protected:
    bool is_computed_;
    matrix_type lu_store_;

    // Get the type of parameter in array subscript operator of given type U. It
    // is assumed that type U is either a container or an iterator.
    template <typename U>
    using get_size_type = typename U::size_type;
    template <typename U>
    using get_difference_type = typename U::difference_type;
    template <typename U>
    using ind_type_ =
        typename util::lazy_conditional<util::is_iterable<U>::value,
                                        get_size_type,
                                        get_difference_type,
                                        U>::type;
};

template <typename>
class BandLU;

template <typename... Ts>
class BandLU<BandMatrix<Ts...>> : public BandLUBase<BandLU, BandMatrix<Ts...>> {
   public:
    using base_type = BandLUBase<intp::BandLU, BandMatrix<Ts...>>;
    using matrix_type = typename base_type::matrix_type;
    using size_type = typename matrix_type::size_type;

   private:
    friend base_type;
    using base_type::base_type;
    using base_type::lu_store_;

    void compute_impl() {
        const size_type n = lu_store_.dim();
        const size_type p = lu_store_.lower_band_width();
        const size_type q = lu_store_.upper_band_width();

        for (size_type k = 0; k < n - 1; ++k) {
            for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                lu_store_(i, k) /= lu_store_(k, k);
            }
            for (size_type j = k + 1; j < std::min(k + q + 1, n); ++j) {
                for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                    lu_store_(i, j) -= lu_store_(i, k) * lu_store_(k, j);
                }
            }
        }
    }

    template <typename Iter>
    void solve_in_place_impl(Iter& iter) const {
        const size_type n = lu_store_.dim();
        const size_type p = lu_store_.lower_band_width();
        const size_type q = lu_store_.upper_band_width();

        using ind_type =
            typename base_type::template ind_type_<util::remove_cvref_t<Iter>>;
        // applying l matrix
        for (size_type j = 0; j < n; ++j) {
            for (size_type i = j + 1; i < std::min(j + p + 1, n); ++i) {
                iter[static_cast<ind_type>(i)] -=
                    lu_store_(i, j) * iter[static_cast<ind_type>(j)];
            }
        }
        // applying u matrix
        for (size_type j = n - 1; j < n; --j) {
            iter[static_cast<ind_type>(j)] /= lu_store_(j, j);
            for (size_type i = j < q ? 0 : j - q; i < j; ++i) {
                iter[static_cast<ind_type>(i)] -=
                    lu_store_(i, j) * iter[static_cast<ind_type>(j)];
            }
        }
    }
};

template <typename... Ts>
class BandLU<ExtendedBandMatrix<Ts...>>
    : public BandLUBase<BandLU, ExtendedBandMatrix<Ts...>> {
   public:
    using base_type = BandLUBase<intp::BandLU, ExtendedBandMatrix<Ts...>>;
    using matrix_type = typename base_type::matrix_type;
    using size_type = typename matrix_type::size_type;

   private:
    friend base_type;
    using base_type::base_type;
    using base_type::lu_store_;

    void compute_impl() {
        const size_type n = lu_store_.dim();
        const size_type p = lu_store_.lower_band_width();
        const size_type q = lu_store_.upper_band_width();

        for (size_type k = 0; k < n - 1; ++k) {
            // update main bands
            for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                lu_store_.main_bands_val(i, k) /=
                    lu_store_.main_bands_val(k, k);
            }

            // update bottom side bands
            for (size_type i = std::max(n - q, k + p + 1); i < n; ++i) {
                lu_store_.side_bands_val(i, k) /=
                    lu_store_.main_bands_val(k, k);
            }

            // update main bands
            for (size_type j = k + 1; j < std::min(k + q + 1, n); ++j) {
                for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                    lu_store_.main_bands_val(i, j) -=
                        lu_store_.main_bands_val(i, k) *
                        lu_store_.main_bands_val(k, j);
                }
            }

            // update upper right corner due to right side bands
            for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                for (size_type j = std::max(n - p, k + q + 1); j < n; ++j) {
                    lu_store_(i, j) -= lu_store_.main_bands_val(i, k) *
                                       lu_store_.side_bands_val(k, j);
                }
            }

            // update lower left corner due to bottom side bands
            for (size_type j = k + 1; j < std::min(k + q + 1, n); ++j) {
                for (size_type i = std::max(n - q, k + p + 1); i < n; ++i) {
                    lu_store_(i, j) -= lu_store_.side_bands_val(i, k) *
                                       lu_store_.main_bands_val(k, j);
                }
            }

            // update main bands due to side bands
            if (k < std::max(n - p - 1, n - q - 1)) {
                for (size_type i = std::max(n - q, k + p + 1); i < n; ++i) {
                    for (size_type j = std::max(n - p, k + q + 1); j < n; ++j) {
                        lu_store_.main_bands_val(i, j) -=
                            lu_store_.side_bands_val(i, k) *
                            lu_store_.side_bands_val(k, j);
                    }
                }
            }
        }
    }

    template <typename Iter>
    void solve_in_place_impl(Iter& iter) const {
        const size_type n = lu_store_.dim();
        const size_type p = lu_store_.lower_band_width();
        const size_type q = lu_store_.upper_band_width();

        using ind_type =
            typename base_type::template ind_type_<util::remove_cvref_t<Iter>>;

        // apply l matrix
        for (size_type j = 0; j < n; ++j) {
            for (size_type i = j + 1; i < std::min(j + p + 1, n); ++i) {
                iter[static_cast<ind_type>(i)] -=
                    lu_store_.main_bands_val(i, j) *
                    iter[static_cast<ind_type>(j)];
            }

            // bottom side bands
            if (j < n - p - 1) {
                for (size_type i = std::max(n - q, j + p + 1); i < n; ++i) {
                    iter[static_cast<ind_type>(i)] -=
                        lu_store_.side_bands_val(i, j) *
                        iter[static_cast<ind_type>(j)];
                }
            }
        }
        // apply u matrix
        for (size_type j = n - 1; j < n; --j) {
            iter[static_cast<ind_type>(j)] /= lu_store_.main_bands_val(j, j);
            for (size_type i = j < q ? 0 : j - q; i < j; ++i) {
                iter[static_cast<ind_type>(i)] -=
                    lu_store_.main_bands_val(i, j) *
                    iter[static_cast<ind_type>(j)];
            }

            // right side bands
            if (j > n - p - 1) {
                for (size_type i = 0; i < j - q; ++i) {
                    iter[static_cast<ind_type>(i)] -=
                        lu_store_.side_bands_val(i, j) *
                        iter[static_cast<ind_type>(j)];
                }
            }
        }
    }
};

}  // namespace intp
