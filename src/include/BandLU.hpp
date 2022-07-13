#pragma once

#include <type_traits>

#include "BandMatrix.hpp"
#include "util.hpp"

namespace intp {

/**
 * @brief Base class of band matrix LU solver. The actual solver should inherits
 * from it and implement corresponding impl methods, per the CRTP.
 *
 * @tparam Solver
 * @tparam Matrix
 */
template <template <typename> class Solver, typename Matrix>
class BandLUBase : public util::CRTP<Solver<Matrix>> {
   public:
    using matrix_type = Matrix;

    BandLUBase() noexcept : __is_computed(false) {}

    template <typename _Mat>
    BandLUBase(_Mat&& mat) : __lu_store(std::forward<_Mat>(mat)) {
        this->cast().compute_impl();
        __is_computed = true;
    }

    template <typename _Mat>
    void compute(_Mat&& mat) {
        static_assert(
            std::is_same<util::remove_cvref_t<_Mat>, matrix_type>::value);
        if (!__is_computed) {
            __lu_store = std::forward<matrix_type>(mat);
            this->cast().compute_impl();
            __is_computed = true;
        }
    }

    template <typename Vec>
    util::remove_cvref_t<Vec> solve(Vec&& vec) const {
        util::remove_cvref_t<Vec> vec_tmp(std::forward<Vec>(vec));
        solve_in_place(vec_tmp);
        return vec_tmp;  // Thanks to NRVO, no copy/move will be performed here.
        // But `return std::move(solve_in_place(vec_tmp));` will do extra moves.
    }

    template <typename Vec>
    void solve_in_place(Vec& vec) const {
        this->cast().solve_in_place_impl(vec);
    }

   protected:
    bool __is_computed;
    matrix_type __lu_store;
};

template <typename>
class BandLU;

template <typename T>
class BandLU<BandMatrix<T>> : public BandLUBase<BandLU, BandMatrix<T>> {
   public:
    using base_type = BandLUBase<intp::BandLU, BandMatrix<T>>;
    using matrix_type = BandMatrix<T>;
    using size_type = typename matrix_type::size_type;

   private:
    friend base_type;
    using base_type::__lu_store;
    using base_type::base_type;

    void compute_impl() {
        size_type n = __lu_store.dim();
        size_type p = __lu_store.lower_band_width();
        size_type q = __lu_store.upper_band_width();

        for (size_type k = 0; k < n - 1; ++k) {
            for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                __lu_store(i, k) /= __lu_store(k, k);
            }
            for (size_type j = k + 1; j < std::min(k + q + 1, n); ++j) {
                for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                    __lu_store(i, j) -= __lu_store(i, k) * __lu_store(k, j);
                }
            }
        }
    }

    template <typename Vec>
    void solve_in_place_impl(Vec& vec) const {
        size_type n = __lu_store.dim();
        size_type p = __lu_store.lower_band_width();
        size_type q = __lu_store.upper_band_width();

        // applying l matrix
        for (size_type j = 0; j < n; ++j) {
            for (size_type i = j + 1; i < std::min(j + p + 1, n); ++i) {
                vec[i] -= __lu_store(i, j) * vec[j];
            }
        }
        // applying u matrix
        for (size_type j = n - 1; j < n; --j) {
            vec[j] /= __lu_store(j, j);
            for (size_type i = j < q ? 0 : j - q; i < j; ++i) {
                vec[i] -= __lu_store(i, j) * vec[j];
            }
        }
    }
};

template <typename T>
class BandLU<ExtendedBandMatrix<T>>
    : public BandLUBase<BandLU, ExtendedBandMatrix<T>> {
   public:
    using base_type = BandLUBase<intp::BandLU, ExtendedBandMatrix<T>>;
    using matrix_type = ExtendedBandMatrix<T>;
    using size_type = typename matrix_type::size_type;

   private:
    friend base_type;
    using base_type::__lu_store;
    using base_type::base_type;

    void compute_impl() {
        size_type n = __lu_store.dim();
        size_type p = __lu_store.lower_band_width();
        size_type q = __lu_store.upper_band_width();

        for (size_type k = 0; k < n - 1; ++k) {
            // update main bands
            for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                __lu_store.main_bands_val(i, k) /=
                    __lu_store.main_bands_val(k, k);
            }

            // update bottom side bands
            for (size_type i = std::max(n - q, k + p + 1); i < n; ++i) {
                __lu_store.side_bands_val(i, k) /=
                    __lu_store.main_bands_val(k, k);
            }

            // update main bands
            for (size_type j = k + 1; j < std::min(k + q + 1, n); ++j) {
                for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                    __lu_store.main_bands_val(i, j) -=
                        __lu_store.main_bands_val(i, k) *
                        __lu_store.main_bands_val(k, j);
                }
            }

            // update upper right corner due to right side bands
            for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                for (size_type j = std::max(n - p, k + q + 1); j < n; ++j) {
                    __lu_store(i, j) -= __lu_store.main_bands_val(i, k) *
                                        __lu_store.side_bands_val(k, j);
                }
            }

            // update lower left corner due to bottom side bands
            for (size_type j = k + 1; j < std::min(k + q + 1, n); ++j) {
                for (size_type i = std::max(n - q, k + p + 1); i < n; ++i) {
                    __lu_store(i, j) -= __lu_store.side_bands_val(i, k) *
                                        __lu_store.main_bands_val(k, j);
                }
            }

            // update main bands due to side bands
            if (k < std::max(n - p - 1, n - q - 1)) {
                for (size_type i = std::max(n - q, k + p + 1); i < n; ++i) {
                    for (size_type j = std::max(n - p, k + q + 1); j < n; ++j) {
                        __lu_store.main_bands_val(i, j) -=
                            __lu_store.side_bands_val(i, k) *
                            __lu_store.side_bands_val(k, j);
                    }
                }
            }
        }
    }

    template <typename Vec>
    void solve_in_place_impl(Vec& vec) const {
        size_type n = __lu_store.dim();
        size_type p = __lu_store.lower_band_width();
        size_type q = __lu_store.upper_band_width();

        // apply l matrix
        for (size_type j = 0; j < n; ++j) {
            for (size_type i = j + 1; i < std::min(j + p + 1, n); ++i) {
                vec[i] -= __lu_store.main_bands_val(i, j) * vec[j];
            }

            // bottom side bands
            if (j < n - p - 1) {
                for (size_type i = std::max(n - q, j + p + 1); i < n; ++i) {
                    vec[i] -= __lu_store.side_bands_val(i, j) * vec[j];
                }
            }
        }
        // apply u matrix
        for (size_type j = n - 1; j < n; --j) {
            vec[j] /= __lu_store.main_bands_val(j, j);
            for (size_type i = j < q ? 0 : j - q; i < j; ++i) {
                vec[i] -= __lu_store.main_bands_val(i, j) * vec[j];
            }

            // right side bands
            if (j > n - p - 1) {
                for (size_type i = 0; i < j - q; ++i) {
                    vec[i] -= __lu_store.side_bands_val(i, j) * vec[j];
                }
            }
        }
    }
};

}  // namespace intp
