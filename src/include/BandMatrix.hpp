#pragma once

#include <iostream>
#include <type_traits>
#include <vector>

#include "Mesh.hpp"

namespace intp {

/**
 * @brief Band square matrix A is stored in n by (1+p+q) matrix B, with A_{i,j}=
 * B_{j,i+q-j}, i.e. each diagonal is stored as a row in B, aligned by column
 * index.
 *
 * @tparam T value type of matrix element
 */
template <typename T>
class BandMatrix {
   public:
    using size_type = size_t;
    using val_type = T;
    using mat_type = BandMatrix<val_type>;

    // Create a zero band matrix with given dimension, lower and upper
    // bandwidth.
    BandMatrix(size_type dim, size_type lower, size_type upper)
        : n(dim), p(lower), q(upper), __bands{n, 1 + p + q} {}

    BandMatrix() : BandMatrix(0, 0, 0) {}

    // properties

    size_type dim() const noexcept { return n; }

    size_type lower_band_width() const noexcept { return p; }

    size_type upper_band_width() const noexcept { return q; }

    /**
     * @brief Return read/write reference to matrix element, indices are
     * zero-based.
     *
     * @param i row index
     * @param j column index
     * @return val_type&
     */
    val_type& operator()(size_type i, size_type j) {
        return __bands(j, i + q - j);
    }

    val_type operator()(size_type i, size_type j) const {
        return __bands(j, i + q - j);
    }

    /**
     * @brief Matrix-Vector multiplication
     *
     * @tparam Vec
     * @param x
     */
    template <typename Vec>
    util::remove_cvref_t<Vec> operator*(const Vec& x) const {
        util::remove_cvref_t<Vec> xx(x.size());
        for (size_type i = 0; i < x.size(); ++i) {
            for (size_type j = p > i ? p - i : 0, k = i > p ? i - p : 0;
                 j < std::min(p + q + 1, n + p - i); ++j, ++k) {
                xx[i] += __bands(i, j) * x[k];
            }
        }
        return xx;
    }

   protected:
    size_type n;
    size_type p, q;
    Mesh<val_type, 2> __bands;
};

template <typename T>
class ExtendedBandMatrix : public BandMatrix<T> {
   public:
    using base_type = BandMatrix<T>;
    using size_type = typename base_type::size_type;
    using val_type = typename base_type::val_type;

    ExtendedBandMatrix(size_type dim, size_type lower, size_type upper)
        : base_type(dim, lower, upper),
          __right_side_bands{dim - upper - 1, lower},
          __bottom_side_bands{dim - lower - 1, upper} {}

    ExtendedBandMatrix() : ExtendedBandMatrix(1, 0, 0) {}

    val_type& main_bands_val(size_type i, size_type j) {
        return base_type::operator()(i, j);
    }

    val_type main_bands_val(size_type i, size_type j) const {
        return base_type::operator()(i, j);
    }

    val_type& side_bands_val(size_type i, size_type j) {
        return j > i + q ? __right_side_bands(i, j + p - n)
                         : __bottom_side_bands(j, i + q - n);
    }

    val_type side_bands_val(size_type i, size_type j) const {
        return j > i + q ? __right_side_bands(i, j + p - n)
                         : __bottom_side_bands(j, i + q - n);
    }

    val_type& operator()(size_type i, size_type j) {
        return j > i + q || i > j + p ? side_bands_val(i, j)
                                      : main_bands_val(i, j);
    }

    val_type operator()(size_type i, size_type j) const {
        return j > i + q || i > j + p ? side_bands_val(i, j)
                                      : main_bands_val(i, j);
    }

    template <typename Vec>
    util::remove_cvref_t<Vec> operator*(const Vec& x) const {
        util::remove_cvref_t<Vec> xx(x.size());
        for (size_type i = 0; i < x.size(); ++i) {
            for (size_type j = i > p ? i - p : 0; j < std::min(i + q + 1, n);
                 ++j) {
                xx[i] += main_bands_val(i, j) * x[j];
            }

            // right side bands
            if (i < n - q - 1) {
                for (size_type j = std::max(n - p, i + q + 1); j < n; ++j) {
                    xx[i] += side_bands_val(i, j) * x[j];
                }
            }
        }

        // bottom side bands
        for (size_type j = 0; j < n - p - 1; ++j) {
            for (size_type i = std::max(n - q, j + p + 1); i < n; ++i) {
                xx[i] += side_bands_val(i, j) * x[j];
            }
        }
        return xx;
    }

   private:
    Mesh<val_type, 2> __right_side_bands;
    Mesh<val_type, 2> __bottom_side_bands;

    using base_type::n;
    using base_type::p;
    using base_type::q;
};

}  // namespace intp
