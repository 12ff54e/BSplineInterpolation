#pragma once

#include <iostream>
#include <type_traits>
#include <vector>

#include "Mesh.hpp"

namespace intp {

/**
 * @brief Band square matrix A is stored in n by (1+p+q) matrix B, with A_{i,j}=
 * B_{j,i+q-j}, i.e. each diagonal is stored as a column in B, aligned by column
 * index in A.
 *
 * @tparam T value type of matrix element
 */
template <typename T, typename Alloc = std::allocator<T>>
class BandMatrix {
   public:
    using size_type = size_t;
    using val_type = T;
    using allocator_type = Alloc;
    using matrix_type = BandMatrix<val_type, allocator_type>;

    // Create a zero band matrix with given dimension, lower and upper
    // bandwidth.
    BandMatrix(size_type n, size_type p, size_type q)
        : n_(n), p_(p), q_(q), bands_{n_, 1 + p_ + q_} {}

    BandMatrix() : BandMatrix(0, 0, 0) {}

    // properties

    size_type dim() const noexcept { return n_; }

    size_type lower_band_width() const noexcept { return p_; }

    size_type upper_band_width() const noexcept { return q_; }

    /**
     * @brief Return read/write reference to matrix element, indices are
     * zero-based.
     *
     * @param i row index
     * @param j column index
     * @return val_type&
     */
    val_type& operator()(size_type i, size_type j) {
        CUSTOM_ASSERT(j + p_ >= i && i + q_ >= j,
                      "Given i and j not in main bands.");
        return bands_(j, i + q_ - j);
    }

    val_type operator()(size_type i, size_type j) const {
        CUSTOM_ASSERT(j + p_ >= i && i + q_ >= j,
                      "Given i and j not in main bands.");
        return bands_(j, i + q_ - j);
    }

    /**
     * @brief Matrix-Vector multiplication
     *
     * @param x vector to be multiplied
     */
    template <typename Iter>
    util::remove_cvref_t<Iter> operator*(const Iter& x) const {
        util::remove_cvref_t<Iter> xx(x.size());
        for (size_type i = 0; i < x.size(); ++i) {
            for (size_type j = p_ > i ? p_ - i : 0, k = i > p_ ? i - p_ : 0;
                 j < std::min(p_ + q_ + 1, n_ + p_ - i); ++j, ++k) {
                xx[i] += bands_(i, j) * x[k];
            }
        }
        return xx;
    }

    /**
     * @brief Insertion operator, used for debug mostly.
     *
     */
    friend std::ostream& operator<<(std::ostream& os, const BandMatrix& mat) {
        for (size_t i = 0; i < mat.n_; ++i) {
            for (size_t j = i > mat.p_ ? i - mat.p_ : 0;
                 j<i + mat.q_ + 1> mat.n_ ? mat.n_ : i + mat.q_ + 1; ++j) {
                os << "{" << i << ", " << j << "}->" << mat(i, j) << '\n';
            }
        }
        return os;
    }

   protected:
    size_type n_;
    size_type p_, q_;
    Mesh<val_type, 2, allocator_type> bands_;
};

template <typename T, typename Alloc = std::allocator<T>>
class ExtendedBandMatrix : public BandMatrix<T, Alloc> {
   public:
    using base_type = BandMatrix<T, Alloc>;
    using size_type = typename base_type::size_type;
    using val_type = typename base_type::val_type;
    using allocator_type = typename base_type::allocator_type;

    ExtendedBandMatrix(size_type dim, size_type lower, size_type upper)
        : base_type(dim, lower, upper),
          right_side_bands_{dim - upper - 1, lower},
          bottom_side_bands_{dim - lower - 1, upper} {}

    ExtendedBandMatrix() : ExtendedBandMatrix(1, 0, 0) {}

    val_type& main_bands_val(size_type i, size_type j) {
        return base_type::operator()(i, j);
    }

    val_type main_bands_val(size_type i, size_type j) const {
        return base_type::operator()(i, j);
    }

    val_type& side_bands_val(size_type i, size_type j) {
        CUSTOM_ASSERT(j >= std::max(n_ - p_, i + q_ + 1) ||
                          i >= std::max(n_ - q_, j + p_ + 1),
                      "Given i and j not in side bands.");
        return j > i + q_ ? right_side_bands_(i, j + p_ - n_)
                          : bottom_side_bands_(j, i + q_ - n_);
    }

    val_type side_bands_val(size_type i, size_type j) const {
        CUSTOM_ASSERT(j >= std::max(n_ - p_, i + q_ + 1) ||
                          i >= std::max(n_ - q_, j + p_ + 1),
                      "Given i and j not in side bands.");
        return j > i + q_ ? right_side_bands_(i, j + p_ - n_)
                          : bottom_side_bands_(j, i + q_ - n_);
    }

    val_type& operator()(size_type i, size_type j) {
        return j > i + q_ || i > j + p_ ? side_bands_val(i, j)
                                        : main_bands_val(i, j);
    }

    val_type operator()(size_type i, size_type j) const {
        return j > i + q_ || i > j + p_ ? side_bands_val(i, j)
                                        : main_bands_val(i, j);
    }

    template <typename Iter>
    util::remove_cvref_t<Iter> operator*(const Iter& x) const {
        util::remove_cvref_t<Iter> xx(x.size());
        for (size_type i = 0; i < x.size(); ++i) {
            for (size_type j = i > p_ ? i - p_ : 0;
                 j < std::min(i + q_ + 1, n_); ++j) {
                xx[i] += main_bands_val(i, j) * x[j];
            }

            // right side bands
            if (i < n_ - q_ - 1) {
                for (size_type j = std::max(n_ - p_, i + q_ + 1); j < n_; ++j) {
                    xx[i] += side_bands_val(i, j) * x[j];
                }
            }
        }

        // bottom side bands
        for (size_type j = 0; j < n_ - p_ - 1; ++j) {
            for (size_type i = std::max(n_ - q_, j + p_ + 1); i < n_; ++i) {
                xx[i] += side_bands_val(i, j) * x[j];
            }
        }
        return xx;
    }

   private:
    Mesh<val_type, 2, allocator_type> right_side_bands_;
    Mesh<val_type, 2, allocator_type> bottom_side_bands_;

    using base_type::n_;
    using base_type::p_;
    using base_type::q_;
};

}  // namespace intp
