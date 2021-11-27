#pragma once

#include <iostream>
#include <type_traits>
#include <vector>

/**
 * @brief Row-major matrix
 *
 * @tparam T
 */
template <typename T>
class Matrix {
   public:
    using size_type = unsigned int;
    using val_type = T;

   private:
    std::vector<val_type> _data;
    size_type row, col;

   public:
    Matrix(size_type row, size_type col)
        : row(row), col(col), _data(row * col, val_type{0}) {}
    Matrix(size_type dim) : Matrix(dim, dim) {}

    // element access

    val_type& operator()(size_type r, size_type c) {
        return _data[r * col + c];
    }

    // properties

    std::pair<size_type, size_type> dim() const {
        return std::make_pair(row, col);
    }
    bool is_square_matrix() const { return row == col; }
};

/**
 * @brief Band square matrix A is stored in (1+p+q) by n matrix B, with A_{i,j}=
 * B_{i,j-i+p}, i.e. each diagonal is stored as a column in B, aligned by row
 * index.
 *
 * @tparam T value type of matrix element
 */
template <typename T>
class BandMatrix {
   public:
    using size_type = unsigned int;
    using val_type = T;
    using mat_type = BandMatrix<val_type>;

   private:
    Matrix<val_type> _bands;

   public:
    const size_type n;
    const size_type p, q;

    // Create a zero band matrix with given dimension, lower and upper
    // bandwidth.
    BandMatrix(size_type dim, size_type lower, size_type upper)
        : n(dim), p(lower), q(upper), _bands(dim, 1 + lower + upper) {}

    /**
     * @brief Return read/write reference to matrix element, indexes are
     * zero-based.
     *
     * @param i row index
     * @param j column index
     * @return val_type&
     */
    val_type& operator()(size_type i, size_type j) {
        return _bands(i, j - i + p);
    }

    /**
     * @brief Return a new matrix with its (i>j) elements containing L matrix
     * and the rest containing U matrix.
     *
     * @return BandMatrix&
     */
    mat_type LU_decompose() & {
        mat_type tmp(*this);

        return std::move(tmp).LU_decompose();
    }
    /**
     * @brief In place LU decomposition, overwrites the matrix itself with its
     * (i>j) elements containing L matrix and the rest containing U matrix.
     *
     * @return BandMatrix
     */
    mat_type LU_decompose() && {
        for (size_type k = 0; k < n - 1; ++k) {
            for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                operator()(i, k) /= operator()(k, k);
            }
            for (size_type j = k + 1; j < std::min(k + q + 1, n); ++j) {
                for (size_type i = k + 1; i < std::min(k + p + 1, n); ++i) {
                    operator()(i, j) -= operator()(i, k) * operator()(k, j);
                }
            }
        }
        return *this;
    }

    /**
     * @brief Solve linear system A.x = b, the matrix will be LU decomposed
     * in-place.
     *
     * @tparam Vec Used for match lvalue/rvalue reference of b
     * @param b will be rewritten to x if it is a rvalue
     * @return Vec without reference
     */
    template <typename Vec>
    typename std::decay<Vec>::type linear_solve(Vec&& b) && {
        auto lu = std::move(*this).LU_decompose();
        typename std::decay<Vec>::type bb{std::forward<Vec>(b)};

        for (int j = 0; j < n; ++j) {
            for (int i = j + 1; i < std::min(j + p + 1, n); ++i) {
                bb[i] -= lu(i, j) * bb[j];
            }
        }
        for (int j = n - 1; j >= 0; --j) {
            bb[j] /= lu(j, j);
            for (int i = j < q ? 0 : j - q; i < j; ++i) {
                bb[i] -= lu(i, j) * bb[j];
            }
        }

        return bb;
    }

    /**
     * @brief Solve linear system A.x = b
     *
     * @tparam Vec Used for match lvalue/rvalue reference of b
     * @param b will be rewritten to x if it is a rvalue
     * @return Vec without reference
     */
    template <typename Vec>
    typename std::decay<Vec>::type linear_solve(Vec&& b) & {
        mat_type tmp(*this);
        return std::move(tmp).linear_solve(std::forward<Vec>(b));
    }

    template <typename Vec>
    typename std::decay<Vec>::type operator*(const Vec& x) {
        typename std::decay<Vec>::type xx(x.size());
        for (int i = 0; i < x.size(); ++i) {
            for (int j = p > i ? p - i : 0, k = i > p ? i - p : 0;
                 j < std::min(p + q + 1, n + p - i); ++j, ++k) {
                xx[i] += _bands(i, j) * x[k];
            }
        }
        return xx;
    }
};
