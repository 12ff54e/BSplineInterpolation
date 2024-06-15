#pragma once

#include <array>
#include <cmath>
#include <type_traits>
#include <utility>

/**
 * @brief A vector in Euclidean space.
 *
 * @tparam D Dimension of the underlying space.
 * @tparam T Type of coordinate.
 */
template <std::size_t D, typename T>
class Vec;

template <std::size_t D, typename T>
class VecBase {
    static_assert(D > 0, "D cannot be 0.");

   public:
    static constexpr std::size_t dim = D;
    using value_type = T;
    using vec_type = Vec<dim, value_type>;

   protected:
    std::array<value_type, dim> coord;

   private:
    template <std::size_t... indices>
    static bool equal_aux_(intp::util::index_sequence<indices...>,
                           const VecBase& lhs,
                           const VecBase& rhs) {
#if __cplusplus >= 201703L
        return (... && (lhs[indices] == rhs[indices]));
#else
        for (std::size_t i = 0; i < dim; ++i) {
            if (!(lhs[i] == rhs[i])) { return false; }
        }
        return true;
#endif
    }

   public:
    constexpr VecBase() = default;

    template <typename... Ts,
              typename = typename std::enable_if<sizeof...(Ts) <= dim>::type>
    constexpr VecBase(Ts... vs) noexcept
        : coord{static_cast<value_type>(vs)...} {}

    /**
     * @brief Conversion constructor from supported type
     *
     * @param other another vec
     */
    template <typename U,
              typename = typename std::enable_if<
                  std::is_convertible<U, T>::value>::type>
    VecBase(VecBase<dim, U> other) noexcept {
        for (std::size_t i = 0; i < dim; ++i) {
            coord[i] = static_cast<value_type>(other[i]);
        }
    }

    // zero element

    static constexpr vec_type zero() noexcept {
        return VecBase<dim, value_type>{};
    }

    // element access

    T& operator[](std::size_t i) noexcept {
        return coord[i];
    }

    const T& operator[](std::size_t i) const noexcept {
        return coord[i];
    }

    // comparison operations

    friend bool operator==(const vec_type& lhs, const vec_type& rhs) {
        return equal_aux_(intp::util::make_index_sequence<dim>{}, lhs, rhs);
    }

    friend bool operator!=(const vec_type& lhs, const vec_type& rhs) {
        return !operator==(lhs, rhs);
    }

    // arithmetic operations

    friend vec_type& operator+=(vec_type& lhs, const vec_type& rhs) noexcept {
        for (std::size_t i = 0; i < dim; ++i) { lhs[i] += rhs[i]; }
        return lhs;
    }
    friend vec_type& operator-=(vec_type& lhs, const vec_type& rhs) noexcept {
        for (std::size_t i = 0; i < dim; ++i) { lhs[i] -= rhs[i]; }
        return lhs;
    }

    friend vec_type& operator*=(vec_type& lhs, const T& scalar) noexcept {
        for (std::size_t i = 0; i < dim; ++i) { lhs[i] *= scalar; }
        return lhs;
    }

    friend vec_type& operator/=(vec_type& lhs, const T& scalar) noexcept {
        for (std::size_t i = 0; i < dim; ++i) { lhs[i] /= scalar; }
        return lhs;
    }

    friend vec_type operator+(vec_type lhs, const vec_type& rhs) noexcept {
        return lhs += rhs;
    }
    friend vec_type operator-(vec_type lhs, const vec_type& rhs) noexcept {
        return lhs -= rhs;
    }

    friend vec_type operator*(vec_type vec, const T& scalar) noexcept {
        return vec *= scalar;
    }
    friend vec_type operator*(const T& scalar, vec_type vec) noexcept {
        return vec *= scalar;
    }

    friend vec_type operator/(vec_type vec, const T& scalar) noexcept {
        return vec /= scalar;
    }

    // conversion operator

    // convert to the underlying coordinate array
    constexpr operator std::array<value_type, dim>() const {
        return coord;
    }

    // convert to derived class
    constexpr operator vec_type() const& {
        return static_cast<vec_type>(*this);
    }
    operator vec_type() && {
        return static_cast<vec_type&&>(std::move(*this));
    }

    // properties

    T L2_norm_square() const {
        T norm{};
        for (auto& c : coord) { norm += c * c; }
        return norm;
    }

    T mag() const {
        return std::sqrt(L2_norm_square());
    }
};

/**
 * @brief Generic Vector type of any dimension, default to consisting double as
 * coordinates
 *
 * @tparam D Dimension
 * @tparam T Underlying types of each dimension
 */
template <std::size_t D, typename T = double>
class Vec : public VecBase<D, T> {};

template <typename T>
class Vec<2, T> : public VecBase<2, T> {
   public:
    using VecBase<2, T>::VecBase;

    T& x() noexcept { return this->coord[0]; }
    T& y() noexcept { return this->coord[1]; }

    const T& x() const noexcept { return this->coord[0]; }
    const T& y() const noexcept { return this->coord[1]; }
};
