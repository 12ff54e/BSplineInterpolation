#pragma once

namespace util {

/**
 * @brief Polyfill for C++14 interger_sequence, but with [T = unsigned int] only
 *
 * @tparam Indices
 */
template <unsigned... Indices>
struct index_sequence {
    using val_type = unsigned;
    const static unsigned size = sizeof...(Indices);
};

template <unsigned N, unsigned... Indices>
struct make_index_sequence_impl
    : make_index_sequence_impl<N - 1, N - 1, Indices...> {};

template <unsigned... Indices>
struct make_index_sequence_impl<0, Indices...> {
    using type = index_sequence<Indices...>;
};

template <unsigned N>
using make_index_sequence = typename make_index_sequence_impl<N>::type;

template <typename... T>
using make_index_sequence_for = make_index_sequence<sizeof...(T)>;

/**
 * @brief Compile time power for unsigned exponent
 *
 * @param base
 * @param exp
 * @return base^exp in type of base
 */
template <typename T1, typename T2>
constexpr typename std::enable_if<std::is_unsigned<T2>::value, T1>::type pow(
    T1 base,
    T2 exp) {
    return exp == 0 ? T1{1} : base * pow(base, exp - 1);
}

}  // namespace util
