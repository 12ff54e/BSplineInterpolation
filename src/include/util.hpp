#pragma once

#include <array>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace intp {

namespace util {

/**
 * @brief Polyfill for C++14 intergers_sequence, but with [T = size_t] only
 *
 * @tparam Indices
 */
template <size_t... Indices>
struct index_sequence {
    using val_type = size_t;
    const static size_t size = sizeof...(Indices);
};

template <size_t N, size_t... Indices>
struct make_index_sequence_impl
    : make_index_sequence_impl<N - 1, N - 1, Indices...> {};

template <size_t... Indices>
struct make_index_sequence_impl<0, Indices...> {
    using type = index_sequence<Indices...>;
};

template <size_t N>
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

template <typename Func, typename... Args, unsigned... indices>
void dispatch_indexed_helper(index_sequence<indices...>,
                             Func& func,
                             Args&&... args) {
    // polyfill of C++17 fold expression over comma
    std::array<std::nullptr_t, sizeof...(Args)>{
        (func(indices, std::forward<Args>(args)), nullptr)...};
}

/**
 * @brief dispatch_indexed(f, x0, x1, ...) invokes f(0, x0), f(1, x1), ..., and
 * ignores their return values.
 *
 * @param func
 * @param args
 */
template <typename Func, typename... Args>
void dispatch_indexed(Func&& func, Args&&... args) {
    dispatch_indexed_helper(util::make_index_sequence_for<Args...>{}, func,
                            std::forward<Args>(args)...);
}

#ifdef STACK_ALLOCATOR

/**
 * @brief A simple stack allocator, with fixed size and a LIFO allocation
 * strategy, by courtesy of Charles Salvia, in his SO answer
 * https://stackoverflow.com/a/28574062/7255197 .
 *
 * @tparam T allocation object type
 * @tparam N buffer size
 */
template <typename T, std::size_t N>
class stack_allocator {
   public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;

    using const_void_pointer = const void*;

   private:
    pointer m_begin;
    pointer m_end;
    pointer m_stack_pointer;

   public:
    explicit stack_allocator(pointer buffer)
        : m_begin(buffer), m_end(buffer + N), m_stack_pointer(buffer){};

    template <typename U>
    stack_allocator(const stack_allocator<U, N>& other)
        : m_begin(other.m_begin),
          m_end(other.m_end),
          m_stack_pointer(other.m_stack_pointer) {}

    constexpr static size_type capacity() { return N; }

    pointer allocate(size_type n,
                     const_void_pointer hint = const_void_pointer()) {
        if (n <= size_type(m_end - m_stack_pointer)) {
            pointer result = m_stack_pointer;
            m_stack_pointer += n;
            return result;
        }
        throw std::bad_alloc{};
    }

    void deallocate(pointer ptr, size_type n) { m_stack_pointer -= n; }

    size_type max_size() const noexcept { return N; }

    pointer address(reference x) const noexcept { return std::addressof(x); }

    const_pointer address(const_reference x) const noexcept {
        return std::addressof(x);
    }

    template <typename U>
    struct rebind {
        using other = stack_allocator<U, N>;
    };

    pointer buffer() const noexcept { return m_begin; }
};

template <typename T, std::size_t N, typename U>
bool operator==(const stack_allocator<T, N>& lhs,
                const stack_allocator<U, N>& rhs) noexcept {
    return lhs.buffer() == rhs.buffer();
}

template <typename T, std::size_t N, typename U>
bool operator!=(const stack_allocator<T, N>& lhs,
                const stack_allocator<U, N>& rhs) noexcept {
    return !(lhs == rhs);
}

#endif

struct _is_iteratable_impl {
    template <typename _T,
              typename = typename std::enable_if<std::is_convertible<
                  typename std::iterator_traits<
                      decltype(std::declval<_T&>().begin())>::iterator_category,
                  std::input_iterator_tag>::value>::type,
              typename = typename std::enable_if<std::is_convertible<
                  typename std::iterator_traits<
                      decltype(std::declval<_T&>().end())>::iterator_category,
                  std::input_iterator_tag>::value>::type>
    static std::true_type __test(int);

    template <typename>
    static std::false_type __test(...);
};

struct _is_indexed_impl {
    template <typename _T,
              typename = decltype(std::declval<const _T&>().operator[](0))>
    static std::true_type __test(int);

    template <typename>
    static std::false_type __test(...);
};

/**
 * @brief Check if a type has begin() and end() method that returns iterator
 *
 * @tparam T a type to check
 */
template <typename T>
struct is_iteratable : _is_iteratable_impl {
    static constexpr bool value = decltype(__test<T>(0))::value;
};

template <typename T>
struct is_indexed : _is_indexed_impl {
    static constexpr bool value = decltype(__test<T>(0))::value;
};

/**
 * @brief Polyfill for C++20 stl function with the same name
 *
 * @tparam Vec
 */
template <typename Vec>
using remove_cvref_t =
    typename std::remove_cv<typename std::remove_reference<Vec>::type>::type;

/**s
 * @brief CRTP helper, used for downward casting.
 *
 */
template <typename T, typename...>
struct CRTP {
    T& cast() { return static_cast<T&>(*this); }
    const T& cast() const { return static_cast<const T&>(*this); }
};

void custom_assert(bool assertion, const char* msg) {
#ifdef _DEBUG
    if (!assertion) { throw std::runtime_error(msg); }
#endif
}

}  // namespace util

}  // namespace intp
