#ifndef INTP_BENCH_H_
#define INTP_BENCH_H_

#include <cstddef>

template <typename T>
void do_not_optimize(const T& x) {
#if (defined(__GNUC__) || defined(__clang__)) && !defined(_MSC_VER)
    asm volatile("" : : "r,m"(x));
#else
    // Volatile reads keep every byte of the result observable.
    const volatile unsigned char* bytes =
        reinterpret_cast<const volatile unsigned char*>(&x);
    for (std::size_t i = 0; i < sizeof(T); ++i) { (void)bytes[i]; }
#endif
}

#endif  // INTP_BENCH_H_
