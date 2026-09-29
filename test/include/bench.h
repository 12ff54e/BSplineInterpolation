#ifndef INTP_BENCH_H_
#define INTP_BENCH_H_

template <typename T>
void do_not_optimize(const T& x) {
    asm volatile("" : : "r,m"(x));
}

#endif  // INTP_BENCH_H_
