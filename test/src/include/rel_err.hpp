#ifndef INTP_TEST_RELERR
#define INTP_TEST_RELERR

#include <cmath>
#include <iostream>

#ifdef _DEBUG
#include <iomanip>
#endif

template <typename Func, typename InputIterPt, typename InputIterVal>
double rel_err(const Func& interp,
               std::pair<InputIterPt, InputIterPt> pts,
               std::pair<InputIterVal, InputIterVal> vals) {
    double err{}, l2{};
#ifdef _DEBUG
    auto default_prec = std::cout.precision(17);
    std::cout << "\n[DEBUG] Spline Value           \tExpected\n";
    std::size_t idx = 0;
    constexpr std::size_t print_limit = 5;
#endif
    auto pt_it = pts.first;
    auto val_it = vals.first;
    for (; pt_it != pts.second && val_it != vals.second; ++pt_it, ++val_it) {
        double f = interp(*pt_it);
        err += (f - *val_it) * (f - *val_it);
        l2 += (*val_it) * (*val_it);
#ifdef _DEBUG
        if (idx < print_limit) {
            std::cout << "[DEBUG] " << std::setw(20) << f << ",\t"
                      << std::setw(20) << *val_it << '\n';
        }
#ifndef _TRACE
        ++idx;
#endif
#endif
    }
#ifdef _DEBUG
    std::cout.precision(default_prec);
#ifndef _TRACE
    if (idx > print_limit) {
        std::cout << "[DEBUG] (and " << idx - print_limit
                  << " more entries ...)\n";
    }
#endif
#endif

    return std::sqrt(err / l2);
}

#endif
