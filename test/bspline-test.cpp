#include "../src/include/BSpline.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>

using namespace std;

int main() {
    auto knots = {1., 1., 1., 1., 3., 5., 5., 5., 5.};
    auto cp = {1., 16. / 3, -16. / 3, 4. / 3, 3.};

    BSpline<double> spline3{3, knots.begin(), knots.end(), cp.begin(),
                            cp.end()};

    // Values pre-computed by MMA
    auto val = {1.,
                2.0239999999999996,
                2.552,
                2.6679999999999997,
                2.4560000000000004,
                2.,
                1.3839999999999995,
                0.6920000000000003,
                0.00799999999999973,
                -0.5839999999999992,
                -1.,
                -1.176,
                -1.1279999999999997,
                -0.8919999999999999,
                -0.5040000000000004,
                0.,
                0.5840000000000004,
                1.212000000000001,
                1.847999999999999,
                2.4559999999999995};

    cout << "\nBSpline Test:\n";
    double d;
    auto iter = val.begin();
    for (double x = 1.; x < 5.; x += 0.2, ++iter) {
        d += (spline3(x) - *iter) * (spline3(x) - *iter);
        cout << x << " " << spline3(x) << '\n';
    }

    return std::sqrt(d / val.size()) < 1e-14 ? 0 : 1;
}