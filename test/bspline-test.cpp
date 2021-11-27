#include "../src/include/BSpline.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>

using namespace std;

int main() {
    auto knots = {1., 1., 1., 1., 3., 5., 5., 5., 5.};
    auto cp = {1., 5.333333, -5.33333, 1.33333, 3.};

    BSpline<double> spline3{3, knots.begin(), knots.end(), cp.begin(),
                            cp.end()};

    // Values pre-computed by MMA
    auto val = {1.,       2.024, 2.552,  2.668, 2.456,  2.,     1.384,
                0.692001, 0.008, -0.584, -1.,   -1.176, -1.128, -0.892,
                -0.504,   0.,    0.584,  1.212, 1.848,  2.456};

    cout << "\nBSpline Test:\n";
    double d;
    auto iter = val.begin();
    for (double x = 1.; x < 5.; x += 0.2, ++iter) {
        d += (spline3(x) - *iter) * (spline3(x) - *iter);
        cout << x << " " << spline3(x) << '\n';
    }

    return std::sqrt(d / val.size()) < 1e-6 ? 0 : 1;
}