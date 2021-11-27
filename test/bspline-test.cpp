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

    cout << "\nBSpline Test:\n";
    for (double x = 1.; x < 5.; x += 0.2) {
        cout << x << " " << spline3(x) << '\n';
    }
    return 0;
}