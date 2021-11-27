#include "../src/include/Interpolation.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>

using namespace std;

int main() {
    vector<double> f{{0.905515, 0.894638, -0.433134, 0.43131, -0.131052,
                      0.262974, 0.423888, -0.562671, -0.915567, -0.261017,
                      -0.47915, -0.00939326, -0.445962}};
    InterpolationFunction<double> interp(
        3, std::make_pair(1., (double)(f.size())), f.begin(), f.end());
    cout << "\nInterpolation Test:\n";
    for (double x = 1.; x < (double)(f.size()); x += 0.2) {
        cout << x << " " << interp(x) << '\n';
    }
    return 0;
}