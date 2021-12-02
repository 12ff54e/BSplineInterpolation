#include "../src/include/BSpline.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>

using namespace std;

int main() {
    int err = 0;

    // 1D B-Spline

    auto knots = {0., 0., 0., 0., .5, 1., 1., 1., 1.};
    auto cp = {1., 16. / 3, -16. / 3, 4. / 3, 3.};

    BSpline<double, 1> spline3(3, Mesh<double, 1>(cp),
                               std::make_pair(knots.begin(), knots.end()));

    // Values pre-computed by MMA
    auto val_1d = {1.,
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

    cout << "\n1D B-Spline Test:\n";
    double d = 0;
    auto iter = val_1d.begin();
    for (double x = 0.; x < 1.; x += 0.05, ++iter) {
        const double f = spline3(x);
        d += (f - *iter) * (f - *iter);
        cout << x << " " << f << '\n';
    }

    err = err == 0 && std::sqrt(d / val_1d.size()) < 1e-14 ? 0 : 1;

    // 2D B-Spline

    // random 5x5 control points
    constexpr array<array<double, 5>, 5> cp2{{{-4, -4, 4, 1, 1},
                                              {0, 2, -2, 0, 5},
                                              {5, -4, -4, -3, -5},
                                              {-2, 4, 2, 2, -3},
                                              {1, 0, 5, -3, 5}}};

    Mesh<double, 2> cp2d(5, 5);
    for (unsigned i = 0; i < 5; ++i) {
        for (unsigned j = 0; j < 5; ++j) { cp2d(i, j) = cp2[i][j]; }
    }

    BSpline<double, 2> spline2d3(3, cp2d, make_pair(knots.begin(), knots.end()),
                                 make_pair(knots.begin(), knots.end()));

    cout << "\n2D B-Spline Test:\n";

    // some random points
    constexpr array<pair<double, double>, 10> coords{
        {{0.049800826265867126, 0.8612651106123255},
         {0.575715550141984, 0.5153782982448287},
         {0.051389608765929795, 0.20636530700031974},
         {0.3406193462245328, 0.6559088769326575},
         {0.6562156301169633, 0.5229827075599092},
         {0.7065591320791091, 0.18017930762978818},
         {0.5929532206588517, 0.015110486908892717},
         {0.35737449937918164, 0.2771152465740585},
         {0.04848381891207554, 0.9244449256192726},
         {0.9740876351532088, 0.9174673912595357}}};

    // and their spline value, pre-computed by MMA
    constexpr array<double, 10> vals_2d{
        {1.3002562089978034, -0.9049532347293099, -1.5715077611939692,
         -1.3555951852638821, -0.19537752208751952, 0.7893147505912226,
         1.3459399154210203, -0.7670721516764836, 1.495128187689829,
         1.5162241194068322}};

    d = 0;
    for (unsigned i = 0; i < vals_2d.size(); ++i) {
        const double f = spline2d3(coords[i].first, coords[i].second);
        d += (f - vals_2d[i]) * (f - vals_2d[i]);
        cout << coords[i].first << " " << coords[i].second << " " << f << '\n';
    }
    err = err == 0 && std::sqrt(d / vals_2d.size()) < 1e-14 ? 0 : 1;

    return err;
}