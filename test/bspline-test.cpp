#include "../src/include/BSpline.hpp"
#include "Assertion.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>

using namespace std;

int main() {
    Assertion assertion;

    // 1D B-Spline

    auto knots = {0., 0., 0., 0., .5, 1., 1., 1., 1.};
    auto cp = {1., 16. / 3, -16. / 3, 4. / 3, 3.};

    BSpline<double, 1> spline_1d_3(3, Mesh<double, 1>(cp),
                                   std::make_pair(knots.begin(), knots.end()));

    // Values pre-computed by MMA
    auto vals_1d = {1.,
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

    std::cout << "\n1D B-Spline Test:\n";
    double d = 0;
    auto iter = vals_1d.begin();
    for (double x = 0.; x < 1.; x += 0.05, ++iter) {
        const double f = spline_1d_3(x);
        d += (f - *iter) * (f - *iter);
        std::cout << x << " " << f << '\n';
    }

    assertion(std::sqrt(d / vals_1d.size()) < 1e-14);
    std::cout << "\n1D test "
              << (assertion.status() == 0 ? "succeed" : "failed") << '\n';

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

    BSpline<double, 2> spline_2d_3(3, cp2d,
                                   make_pair(knots.begin(), knots.end()),
                                   make_pair(knots.begin(), knots.end()));

    std::cout << "\n2D B-Spline Test:\n";

    // some random points
    constexpr array<pair<double, double>, 10> coords_2d{
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
        const double f = spline_2d_3(coords_2d[i].first, coords_2d[i].second);
        d += (f - vals_2d[i]) * (f - vals_2d[i]);
        std::cout << coords_2d[i].first << " " << coords_2d[i].second << " "
                  << f << '\n';
    }
    assertion(std::sqrt(d / vals_2d.size()) < 1e-14);
    std::cout << "\n2D test "
              << (assertion.status() == 0 ? "succeed" : "failed") << '\n';

    // 3D Spline

    // random 5x5x5 control points
    constexpr array<array<array<double, 5>, 5>, 5> cp3{{{{{-4, -4, -4, -2, -5},
                                                          {-4, 0, -5, 5, 4},
                                                          {4, -5, -3, 5, -4},
                                                          {3, 4, 0, -3, -1},
                                                          {-5, 2, 1, -5, -2}}},
                                                        {{{5, 4, -3, -5, 5},
                                                          {2, -3, 0, 1, 4},
                                                          {5, 0, -4, 3, -1},
                                                          {2, -2, 3, 0, 2},
                                                          {-1, 1, -4, -3, -2}}},
                                                        {{{-1, 3, -1, 2, 2},
                                                          {-1, -5, -3, -5, -1},
                                                          {3, -3, -5, 2, 0},
                                                          {2, 1, 1, 3, 3},
                                                          {2, 2, 4, 4, -1}}},
                                                        {{{-4, 0, 2, -5, 2},
                                                          {-2, -3, 1, 2, -2},
                                                          {4, -4, 5, -3, 0},
                                                          {-4, -5, 3, -4, 2},
                                                          {2, -4, -2, 1, 5}}},
                                                        {{{-3, 2, 1, 0, -2},
                                                          {2, 3, 2, -5, 5},
                                                          {-1, 2, 3, 0, -5},
                                                          {1, 2, -1, 1, -1},
                                                          {0, -1, 5, 5, 4}}}}};

    Mesh<double, 3> cp3d(5, 5, 5);
    for (unsigned i = 0; i < 5; ++i) {
        for (unsigned j = 0; j < 5; ++j) {
            for (unsigned k = 0; k < 5; ++k) { cp3d(i, j, k) = cp3[i][j][k]; }
        }
    }

    BSpline<double, 3> spline_3d_3(3, cp3d,
                                   make_pair(knots.begin(), knots.end()),
                                   make_pair(knots.begin(), knots.end()),
                                   make_pair(knots.begin(), knots.end()));

    std::cout << "\n3D B-Spline Test:\n";

    // some random points
    constexpr array<array<double, 3>, 20> coords_3d{
        {{0.6771169667756889, 0.9215782737838114, 0.4073093854149159},
         {0.5587048417879672, 0.27660875852915034, 0.480940182622974},
         {0.22166100919838327, 0.20423034083657687, 0.6828162120660075},
         {0.946117300279044, 0.6987537482071375, 0.988288195761676},
         {0.6657135834502097, 0.6362329447266486, 0.6098637234378605},
         {0.4831720672055162, 0.04449885480903437, 0.10794187260224075},
         {0.01970883141191182, 0.707507573262028, 0.3842127548987244},
         {0.2708700179218677, 0.9968009386232979, 0.7784159674453355},
         {0.34911366851070724, 0.37498001215861887, 0.5313335523213891},
         {0.7831337922295372, 0.3523807977438491, 0.5616443755053859},
         {0.5651837447297534, 0.3521803904983014, 0.3766317887491841},
         {0.42938234954177146, 0.03216063762167898, 0.2257110836573426},
         {0.5320732583223098, 0.2563335748568327, 0.4284182446681921},
         {0.27509099660894876, 0.3977734281932175, 0.7406749938043244},
         {0.2645390059243993, 0.8727863468285368, 0.8960358519530291},
         {0.8338988630403117, 0.2793369171836766, 0.6474522893994397},
         {0.01796573121296885, 0.8821327182314185, 0.023056755634403236},
         {0.29524574329222886, 0.7820003453674633, 0.1953107714234823},
         {0.27122873705368566, 0.050574347888162174, 0.8570032085001116},
         {0.3328181108995292, 0.22075349315084103, 0.7955065113795101}}};

    // and their spline value, pre-computed by MMA
    constexpr array<double, 20> vals_3d{
        {-0.1680401055752912,  -1.5184595156006049,  -0.6181298890131312,
         -0.7634699499420325,  -0.11071964593272185, 0.36484761815950045,
         -0.13606044311702267, -1.066384805143148,   -1.3030396801520705,
         -0.1366609994525347,  -1.845924629958203,   0.7750141662857806,
         -1.7374305394795615,  0.1278576545692187,   0.12587781246729338,
         -0.16274434441451202, -0.04873351906943385, 0.11126350145310085,
         -0.24228756896767534, -0.35239496660939207}};

    d = 0;
    for (unsigned i = 0; i < vals_3d.size(); ++i) {
        const double f =
            spline_3d_3(coords_3d[i][0], coords_3d[i][1], coords_3d[i][2]);
        d += (f - vals_3d[i]) * (f - vals_3d[i]);
        std::cout << coords_3d[i][0] << " " << coords_3d[i][1] << " "
                  << coords_3d[i][2] << " " << f << '\n';
    }

    assertion(std::sqrt(d / vals_3d.size()) < 1e-14);
    std::cout << "\n3D test "
              << (assertion.status() == 0 ? "succeed" : "failed") << '\n';

    // 2D with one dimension being periodic

    std::cout << "\n2D B-Spline with periodic boundary Test:\n";

    auto knots2 = {-0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6};
    BSpline<double, 2> spline_2d_3_periodic(
        3, {false, true}, cp2d, make_pair(knots.begin(), knots.end()),
        make_pair(knots2.begin(), knots2.end()));

    constexpr array<double, 10> vals_2d_periodic{
        -2.4226917831367354, -1.6666970964438268, 1.292431192225407,
        0.7418964927133986,  -1.6357978631810504, 0.32021767652671507,
        -0.0198531979338478, -1.6728476167273418, -2.414215484638793,
        0.7421859678077125};

    d = 0;
    for (unsigned i = 0; i < coords_2d.size(); ++i) {
        const double f =
            spline_2d_3_periodic(coords_2d[i].first, coords_2d[i].second);
        d += (f - vals_2d_periodic[i]) * (f - vals_2d_periodic[i]);
        std::cout << coords_2d[i].first << " " << coords_2d[i].second << " "
                  << f << '\n';
    }
    assertion(std::sqrt(d / coords_2d.size()) < 1e-14);
    std::cout << "\n2D test with periodic boundary "
              << (assertion.status() == 0 ? "succeed" : "failed") << '\n';

    return assertion.status();
}