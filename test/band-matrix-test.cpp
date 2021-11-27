#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>

#include "../src/include/BandMatrix.hpp"

using namespace std;

int main() {
    BandMatrix<double> mat{5, 1, 1};
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (i == j) mat(i, j) = -2;
            if (i - j == 1) mat(i, j) = 1;
            if (i - j == -1) mat(i, j) = 1;
        }
    }
    vector<double> b = {0.382505, 0.467217, 0.82818, -0.266586, -0.385533};
    auto x = mat.linear_solve(b);
    auto bb = mat * x;

    cout << "Band matrix solver test A.x = b (A models Dirichlet BC 1D Poisson "
            "equation):\nb = ";
    copy(b.begin(), b.end(), ostream_iterator<double>(cout, " "));
    cout << "\nx = ";
    copy(x.begin(), x.end(), ostream_iterator<double>(cout, " "));
    cout << "\nAnd A.x using the solution x:\nA.x = ";
    copy(bb.begin(), bb.end(), ostream_iterator<double>(cout, " "));

    double d = 0;
    for (int i = 0; i < b.size(); ++i) { d += (b[i] - bb[i]) * (b[i] - bb[i]); }
    d = std::sqrt(d / b.size());

    return d < 1e-10 ? 0 : 1;
}