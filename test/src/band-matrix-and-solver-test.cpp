#include <BandLU.hpp>
#include <BandMatrix.hpp>
#include "include/Assertion.hpp"
#include "include/rel_err.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

template <typename Mat, typename Vec>
void check_solver(Mat&& mat, const Vec& b, Assertion& assertion) {
    using namespace intp;

    // make a copy of b for later comparision, in case b is a iterator/pointer
    std::vector<double> b_vec;
    for (size_t i = 0; i < mat.dim(); ++i) { b_vec.push_back(b[i]); }

    BandLU<util::remove_cvref_t<Mat>> solver{mat};
    auto x = solver.solve(b);
    // matrix-vector multiplication only support container with index operator
    std::vector<double> x_vec;
    for (size_t i = 0; i < mat.dim(); ++i) { x_vec.push_back(x[i]); }
    auto bb_vec = mat * x_vec;

    double d = rel_err([](double x_) { return x_; },
                       std::make_pair(bb_vec.begin(), bb_vec.end()),
                       std::make_pair(b_vec.begin(), b_vec.end()));

    assertion(d < 1e-10);
    std::cout << "\n||b - A . x||/||b|| = " << d << '\n';
}

int main() {
    using namespace intp;

    Assertion assertion;

    // matrix dimension
    constexpr size_t n = 1 << 6;

    // rhs vector
    std::vector<double> b;
    b.reserve(n);
    {
        std::mt19937 rand_gen(
            static_cast<unsigned int>(std::chrono::high_resolution_clock::now()
                                          .time_since_epoch()
                                          .count()));
        std::uniform_real_distribution<> rand_dist(-1., 1.);
        for (size_t i = 0; i < n; ++i) { b.push_back(rand_dist(rand_gen)); }
    }

    // normal band matrix
    {
        // band matrix with lower and upper bandwidth = 1
        // Laplacian operator with zero boundary condition
        BandMatrix<double> mat{n, 1, 1};
        for (size_t i = 0; i < n; ++i) {
            mat(i, i) = -2;
            if (i > 0) mat(i, i - 1) = 1;
            if (i < n - 1) mat(i, i + 1) = 1;
        }
        std::cout << "\nBand matrix, Laplacian operator\n";
        check_solver(std::move(mat), b.data(), assertion);
    }

    // band matrix with extra side bands
    {
        // band matrix with lower and upper bandwidth = 1
        // Basic spline value of order 2
        ExtendedBandMatrix<double> mat1{n, 1, 1};
        for (size_t i = 0; i < n; ++i) {
            mat1(i, i) = 3. / 4.;
            mat1(i, (i + n - 1) % n) = 1. / 8.;
            mat1(i, (i + 1) % n) = 1. / 8.;
        }
        std::cout << "\nBand matrix, extended, 2nd order BSpline\n";
        check_solver(std::move(mat1), b, assertion);

        // band matrix with lower and upper bandwidth = 2
        // Basic spline value of order 4
        ExtendedBandMatrix<double> mat2{n, 2, 2};
        for (size_t i = 0; i < n; ++i) {
            mat2(i, i) = 115. / 192.;
            mat2(i, (i + n - 1) % n) = 19. / 96.;
            mat2(i, (i + 1) % n) = 19. / 96.;
            mat2(i, (i + n - 2) % n) = 1. / 384.;
            mat2(i, (i + 2) % n) = 1. / 384.;
        }
        std::cout << "\nBand matrix, extended, 4th order BSpline\n";
        check_solver(std::move(mat2), b, assertion);

        // band matrix with lower bandwidth = 1 and upper bandwidth = 3
        // Arbitrary values
        ExtendedBandMatrix<double> mat3{n, 1, 3};
        for (size_t i = 0; i < n; ++i) {
            mat3(i, (i + n - 1) % n) = 2889. / 16000.;
            mat3(i, i) = 1701. / 3200.;
            mat3(i, (i + 1) % n) = 33. / 128.;
            mat3(i, (i + 2) % n) = 729. / 20000.;
            mat3(i, (i + 3) % n) = 9. / 20000.;
        }
        std::cout << "\nBand matrix, extended, 4td order uneven BSpline\n";
        check_solver(std::move(mat3), b, assertion);
    }

    return assertion.status();
}
