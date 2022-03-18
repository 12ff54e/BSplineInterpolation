#include <iterator>

#include "../src/include/Mesh.hpp"
#include "Assertion.hpp"

using namespace std;

int main() {
    Assertion assertion;
    Mesh<int, 3> mesh{10, 20, 30};

    assertion(mesh.dim_acc_size()[0] == 1);
    assertion(mesh.dim_acc_size()[1] == 30);
    assertion(mesh.dim_acc_size()[2] == 20 * 30);
    assertion(mesh.dim_acc_size()[3] == 10 * 20 * 30);

    assertion(mesh.dim_size(0) == 10);
    assertion(mesh.dim_size(1) == 20);
    assertion(mesh.dim_size(2) == 30);

    assertion(mesh.size() == 10 * 20 * 30);

    mesh(3, 4, 5) = 1;
    assertion(*(mesh.data() + 3 * 600 + 4 * 30 + 5) == 1.);

    auto it = mesh.begin();
    advance(it, 3 * 600 + 4 * 30 + 5);
    auto indices = mesh.iter_indices(it);

    assertion(indices[0] == 3);
    assertion(indices[1] == 4);
    assertion(indices[2] == 5);

    mesh.resize({30, 20, 5});

    assertion(mesh.dim_acc_size()[0] == 1);
    assertion(mesh.dim_acc_size()[1] == 5);
    assertion(mesh.dim_acc_size()[2] == 20 * 5);
    assertion(mesh.dim_acc_size()[3] == 30 * 20 * 5);

    // test iterator along one dimension

    auto dim_it = mesh.begin(0, {0, 6, 2});
    for (auto it = dim_it; it != mesh.end(0, {0, 6, 2}); ++it) {
        *it = it - dim_it;
    }

    dim_it = mesh.begin(1, {4, 0, 4});
    for (auto it = dim_it; it != mesh.end(1, {4, 0, 4}); ++it) {
        *it = it - dim_it;
    }

    dim_it = mesh.begin(2, {9, 7, 0});
    for (auto it = dim_it; it != mesh.end(2, {9, 7, 0}); ++it) {
        *it = it - dim_it;
    }

    for (int dim_ind = 0; dim_ind < 3; ++dim_ind) {
        for (unsigned i = 0; i < mesh.dim_size(dim_ind); ++i) {
            int val;
            if (dim_ind == 0) {
                assertion((val = mesh(i, 6, 2)) == i);
            } else if (dim_ind == 1) {
                assertion((val = mesh(4, i, 4)) == i);
            } else {
                assertion((val = mesh(9, 7, i)) == i);
            }
            if (assertion.last_status() != 0) {
                std::cout << "Dimension-wise iterator fails at dimension-"
                          << dim_ind << ", index = " << i
                          << ", and mesh value = " << val << '\n';
            }
        }
    }

    const Mesh<int, 3>& mesh_const_ref = mesh;
    auto iter2 = mesh_const_ref.begin(1, {4, 0, 4});
    iter2 += 5;
    assertion(*(iter2 - 4) == 1, "Dimension-wise iterator operator-(it, n)");
    assertion(*(iter2 + 4) == 9, "Dimension-wise iterator operator+(it, n)");
    assertion(*(4 + iter2) == 9, "Dimension-wise iterator operator+(n, it)");
    assertion(iter2[2] == 7, "Dimension-wise iterator operator[]");
    assertion(iter2[-3] == 2, "Dimension-wise iterator operator[]");

    // test special constructor for 1D case

    std::vector<double> vec{1, 1, 2, 3, 5, 8, 13, 21};
    Mesh<double, 1> mesh_1d_from_container(vec);
    Mesh<double, 1> mesh_1d_from_iterator(vec.begin(), vec.end());

    for (unsigned i = 0; i < vec.size(); ++i) {
        assertion(mesh_1d_from_container(i) == vec[i] &&
                  mesh_1d_from_iterator(i) == vec[i]);
        if (assertion.last_status() != 0) {
            std::cout << "Mesh 1D from container test failed.\n";
        }
    }

    // test equal-length-on-each-dimension constructor

    Mesh<int, 4> mesh_4d(5);
    assertion(mesh_4d.size() == util::pow(5, 4u));

    return assertion.status();
}