#include <iterator>

#include "../src/include/BSpline.hpp"
#include "Assertion.hpp"

using namespace std;

int main() {
    Assertion assertion;
    Mesh<double, 3> mesh(10, 20, 30);

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

    return assertion.status();
}