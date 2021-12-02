#include <utility>

#include "../src/include/BSpline.hpp"

using namespace std;

int main() {
    int err = 0;
    Mesh<double, 3> mesh(10, 20, 30);

    err = err == 0 && mesh.dim_acc_size()[0] == 1 ? 0 : 1;
    err = err == 0 && mesh.dim_acc_size()[1] == 30 ? 0 : 1;
    err = err == 0 && mesh.dim_acc_size()[2] == 20 * 30 ? 0 : 1;
    err = err == 0 && mesh.dim_acc_size()[3] == 10 * 20 * 30 ? 0 : 1;

    err = err == 0 && mesh.dim_size(0) == 10 ? 0 : 1;
    err = err == 0 && mesh.dim_size(1) == 20 ? 0 : 1;
    err = err == 0 && mesh.dim_size(2) == 30 ? 0 : 1;

    err = err == 0 && mesh.size() == 10 * 20 * 30 ? 0 : 1;

    mesh(3, 4, 5) = 1;
    err = err == 0 && *(mesh.data() + 3 * 600 + 4 * 30 + 5) == 1. ? 0 : 1;
    return err;
}