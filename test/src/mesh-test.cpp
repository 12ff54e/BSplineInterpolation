#include <iterator>

#include <Mesh.hpp>
#include "include/Assertion.hpp"

int main() {
    using namespace std;
    using namespace intp;

    Assertion assertion;
    Mesh<int, 3> mesh(10, 20, 30);

    assertion(mesh.dim_size(0) == 10, "Dimension wise size wrong.");
    assertion(mesh.dim_size(1) == 20, "Dimension wise size wrong.");
    assertion(mesh.dim_size(2) == 30, "Dimension wise size wrong.");

    assertion(mesh.size() == 10 * 20 * 30, "Total size wrong.");

    mesh(3, 4, 5) = 1;
    assertion(*(mesh.data() + 3 * 600 + 4 * 30 + 5) == 1.,
              "Modify data by index failed.");

#ifdef INTP_ENABLE_ASSERTION
    bool f{};
    try {
        mesh(3, 25, 5);
    } catch (std::exception& e) { f = true; }
    assertion(f, "Out of boundary access detection failed.");
#endif  // INTP_ENABLE_ASSERTION

    auto mesh_it = mesh.begin();
    advance(mesh_it, 3 * 600 + 4 * 30 + 5);
    auto indices = mesh.iter_indices(mesh_it);

    assertion(indices[0] == 3, "Iterator indexing failed.");
    assertion(indices[1] == 4, "Iterator indexing failed.");
    assertion(indices[2] == 5, "Iterator indexing failed.");

    mesh.resize({30, 20, 5});

    // test iterator along one dimension

    auto dim_it = mesh.begin(0, {0, 6, 2});
    for (auto it = dim_it; it != mesh.end(0, {0, 6, 2}); ++it) {
        *it = static_cast<int>(it - dim_it);
    }

    dim_it = mesh.begin(1, {4, 0, 4});
    for (auto it = dim_it; it != mesh.end(1, {4, 0, 4}); ++it) {
        *it = static_cast<int>(it - dim_it);
    }

    dim_it = mesh.begin(2, {9, 7, 0});
    for (auto it = dim_it; it != mesh.end(2, {9, 7, 0}); ++it) {
        *it = static_cast<int>(it - dim_it);
    }

    for (std::size_t dim_ind = 0; dim_ind < 3; ++dim_ind) {
        for (std::size_t i = 0; i < mesh.dim_size(dim_ind); ++i) {
            std::size_t val;
            if (dim_ind == 0) {
                assertion((val = static_cast<std::size_t>(mesh(i, 6, 2))) == i);
            } else if (dim_ind == 1) {
                assertion((val = static_cast<std::size_t>(mesh(4, i, 4))) == i);
            } else {
                assertion((val = static_cast<std::size_t>(mesh(9, 7, i))) == i);
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

    std::vector<int> vec{1, 1, 2, 3, 5, 8, 13, 21};
    Mesh<int, 1> mesh_1d_from_iterator(std::make_pair(vec.begin(), vec.end()));

    for (unsigned i = 0; i < vec.size(); ++i) {
        assertion(mesh_1d_from_iterator(i) == vec[i]);
        if (assertion.last_status() != 0) {
            std::cout << "Mesh 1D from container test failed.\n";
        }
    }

    // test equal-length-on-each-dimension constructor

    Mesh<int, 4> mesh_4d(5);
    assertion(mesh_4d.size() == util::pow(5u, 4u),
              "Equal length on each dimension mesh size wrong.");

    return assertion.status();
}
