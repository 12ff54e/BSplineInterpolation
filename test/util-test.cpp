#include <iostream>
#include <type_traits>
#include <vector>

#include "../src/include/util.hpp"
#include "Assertion.hpp"

int main() {
    using namespace intp;

    Assertion assertion;

    assertion(std::is_same<util::make_index_sequence<3>,
                           util::index_sequence<0u, 1u, 2u>>::value,
              "Issues on index_sequence.");

    assertion(util::pow(2, 3u) == 8, "Issues on power function.");

#ifdef STACK_ALLOCATOR
    constexpr unsigned N = 4;
    int buffer[N];
    util::stack_allocator<int, N> alloc(buffer);
    std::vector<int, util::stack_allocator<int, N>> vec({1, 2, 3, 4}, alloc);

    assertion(buffer[0] == 1);
    assertion(buffer[1] == 2);
    assertion(buffer[2] == 3);
    assertion(buffer[3] == 4);

    if (assertion.last_status() != 0) {
        std::cout << "Vector is not allocated on the given allocator.\n";
    }

    vec.resize(2);

    assertion(buffer[0] == 1);
    assertion(buffer[1] == 2);

    if (assertion.last_status() != 0) {
        std::cout << "Deallocate do not work properly.\n";
    }

    vec.push_back(5);
    vec.push_back(7);

    bool exception_thrown = false;
    try {
        vec.push_back(42);
    } catch (const std::exception& e) { exception_thrown = true; }
    assertion(exception_thrown,
              "No exception is thrown when allocator capacity is exceeded.\n");
#endif

    assertion(util::is_iterable<std::vector<int>>::value);
    assertion(!util::is_iterable<double>::value);

    assertion(util::is_indexed<std::vector<int>>::value);
    assertion(!util::is_indexed<std::initializer_list<int>>::value);

    return assertion.status();
}
