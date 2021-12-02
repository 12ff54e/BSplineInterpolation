#include <iostream>
#include <type_traits>

#include "../src/include/util.hpp"

int main() {
    int err = 0;

    err = err == 0 && std::is_same<util::make_index_sequence<3>,
                                   util::index_sequence<0u, 1u, 2u>>::value
              ? 0
              : 1;

    err = err == 0 && util::pow(2, 3u) == 8 ? 0 : 1;

    return err;
}