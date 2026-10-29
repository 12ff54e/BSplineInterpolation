#include "bsintp/Parallel.hpp"

#include "Assertion.hpp"

#include <atomic>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef _MSC_VER
static_assert(__cplusplus == 201103L, "This test must compile in C++11 mode");
#endif

namespace {

void check_ranges(Assertion& assertion) {
    for (std::size_t count : {0u, 1u, 17u, 1003u}) {
        for (std::size_t grain : {1u, 5u, 64u, 2000u}) {
            std::vector<std::atomic<int>> visits(count);
            for (auto& visit : visits) { visit.store(0); }
            std::atomic<bool> valid(true);
            intp::parallel_for(count, grain,
                               [&](std::size_t begin, std::size_t end) {
                                   if (begin >= end || end > count) {
                                       valid = false;
                                       return;
                                   }
                                   for (auto i = begin; i < end; ++i) {
                                       ++visits[i];
                                   }
                               });
            assertion(valid && std::all_of(visits.begin(), visits.end(),
                                           [](const std::atomic<int>& n) {
                                               return n.load() == 1;
                                           }),
                      "Every index must complete exactly once");
        }
    }
    bool caught = false;
    try {
        intp::parallel_for(1, 0, [](std::size_t, std::size_t) {});
    } catch (const std::invalid_argument&) { caught = true; }
    assertion(caught, "Zero grain must be rejected");
}

void check_failure(bool caller, Assertion& assertion) {
    std::atomic<std::size_t> completed(0);
    std::atomic<bool> alive(true);
    bool caught = false;
    try {
        const auto buffer = std::make_shared<int>(42);
        const std::weak_ptr<int> borrowed = buffer;
        intp::parallel_for(129, 8,
                           [&, borrowed](std::size_t begin, std::size_t end) {
                               if (begin == (caller ? 0u : 8u)) {
                                   throw std::runtime_error("callback failure");
                               }
                               const auto live = borrowed.lock();
                               if (!live || *live != 42) { alive = false; }
                               completed += end - begin;
                           });
    } catch (const std::runtime_error& error) {
        caught = std::string(error.what()) == "callback failure";
    }
    assertion(caught && alive && completed == 121,
              "Caller and worker failures must drain all accepted work");
}

}  // namespace

int main() {
    Assertion assertion;
    check_ranges(assertion);
    check_failure(true, assertion);
    check_failure(false, assertion);
    return assertion.status();
}
