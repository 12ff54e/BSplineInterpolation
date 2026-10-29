#ifdef INTP_ALLOCATOR_TEST_SERIAL
#undef INTP_MULTITHREAD
#endif
#define INTP_CUSTOM_ALLOCATOR

#include "bsintp/aligned-allocator.hpp"

#include <algorithm>
#include <cstdint>
#include <future>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace {

void check(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}

struct Resource {
    std::mutex mutex;
    std::map<void*, std::size_t> live;
    std::weak_ptr<Resource> copy_target;
    int fail_after = -1;
    std::size_t allocations = 0;

    ~Resource() {
        if (!live.empty()) { std::terminate(); }
    }
};

thread_local std::shared_ptr<Resource> active;

template <typename T>
class TestAllocator {
   public:
    using value_type = T;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_copy_assignment = std::false_type;
    using propagate_on_container_swap = std::true_type;

    template <typename U>
    struct rebind {
        using other = TestAllocator<U>;
    };

    TestAllocator() : resource(active) {
        if (!resource) { resource = std::make_shared<Resource>(); }
    }
    explicit TestAllocator(std::shared_ptr<Resource> owner)
        : resource(std::move(owner)) {}
    TestAllocator(const TestAllocator&) = default;
    TestAllocator& operator=(const TestAllocator&) = default;
    template <typename U>
    TestAllocator(const TestAllocator<U>& other) noexcept
        : resource(other.resource) {}

    TestAllocator select_on_container_copy_construction() const {
        const auto target = resource->copy_target.lock();
        return TestAllocator(target ? target : resource);
    }

    T* allocate(std::size_t count) {
        if (count > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
            throw std::bad_alloc();
        }
        const std::lock_guard<std::mutex> lock(resource->mutex);
        if (resource->fail_after == 0) { throw std::bad_alloc(); }
        if (resource->fail_after > 0) { --resource->fail_after; }
        const auto alignment =
            std::max(std::size_t{32}, std::size_t{alignof(T)});
        auto* ptr = static_cast<T*>(intp::detail::allocate_aligned_memory(
            alignment, count * sizeof(T)));
        if (!ptr) { throw std::bad_alloc(); }
        try {
            resource->live.emplace(ptr, count * sizeof(T));
        } catch (...) {
            intp::detail::deallocate_aligned_memory(ptr);
            throw;
        }
        ++resource->allocations;
        return ptr;
    }

    void deallocate(T* ptr, std::size_t count) noexcept {
        const std::lock_guard<std::mutex> lock(resource->mutex);
        const auto found = resource->live.find(ptr);
        if (found == resource->live.end() ||
            found->second != count * sizeof(T)) {
            std::terminate();
        }
        resource->live.erase(found);
        intp::detail::deallocate_aligned_memory(ptr);
    }

    std::shared_ptr<Resource> resource;
};

template <typename T, typename U>
bool operator==(const TestAllocator<T>& a, const TestAllocator<U>& b) noexcept {
    return a.resource == b.resource;
}
template <typename T, typename U>
bool operator!=(const TestAllocator<T>& a, const TestAllocator<U>& b) noexcept {
    return !(a == b);
}

}  // namespace

namespace intp {
template <typename T>
using BufferAllocator = TestAllocator<T>;
}  // namespace intp

#include "bsintp/Interpolation.hpp"

#include <cmath>
#include <iostream>

namespace {

using Function = intp::InterpolationFunction<double, 2, 3>;

Function make_function() {
    intp::Mesh<double, 2> values(13, 15);
    for (std::size_t i = 0; i < 13; ++i) {
        for (std::size_t j = 0; j < 15; ++j) {
            const double x = i / 12., y = j / 14.;
            values(i, j) = x * x * x + 2. * y * y;
        }
    }
    return Function(values, std::make_pair(0., 1.), std::make_pair(0., 1.));
}

void check_value(const Function& function) {
    const double x = .37, y = .61;
    check(std::abs(function(x, y) - (x * x * x + 2. * y * y)) < 1e-11,
          "custom allocation changed the interpolated value");
    check(std::abs(function.derivative_at({{x, y}}, {{1, 1}})) < 1e-10,
          "custom allocation changed the mixed derivative");
}

void check_lifetime() {
    auto first = std::make_shared<Resource>();
    auto second = std::make_shared<Resource>();
    active = first;
    {
        auto function = make_function();
        check(first->allocations >= 3, "solve buffers did not use the hook");
        check(reinterpret_cast<std::uintptr_t>(
                  function.spline().control_points().data()) %
                      32 ==
                  0,
              "control points lost their alignment");
        check_value(function);
        active = second;
        function.generate_polynomial_form();
        check(second->allocations > 0, "polynomials did not use the hook");
        check(reinterpret_cast<std::uintptr_t>(
                  function.spline().polynomial_coef().data()) %
                      32 ==
                  0,
              "polynomial coefficients lost their alignment");
        function.release_control_points();
        check(first->live.empty(), "release used the current resource");
        check_value(function);

        // Exercise the adaptor, including a base allocator that returns a
        // different resource from select_on_container_copy_construction().
        second->copy_target = first;
        auto copy = function;
        check(!first->live.empty(), "copy selection lost allocator state");
        copy.scale_coefficients(2.);
        check_value(function);
        Function moved;
        moved = std::move(function);
        check_value(moved);
        active.reset();
        function = copy;
        check(std::abs(function(.37, .61) - copy(.37, .61)) < 1e-11,
              "copy assignment lost coefficients");
    }
    check(first->live.empty() && second->live.empty(), "allocation leaked");

    using Alloc =
        intp::util::default_init_allocator<double, TestAllocator<double>>;
    const Alloc alloc{TestAllocator<double>(first)};
    std::allocator_traits<Alloc>::rebind_alloc<int> rebound(alloc);
    check(rebound.resource == first, "rebind lost the resource");

    std::weak_ptr<Resource> owner;
    std::unique_ptr<Function> retained;
    {
        active = std::make_shared<Resource>();
        owner = active;
        retained.reset(new Function(make_function()));
        retained->generate_polynomial_form();
        retained->release_control_points();
        active.reset();
    }
    check(!owner.expired(), "resource died before its retained coefficients");
    std::thread destroy([&] {
        check_value(*retained);
        retained.reset();
    });
    destroy.join();
    check(owner.expired(), "cross-thread destruction retained the resource");
}

void check_failures() {
    active = std::make_shared<Resource>();
    bool completed = false;
    for (int limit = 0; limit < 16 && !completed; ++limit) {
        active->fail_after = limit;
        try {
            auto function = make_function();
            function.generate_polynomial_form();
            function.release_control_points();
            check_value(function);
            completed = true;
        } catch (const std::bad_alloc&) {}
        check(active->live.empty(), "failed construction leaked a buffer");
    }
    check(completed, "failure injection never reached successful construction");
    active->fail_after = -1;
    {
        auto function = make_function();
        active->fail_after = 0;
        bool failed = false;
        try {
            function.generate_polynomial_form();
        } catch (const std::bad_alloc&) { failed = true; }
        check(failed && !function.spline().has_polynomial(),
              "polynomial allocation failure was not propagated");
        check_value(function);
    }
    check(active->live.empty(), "failed polynomial construction leaked");
    active.reset();
}

void check_dimensions() {
    active = std::make_shared<Resource>();
    {
        std::vector<double> line{0., 1., 8., 27., 64.};
        intp::InterpolationFunction<double, 1, 3> one(
            intp::util::get_range(line), std::make_pair(0., 4.));
        intp::Mesh<double, 3> values(5, 6, 7);
        for (std::size_t i = 0; i < 5; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                for (std::size_t k = 0; k < 7; ++k) {
                    values(i, j, k) = i * i * i + 2. * j * j + k;
                }
            }
        }
        intp::InterpolationFunction<double, 3, 3> three(
            values, std::make_pair(0., 4.), std::make_pair(0., 5.),
            std::make_pair(0., 6.));
        for (bool polynomial : {false, true}) {
            if (polynomial) {
                one.generate_polynomial_form();
                three.generate_polynomial_form();
                one.release_control_points();
                three.release_control_points();
            }
            check(std::abs(one(.5) - .125) < 1e-11 &&
                      std::abs(three(.5, 1., 2.) - 4.125) < 1e-11,
                  "odd-dimensional solve lost its result buffer");
        }
    }
    check(active->live.empty(), "odd-dimensional construction leaked");
    active.reset();
}

void check_alignment() {
    struct alignas(64) Value {
        double value;
    };
    using Alloc =
        intp::util::default_init_allocator<Value, TestAllocator<Value>>;
    active = std::make_shared<Resource>();
    {
        std::vector<Value, Alloc> values(4);
        check(
            reinterpret_cast<std::uintptr_t>(values.data()) % alignof(Value) ==
                0,
            "over-aligned values lost their alignment");
        Alloc alloc;
        bool rejected = false;
        try {
            alloc.allocate(std::numeric_limits<std::size_t>::max());
        } catch (const std::bad_alloc&) { rejected = true; }
        check(rejected, "the custom allocator did not reject size overflow");
    }
    check(active->live.empty(), "over-aligned allocation leaked");
    active.reset();
}

}  // namespace

int main() {
    try {
        check_lifetime();
        check_failures();
        check_dimensions();
        check_alignment();
        // Independent callers must retain their own allocator identities.
        auto first = std::async(std::launch::async, check_lifetime);
        auto second = std::async(std::launch::async, check_lifetime);
        first.get();
        second.get();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
