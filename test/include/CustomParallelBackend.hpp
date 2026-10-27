#ifndef INTP_TEST_CUSTOM_PARALLEL_BACKEND
#define INTP_TEST_CUSTOM_PARALLEL_BACKEND

// This include order lets adapters reuse the pool type before defining hooks.
#include "bsintp/DedicatedThreadPool.hpp"

#include <chrono>
#include <deque>
#include <future>
#include <stdexcept>
#include <utility>

#if !defined(INTP_MULTITHREAD) || !defined(INTP_CUSTOM_PARALLEL_BACKEND)
#error "This adapter requires both parallel configuration macros"
#endif

namespace parallel_test {

// A single worker executes queued tasks only while waiting. A nested spline
// construction would deadlock if any library join bypassed the wait hook.
class Backend {
   public:
    std::size_t submitted = 0;
    std::size_t completed = 0;
    std::size_t waited = 0;
    std::size_t fail_submission = 0;
    std::size_t fail_execution = 0;

    template <typename Func>
    std::future<void> submit(Func&& func) {
        if (fail_submission == submitted + 1) {
            throw std::runtime_error("submission failure");
        }
        std::packaged_task<void()> task;
        if (fail_execution == submitted + 1) {
            task = std::packaged_task<void()>(
                [] { throw std::runtime_error("task failure"); });
        } else {
            task = std::packaged_task<void()>(std::forward<Func>(func));
        }
        auto result = task.get_future();
        tasks_.push_back(std::move(task));
        ++submitted;
        return result;
    }

    void wait(std::future<void>& result) {
        ++waited;
        while (result.wait_for(std::chrono::seconds(0)) !=
               std::future_status::ready) {
            if (tasks_.empty()) {
                throw std::logic_error("no task can complete this future");
            }
            auto task = std::move(tasks_.front());
            tasks_.pop_front();
            task();
            ++completed;
        }
        result.get();
    }

    bool empty() const { return tasks_.empty(); }

    void reset() {
        if (!empty() || submitted != completed) {
            throw std::logic_error("reset while tasks are unfinished");
        }
        submitted = completed = waited = 0;
        fail_submission = fail_execution = 0;
    }

   private:
    std::deque<std::packaged_task<void()>> tasks_;
};

inline Backend& backend() {
    static thread_local Backend instance;
    return instance;
}

}  // namespace parallel_test

namespace intp {

template <typename Func>
std::future<void> parallel_submit(Func&& task) {
    return parallel_test::backend().submit(std::forward<Func>(task));
}

inline void parallel_wait(std::future<void>& task) {
    parallel_test::backend().wait(task);
}

}  // namespace intp

#endif  // INTP_TEST_CUSTOM_PARALLEL_BACKEND
