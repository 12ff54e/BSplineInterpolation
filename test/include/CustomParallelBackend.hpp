#ifndef INTP_TEST_CUSTOM_PARALLEL_BACKEND
#define INTP_TEST_CUSTOM_PARALLEL_BACKEND

// Adapters may include the pool type before defining the range hook.
#include "bsintp/DedicatedThreadPool.hpp"

#include <algorithm>
#include <chrono>
#include <deque>
#include <exception>
#include <future>
#include <stdexcept>
#include <utility>
#include <vector>

#if !defined(INTP_MULTITHREAD) || !defined(INTP_CUSTOM_PARALLEL_BACKEND)
#error "This adapter requires both parallel configuration macros"
#endif

namespace parallel_test {

// A single worker executes queued tasks only while waiting. A nested spline
// construction would deadlock if a range join did not make progress.
class Backend {
   public:
    std::size_t submitted = 0;
    std::size_t completed = 0;
    std::size_t waited = 0;
    std::size_t fail_submission = 0;
    std::size_t fail_execution = 0;

    template <typename Func>
    void parallel_for(std::size_t count,
                      std::size_t grain,
                      const Func& function) {
        if (grain == 0) {
            throw std::invalid_argument("parallel_for grain must be positive");
        }
        if (count == 0) { return; }
        // Deliberately ignore the hint and queue every range, including small
        // ones. The spline must accept the backend's partition boundaries.
        const std::size_t chunk_size = 7;
        const auto task_num = 1 + (count - 1) / chunk_size;
        std::vector<std::future<void>> tasks;
        tasks.reserve(task_num);
        std::exception_ptr error;
        try {
            for (std::size_t i = 0; i < task_num; ++i) {
                const auto begin = i * chunk_size;
                const auto end = begin + std::min(chunk_size, count - begin);
                tasks.push_back(
                    submit([&function, begin, end] { function(begin, end); }));
            }
        } catch (...) { error = std::current_exception(); }
        for (auto& task : tasks) {
            try {
                wait(task);
            } catch (...) {
                if (!error) { error = std::current_exception(); }
            }
        }
        if (error) { std::rethrow_exception(error); }
    }

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
void parallel_for(std::size_t count, std::size_t grain, const Func& function) {
    parallel_test::backend().parallel_for(count, grain, function);
}

}  // namespace intp

#endif  // INTP_TEST_CUSTOM_PARALLEL_BACKEND
