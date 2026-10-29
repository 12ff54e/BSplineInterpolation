#ifndef INTP_PARALLEL
#define INTP_PARALLEL

#ifndef INTP_CUSTOM_PARALLEL_BACKEND
#include "DedicatedThreadPool.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <future>
#include <stdexcept>
#include <utility>
#include <vector>

namespace intp {
namespace detail {

// Keep accepted tasks within the lifetime of their borrowed buffers, including
// failures in submission or in work performed by the calling thread.
class ParallelTasks {
   public:
    explicit ParallelTasks(std::size_t count) { tasks_.reserve(count); }

    ParallelTasks(const ParallelTasks&) = delete;
    ParallelTasks& operator=(const ParallelTasks&) = delete;

    ~ParallelTasks() { drain(); }

    template <typename Func>
    void submit(Func&& task) {
        // Reserve before submission so allocation cannot lose a task's handle.
        if (tasks_.size() == tasks_.capacity()) {
            tasks_.reserve(2 * tasks_.size() + 1);
        }
        tasks_.push_back(DedicatedThreadPool<void>::get_instance(8).queue_task(
            std::forward<Func>(task)));
    }

    void wait() {
        const auto error = drain();
        if (error) { std::rethrow_exception(error); }
    }

   private:
    std::exception_ptr drain() noexcept {
        std::exception_ptr error;
        for (auto& task : tasks_) {
            try {
                task.get();
            } catch (...) {
                if (!error) { error = std::current_exception(); }
            }
        }
        tasks_.clear();
        return error;
    }

    std::vector<std::future<void>> tasks_;
};

}  // namespace detail

// Process [0, count) in disjoint ranges. Finish all accepted work before
// returning or throwing, while the caller's buffers and function remain alive.
template <typename Func>
void parallel_for(std::size_t count, std::size_t grain, const Func& function) {
    if (grain == 0) {
        throw std::invalid_argument("parallel_for grain must be positive");
    }
    if (count == 0) { return; }
    const auto task_num = 1 + (count - 1) / grain;
    detail::ParallelTasks tasks(task_num - 1);
    for (std::size_t i = 1; i < task_num; ++i) {
        const auto begin = i * grain;
        const auto end = begin + std::min(grain, count - begin);
        tasks.submit([&function, begin, end] { function(begin, end); });
    }
    function(0, std::min(grain, count));
    tasks.wait();
}

}  // namespace intp
#endif  // INTP_CUSTOM_PARALLEL_BACKEND

#endif  // INTP_PARALLEL
