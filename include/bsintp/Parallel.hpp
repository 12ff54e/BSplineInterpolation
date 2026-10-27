#ifndef INTP_PARALLEL
#define INTP_PARALLEL

#include "DedicatedThreadPool.hpp"

#include <exception>
#include <future>
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
        tasks_.push_back(intp::parallel_submit(std::forward<Func>(task)));
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
                intp::parallel_wait(task);
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
}  // namespace intp

#endif  // INTP_PARALLEL
