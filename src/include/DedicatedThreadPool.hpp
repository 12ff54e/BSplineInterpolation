#pragma once

#include <functional>
#include <future>
#include <queue>
#include <thread>

namespace intp {

/**
 * @brief A thread pool managing a bunch of threads and a queue of tasks. It can
 * spawn a specified number of threads during construction and joins all the
 * threads in destruction. During its lifetime, it can accepts and dispatch
 * tasks to those threads.
 *
 * @tparam T The return type of tasks
 */
template <typename T>
class DedicatedThreadPool {
   public:
    using return_type = T;

    /**
     * @brief Construct a new Thread Pool
     *
     * @param num Thread number in the pool
     */
    DedicatedThreadPool(std::size_t num = std::thread::hardware_concurrency())
        : thread_num(num == 0 ? DEFAULT_THREAD_NUM : num),
          join_threads(threads) {
        for (std::size_t i = 0; i < thread_num; i++) {
            threads.emplace_back(&DedicatedThreadPool::thread_loop, this);
        }
#ifdef _DEBUG
        std::cout << "[DEBUG] Thread pool initialized with " << num
                  << " threads.\n";
#endif
    }

    ~DedicatedThreadPool() {
        {
            std::unique_lock<std::mutex> lk(queue_mutex);
            should_terminate = true;
        }
        // notify all threads to exit working loop
        cv.notify_all();
    }

    /**
     * @brief Add a task in queue, the Func type should not have any parameters.
     * Return a future containing function return value.
     */
    template <typename Func>
    std::future<return_type> queue_task(Func func) {
        std::packaged_task<return_type()> task(std::move(func));
        auto res = task.get_future();
        {
            std::unique_lock<std::mutex> lk(queue_mutex);
            tasks.push(std::move(task));
        }
        cv.notify_one();
        return res;
    }

    /**
     * @brief Return true if there are no tasks in queue. (But there may be
     * tasks executing by threads.)
     */
    bool is_queue_empty() {
        bool queue_empty;
        {
            std::unique_lock<std::mutex> lk(queue_mutex);
            queue_empty = tasks.empty();
        }
        return queue_empty;
    }

    /**
     * @brief Singleton style instance getter
     *
     * @param num Thread number in the pool, default to be
     * hardware_concurrency()
     */
    static DedicatedThreadPool& get_instance(
        std::size_t num = std::thread::hardware_concurrency()) {
        static DedicatedThreadPool thread_pool(num);
        return thread_pool;
    }

   private:
    static constexpr std::size_t DEFAULT_THREAD_NUM = 8;

    /**
     * @brief Hold a ref of thread vector, use RAII to ensure all the threads is
     * joined.
     *
     */
    struct JoinThreads {
        std::vector<std::thread>& ts_;

        JoinThreads(std::vector<std::thread>& ts) : ts_(ts) {}
        JoinThreads(const JoinThreads&) = delete;
        JoinThreads(JoinThreads&) = delete;
        JoinThreads& operator=(JoinThreads&) = delete;
        ~JoinThreads() {
            for (auto& t : ts_) {
                if (t.joinable()) { t.join(); }
            }
        }
    };

    /**
     * @brief  scheduler function
     *
     */
    void thread_loop() {
        while (true) {
            std::packaged_task<T()> task;
            {
                std::unique_lock<std::mutex> lk(queue_mutex);
                // Wait until there are tasks in queue (awaked by notify_one in
                // queue_task), or the thread pool is being shutdown
                // (awaked by notify_all in destructor)
                cv.wait(lk,
                        [this] { return !tasks.empty() || should_terminate; });
                if (should_terminate) { return; }

                // take a task from queue
                task = std::move(tasks.front());
                tasks.pop();
            }
            task();
        }
    }

    std::size_t thread_num;         // thread number in the pool
    bool should_terminate = false;  // Tells threads to stop looking for tasks
    std::mutex queue_mutex;  // Prevents data races to the task queue, all ops
                             // on tasks queue should be performed after locking
                             // it. This mutex is not needed with the presence
                             // of a thread safe queue for tasks.
    std::condition_variable
        cv;  // Allows threads to wait on new tasks or termination
    std::vector<std::thread> threads;                     // Thread container
    std::queue<std::packaged_task<return_type()>> tasks;  // Queue for tasks
    JoinThreads join_threads;  // Defined last to ensure destruct first
};

}  // namespace intp
