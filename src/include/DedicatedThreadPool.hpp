#ifndef INTP_THREAD_POOL
#define INTP_THREAD_POOL

#include <future>   //unique_lock, packaged_task
#include <mutex>    // mutex
#include <queue>    // deque
#include <thread>   // hardware_concurrency
#include <utility>  // move

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
   private:
    using lock_type = std::unique_lock<std::mutex>;
    using task_type = std::packaged_task<T()>;

    /**
     * @brief A queue of tasks that can be stolen from. Tasks are pushed and
     * popped in back, and being stolen from front.
     *
     */
    struct shared_working_queue {
       public:
#ifdef INTP_DEBUG
        size_t submitted{};
        size_t executed{};
        size_t stolen{};
        size_t stealing{};
#endif

        shared_working_queue() = default;
        void push(task_type&& task) {
            lock_type lk(deq_mutex);
#ifdef INTP_DEBUG
            ++submitted;
#endif
            deq.push_back(std::move(task));
        }
        bool try_pop(task_type& task) {
            lock_type lk(deq_mutex);
            if (deq.empty()) { return false; }
            task = std::move(deq.back());
            deq.pop_back();
#ifdef INTP_DEBUG
            ++executed;
#endif
            return true;
        }
        bool try_steal(task_type& task) {
            lock_type lk(deq_mutex);
            if (deq.empty()) { return false; }
            task = std::move(deq.front());
            deq.pop_front();
#ifdef INTP_DEBUG
            ++stolen;
#endif
            return true;
        }
        bool empty() const {
            bool empty_;
            {
                lock_type lk(deq_mutex);
                empty_ = deq.empty();
            }
            return empty_;
        }

       private:
        std::deque<task_type> deq;
        mutable std::mutex deq_mutex;
    };

    /**
     * @brief Construct a new Thread Pool
     *
     * @param num Thread number in the pool
     */
    DedicatedThreadPool(size_t num = std::thread::hardware_concurrency())
        : join_threads(threads) {
        auto t_num = num == 0 ? DEFAULT_THREAD_NUM : num;
        try {
            for (size_t i = 0; i < t_num; i++) {
                worker_queues.emplace_back(new shared_working_queue{});
            }
            for (size_t i = 0; i < t_num; i++) {
                threads.emplace_back(&DedicatedThreadPool::thread_loop, this,
                                     i);
            }
#ifdef INTP_DEBUG
            std::cout << "[DEBUG] Thread pool initialized with " << t_num
                      << " threads.\n";
#endif
        } catch (...) {  // in this case dtor is not called, so threads
                         // termination should be proper handled here
            should_terminate = true;
            cv.notify_all();
            throw;
        }
    }

   public:
    using return_type = T;

    ~DedicatedThreadPool() {
        should_terminate = true;
        cv.notify_all();

#ifdef INTP_DEBUG
        std::cout << "\n[DEBUG] The thread pool has " << thread_num()
                  << " threads.\n"
                  << "[DEBUG]    Main queue has " << submitted
                  << " submission.\n"
                  << "[DEBUG] Worker queue stats:\n"
                  << "[DEBUG]\t\tSubmitted  Executed  Stolen  Stealing\n";
        for (size_t i = 0; i < worker_queues.size(); ++i) {
            auto ptr = worker_queues[i].get();
            std::cout << "[DEBUG] Thread " << i << ": \t" << ptr->submitted
                      << '\t' << ptr->executed + ptr->stealing << '\t'
                      << ptr->stolen << '\t' << ptr->stealing << '\n';
        }
        std::cout << "[DEBUG] NOTE: Submitted = Executed - Stolen + Stealing\n";
#endif
    }

    /**
     * @brief Add a task in queue, the Func type should not have any parameters
     * (use lambda capture). Return a future containing function return value.
     */
    template <typename Func>
    std::future<return_type> queue_task(Func func) {
        task_type task(std::move(func));
        auto res = task.get_future();
        {
            lock_type lk(main_queue_mutex);
            main_queue.push(std::move(task));
#ifdef INTP_DEBUG
            ++submitted;
#endif
        }
        cv.notify_one();
        return res;
    }

    /**
     * @brief Return true if there are no tasks in queue. (But there may be
     * tasks being executing by threads.)
     */
    bool is_main_queue_empty() {
        bool queue_empty;
        {
            lock_type lk(main_queue_mutex);
            queue_empty = main_queue.empty();
        }
        return queue_empty;
    }

    size_t thread_num() const {
        return threads.size();
    }

    /**
     * @brief Singleton style instance getter
     *
     * @param num Thread number in the pool, default to be
     * hardware_concurrency()
     */
    static DedicatedThreadPool& get_instance(
        size_t num = std::thread::hardware_concurrency()) {
        static DedicatedThreadPool thread_pool(num);
        return thread_pool;
    }

   private:
    static constexpr size_t DEFAULT_THREAD_NUM = 8;
    static constexpr size_t BATCH_SIZE = 16;

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

    bool try_pop_from_main(task_type& task) {
        {
            lock_type lk(main_queue_mutex);
            if (main_queue.empty()) { return false; }
            size_t c{};
            while (c++ < BATCH_SIZE && !main_queue.empty()) {
                worker_queue_ptr->push(std::move(main_queue.front()));
                main_queue.pop();
            }
        }
        return worker_queue_ptr->try_pop(task);
    }

    bool try_steal_from_others(task_type& task) {
        for (size_t i = 0; i < worker_queues.size() - 1; ++i) {
            const auto idx = (thread_idx + i + 1) % worker_queues.size();
            if (worker_queues[idx]->try_steal(task)) {
#ifdef INTP_DEBUG
                ++(worker_queue_ptr->stealing);
#endif
                return true;
            }
        }
        return false;
    }

    /**
     * @brief  scheduler function
     *
     */
    void thread_loop(size_t idx) {
        thread_idx = idx;
        worker_queue_ptr = worker_queues[idx].get();
        task_type task;
        // Fetch task from local queue, main queue and other thread's local
        // queue in order.
        while (!should_terminate) {
            if (worker_queue_ptr->try_pop(task) || try_pop_from_main(task) ||
                try_steal_from_others(task)) {
                task();
            } else {
                lock_type lk(main_queue_mutex);
                // Wait until there are tasks in main queue (awaked by
                // notify_one in queue_task), or the thread pool is being
                // shutdown (awaked by notify_all in destructor)
                cv.wait(lk, [this] {
                    return !main_queue.empty() || should_terminate;
                });
            }
        }
    }

#ifdef INTP_DEBUG
    size_t submitted{};
#endif
    bool should_terminate = false;  // Tells threads to stop looking for tasks
    std::mutex main_queue_mutex;    // Protects main task queue
    std::condition_variable cv;     // Signals for thread sleep/awake
    std::vector<std::thread> threads;  // Thread container

    std::queue<task_type> main_queue;  // Main queue for tasks
    inline static thread_local size_t thread_idx{};
    // Local queue ptr for tasks
    inline static thread_local shared_working_queue* worker_queue_ptr{};
    std::vector<std::unique_ptr<shared_working_queue>> worker_queues;

    JoinThreads join_threads;  // Defined last to ensure destruct first
};

}  // namespace intp

#endif  // INTP_THREAD_POOL
