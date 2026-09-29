#ifndef INTP_TEST_TIMER
#define INTP_TEST_TIMER

#include <chrono>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

class Timer {
   public:
    static Timer& get_timer();
    void start(std::string func_name);
    void pause_and_start(std::string func_name);
    void pause();
    void pause(std::string func_name);
    void reset();
    void print() const;

    std::chrono::high_resolution_clock::duration get_duration(
        std::string) const;

   private:
    Timer() = default;
    Timer(const Timer&) = delete;
    Timer(Timer&&) = delete;
    Timer& operator=(const Timer&) = delete;
    Timer& operator=(Timer&&) = delete;

    std::vector<std::string> entries;
    std::unordered_map<
        std::string,
        std::pair<std::chrono::high_resolution_clock::duration,
                  std::chrono::high_resolution_clock::time_point>>
        time_consuming;

    std::string current_func_name;
};

#endif  // INTP_TEST_TIMER
