// thread_pool.h — a barrier thread pool.
//
// Spawns n_threads − 1 persistent workers at construction. Each run_parallel(fn) runs fn(tid) on every worker
// (tid 1..n−1) while the caller runs fn(0), and returns once all have finished. An exception from any fn(tid) is
// rethrown by run_parallel, after every worker has finished.
#pragma once

#include <condition_variable>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace rigel {

// A thread budget as the configs name it: a positive count, or 0 for every core.
inline int thread_budget(int n_threads) {
    if (n_threads > 0) return n_threads;
    const int hw = static_cast<int>(std::thread::hardware_concurrency());
    return hw > 0 ? hw : 1;
}

class ThreadPool {
public:
    explicit ThreadPool(int n_threads) : n_workers_(n_threads - 1) {
        workers_.reserve(n_workers_);
        for (int t = 0; t < n_workers_; ++t) workers_.emplace_back(&ThreadPool::worker_loop, this, t + 1);
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lk(mtx_);
            shutdown_ = true;
        }
        cv_work_.notify_all();
        for (auto& w : workers_) w.join();
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    // Total thread count, the calling thread (worker 0) included.
    int n_threads() const { return n_workers_ + 1; }

    // Runs fn(tid) for tid 0..n_threads()-1, the caller taking tid 0, and blocks until every worker has finished. The
    // first exception any fn(tid) threw is rethrown here.
    template <typename Fn>
    void run_parallel(Fn&& fn) {
        {
            std::lock_guard<std::mutex> lk(mtx_);
            task_ = std::forward<Fn>(fn);
            active_workers_ = n_workers_;
            error_ = nullptr;
            ++generation_;
        }
        cv_work_.notify_all();
        run_task(0);
        std::unique_lock<std::mutex> lk(mtx_);
        cv_done_.wait(lk, [&] { return active_workers_ == 0; });
        if (error_) std::rethrow_exception(std::exchange(error_, nullptr));
    }

private:
    void run_task(int tid) noexcept {
        try {
            task_(tid);
        } catch (...) {
            std::lock_guard<std::mutex> lk(mtx_);
            if (!error_) error_ = std::current_exception();
        }
    }

    void worker_loop(int tid) {
        uint64_t my_gen = 0;
        std::unique_lock<std::mutex> lk(mtx_);
        while (true) {
            cv_work_.wait(lk, [&] { return generation_ > my_gen || shutdown_; });
            if (shutdown_) return;
            my_gen = generation_;
            lk.unlock();
            // run_parallel replaces task_ only after every worker has finished with it.
            run_task(tid);
            lk.lock();
            if (--active_workers_ == 0) cv_done_.notify_one();
        }
    }

    int n_workers_;
    std::vector<std::thread> workers_;

    std::mutex mtx_;
    std::condition_variable cv_work_;
    std::condition_variable cv_done_;

    std::function<void(int)> task_;
    std::exception_ptr error_;   // the first exception of the current run
    int active_workers_ = 0;
    uint64_t generation_ = 0;
    bool shutdown_ = false;
};

}  // namespace rigel
