/**
 * thread_pool.h — Lightweight barrier-based thread pool for the EM E-step.
 *
 * Spawns n_threads-1 persistent worker threads at construction.  Each call
 * to run_parallel(fn) fans out fn(tid) across all workers (tid 1..n-1)
 * while the caller executes fn(0).  Blocks until all workers complete.
 *
 * This eliminates per-E-step pthread_create/join overhead (~20-50 µs each
 * on macOS) which accumulates to seconds over thousands of SQUAREM
 * iterations on mega-loci.
 */

#ifndef RIGEL_THREAD_POOL_H
#define RIGEL_THREAD_POOL_H

#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace rigel {

class EStepThreadPool {
public:
    explicit EStepThreadPool(int n_threads)
        : n_workers_(n_threads - 1)
    {
        workers_.reserve(n_workers_);
        for (int t = 0; t < n_workers_; ++t) {
            workers_.emplace_back(&EStepThreadPool::worker_loop, this, t + 1);
        }
    }

    ~EStepThreadPool() {
        {
            std::lock_guard<std::mutex> lk(mtx_);
            shutdown_ = true;
        }
        cv_work_.notify_all();
        for (auto& w : workers_) w.join();
    }

    // Non-copyable, non-movable
    EStepThreadPool(const EStepThreadPool&) = delete;
    EStepThreadPool& operator=(const EStepThreadPool&) = delete;

    /// Total thread count (including the calling thread as worker 0).
    int n_threads() const { return n_workers_ + 1; }

    /// Execute fn(tid) for tid 0..n_threads()-1.  The caller runs tid=0.
    /// Blocks until all workers have finished.
    template<typename Fn>
    void run_parallel(Fn&& fn) {
        {
            std::lock_guard<std::mutex> lk(mtx_);
            task_ = std::forward<Fn>(fn);
            active_workers_ = n_workers_;
            ++generation_;
        }
        cv_work_.notify_all();

        // Caller participates as worker 0
        task_(0);

        // Wait for workers 1..n to finish
        std::unique_lock<std::mutex> lk(mtx_);
        cv_done_.wait(lk, [&] { return active_workers_ == 0; });
    }

private:
    void worker_loop(int tid) {
        uint64_t my_gen = 0;
        std::unique_lock<std::mutex> lk(mtx_);
        while (true) {
            cv_work_.wait(lk, [&] {
                return generation_ > my_gen || shutdown_;
            });
            if (shutdown_) return;
            my_gen = generation_;
            auto fn = task_;
            lk.unlock();

            fn(tid);

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
    int active_workers_ = 0;
    uint64_t generation_ = 0;
    bool shutdown_ = false;
};

} // namespace rigel

#endif // RIGEL_THREAD_POOL_H
