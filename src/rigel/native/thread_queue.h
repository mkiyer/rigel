/**
 * thread_queue.h — Bounded SPMC work queue for parallel BAM scanning.
 *
 * Provides BoundedQueue<T>: a single-producer multiple-consumer queue
 * with backpressure.  The work-unit type (QnameGroup) is defined in
 * bam_scanner.cpp alongside the other BAM-specific structures.
 */

#pragma once

#include <condition_variable>
#include <cstdint>
#include <deque>
#include <mutex>

namespace rigel {

// ================================================================
// BoundedQueue — bounded SPMC queue with backpressure
// ================================================================

template <typename T>
class BoundedQueue {
    std::deque<T> items_;
    std::mutex mutex_;
    std::condition_variable not_empty_;
    std::condition_variable not_full_;
    size_t capacity_;
    bool closed_ = false;

public:
    explicit BoundedQueue(size_t capacity) : capacity_(capacity) {}

    /// Producer: push an item, blocking if queue is full.
    void push(T item) {
        std::unique_lock<std::mutex> lock(mutex_);
        not_full_.wait(lock, [this] {
            return items_.size() < capacity_;
        });
        items_.push_back(std::move(item));
        not_empty_.notify_one();
    }

    /// Producer: signal no more items will be pushed.
    void close() {
        std::unique_lock<std::mutex> lock(mutex_);
        closed_ = true;
        not_empty_.notify_all();
    }

    /// Consumer: pop an item. Returns false when queue is closed AND empty.
    bool pop(T& item) {
        std::unique_lock<std::mutex> lock(mutex_);
        not_empty_.wait(lock, [this] {
            return !items_.empty() || closed_;
        });
        if (items_.empty()) return false;  // closed and drained
        item = std::move(items_.front());
        items_.pop_front();
        not_full_.notify_one();
        return true;
    }
};

}  // namespace rigel
