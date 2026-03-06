/**
 * thread_queue.h — Threading infrastructure for parallel BAM scanning.
 *
 * Provides:
 *   - BoundedQueue<T>: bounded SPMC queue (single producer, multiple consumers)
 *   - QnameGroup: work unit containing deep-copied bam1_t* records
 *   - WorkerState: per-thread mutable state for fragment resolution
 */

#pragma once

#include <condition_variable>
#include <cstdint>
#include <deque>
#include <mutex>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "resolve_context.h"

namespace hulk {

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

// ================================================================
// QnameGroup — work unit for one query-name group
// ================================================================

struct QnameGroup {
    std::vector<bam1_t*> records;  // deep-copied via bam_dup1
    int64_t frag_id = 0;

    QnameGroup() = default;

    ~QnameGroup() {
        for (auto* b : records) {
            if (b) bam_destroy1(b);
        }
    }

    // Move-only (owns bam1_t pointers)
    QnameGroup(const QnameGroup&) = delete;
    QnameGroup& operator=(const QnameGroup&) = delete;

    QnameGroup(QnameGroup&& o) noexcept
        : records(std::move(o.records)), frag_id(o.frag_id)
    {
        o.records.clear();
    }

    QnameGroup& operator=(QnameGroup&& o) noexcept {
        if (this != &o) {
            for (auto* b : records)
                if (b) bam_destroy1(b);
            records = std::move(o.records);
            frag_id = o.frag_id;
            o.records.clear();
        }
        return *this;
    }
};

// ================================================================
// WorkerState — defined in bam_scanner.cpp where stat structs
// are available. This header only provides the queue + QnameGroup.
// ================================================================

}  // namespace hulk
