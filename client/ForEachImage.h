#pragma once
// Persistent bounded thread pool for per-image work (NEB image forces and
// projections).

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#endif

namespace eonc {

namespace detail {

// Cores this process may run on: the affinity mask under Slurm or taskset,
// else hardware_concurrency(). A band on 4 allocated cores of a 32-core
// node then gets 4 threads, not 32 that share 4 cores.
inline long usableCores() {
#if defined(__linux__)
  cpu_set_t set;
  CPU_ZERO(&set);
  if (sched_getaffinity(0, sizeof(set), &set) == 0) {
    const int n = CPU_COUNT(&set);
    if (n > 0)
      return n;
  }
#endif
  return std::max(1L, static_cast<long>(std::thread::hardware_concurrency()));
}

// Workers start once and wait for the next band, so a band costs a
// condition-variable round trip instead of creating and joining threads.
// One band runs at a time. A call made while a band runs (from inside a
// work item, or from a second thread) runs serially on its caller, so the
// pool never waits on itself.
class ImagePool {
public:
  static ImagePool &instance() {
    static ImagePool pool;
    return pool;
  }

  ImagePool(const ImagePool &) = delete;
  ImagePool &operator=(const ImagePool &) = delete;

  ~ImagePool() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
    }
    wake_.notify_all();
    for (auto &t : threads_)
      t.join();
  }

  // Calls work(1) .. work(n), the caller taking indices too, and rethrows
  // the first exception after every index has run.
  void run(long n, const std::function<void(long)> &work) {
    if (n <= 0)
      return;
    std::unique_lock<std::mutex> busy(runMutex_, std::try_to_lock);
    if (!busy.owns_lock() || threads_.empty() || n == 1) {
      serial(n, work);
      return;
    }
    {
      std::lock_guard<std::mutex> lock(mutex_);
      work_ = &work;
      n_ = n;
      next_.store(1);
      failure_ = nullptr;
      // Only as many helpers as there are images beyond the caller's.
      helpers_ = std::min<long>(static_cast<long>(threads_.size()), n - 1);
      active_ = helpers_;
      ++generation_;
    }
    wake_.notify_all();
    drain();
    std::unique_lock<std::mutex> lock(mutex_);
    done_.wait(lock, [this] { return active_ == 0; });
    work_ = nullptr;
    if (failure_)
      std::rethrow_exception(failure_);
  }

  [[nodiscard]] std::size_t threads() const noexcept { return threads_.size(); }

private:
  ImagePool() {
    const long helpers = usableCores() - 1;
    for (long t = 0; t < helpers; t++) {
      try {
        threads_.emplace_back([this, t] { loop(t); });
      } catch (...) {
        break; // the threads already running serve every band
      }
    }
  }

  static void serial(long n, const std::function<void(long)> &work) {
    std::exception_ptr failure;
    for (long i = 1; i <= n; i++) {
      try {
        work(i);
      } catch (...) {
        if (!failure)
          failure = std::current_exception();
      }
    }
    if (failure)
      std::rethrow_exception(failure);
  }

  void drain() {
    for (long i = next_.fetch_add(1); i <= n_; i = next_.fetch_add(1)) {
      try {
        (*work_)(i);
      } catch (...) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!failure_)
          failure_ = std::current_exception();
      }
    }
  }

  void loop(long index) {
    std::size_t seen = 0;
    for (;;) {
      {
        std::unique_lock<std::mutex> lock(mutex_);
        wake_.wait(lock, [&] {
          return stop_ || (generation_ != seen && index < helpers_);
        });
        if (stop_)
          return;
        seen = generation_;
      }
      drain();
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (--active_ == 0)
          done_.notify_one();
      }
    }
  }

  std::vector<std::thread> threads_;
  std::mutex runMutex_;
  std::mutex mutex_;
  std::condition_variable wake_;
  std::condition_variable done_;
  const std::function<void(long)> *work_{nullptr};
  long n_{0};
  std::atomic<long> next_{1};
  std::exception_ptr failure_;
  long helpers_{0};
  long active_{0};
  std::size_t generation_{0};
  bool stop_{false};
};

} // namespace detail

// Calls work(1) .. work(n) on at most usableCores() threads of a pool that
// lives for the process. Each thread takes the next index until none are
// left, so a band with more images than cores does not oversubscribe. The
// first exception thrown by any call is rethrown here after every index has
// run.
template <typename Work> inline void forEachImage(long n, Work &&work) {
  const std::function<void(long)> f = std::ref(work);
  detail::ImagePool::instance().run(n, f);
}

} // namespace eonc
