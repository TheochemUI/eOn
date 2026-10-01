#pragma once
// Bounded thread pool for per-image work (NEB image forces).

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace eonc {

// Calls work(1) .. work(n) on at most hardware_concurrency() threads. Each
// thread takes the next index until none are left, so a band with more images
// than cores does not oversubscribe. The first exception thrown by any call is
// rethrown here after every thread has joined; the remaining indices are still
// evaluated. std::thread rather than std::jthread: Apple Clang's libc++ does
// not ship the latter.
template <typename Work> inline void forEachImage(long n, Work &&work) {
  if (n <= 0)
    return;
  const long hw =
      std::max(1L, static_cast<long>(std::thread::hardware_concurrency()));
  const long nThreads = std::min(n, hw);
  std::atomic<long> next{1};
  std::exception_ptr failure;
  std::mutex failureMutex;
  auto worker = [&] {
    for (long i = next.fetch_add(1); i <= n; i = next.fetch_add(1)) {
      try {
        work(i);
      } catch (...) {
        std::lock_guard<std::mutex> lock(failureMutex);
        if (!failure)
          failure = std::current_exception();
      }
    }
  };
  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(nThreads - 1));
  try {
    for (long t = 1; t < nThreads; t++)
      threads.emplace_back(worker);
  } catch (...) {
    // Thread creation failed; the threads already running finish the band.
    std::lock_guard<std::mutex> lock(failureMutex);
    if (!failure)
      failure = std::current_exception();
  }
  worker();
  for (auto &t : threads)
    t.join();
  if (failure)
    std::rethrow_exception(failure);
}

} // namespace eonc
