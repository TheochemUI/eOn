#include "../ForEachImage.h"
#include "catch2/catch_amalgamated.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

TEST_CASE("forEachImage calls every image index once", "[neb][threads]") {
  const long n = 37;
  std::vector<std::atomic<int>> calls(n + 1);
  eonc::forEachImage(n, [&](long i) { calls[static_cast<size_t>(i)]++; });
  REQUIRE(calls[0] == 0);
  for (long i = 1; i <= n; i++)
    REQUIRE(calls[static_cast<size_t>(i)] == 1);
}

TEST_CASE("forEachImage with no images calls nothing", "[neb][threads]") {
  int calls = 0;
  eonc::forEachImage(0, [&](long) { calls++; });
  REQUIRE(calls == 0);
}

TEST_CASE("forEachImage runs at most hardware_concurrency images at once",
          "[neb][threads]") {
  const long hw =
      std::max(1L, static_cast<long>(std::thread::hardware_concurrency()));
  const long n = 4 * hw + 3;
  std::atomic<long> active{0};
  std::atomic<long> peak{0};
  eonc::forEachImage(n, [&](long) {
    long now = ++active;
    long seen = peak.load();
    while (now > seen && !peak.compare_exchange_weak(seen, now)) {
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(2));
    --active;
  });
  REQUIRE(peak.load() >= 1);
  REQUIRE(peak.load() <= hw);
}

TEST_CASE("forEachImage rethrows an image's exception after the others finish",
          "[neb][threads]") {
  const long n = 20;
  std::vector<std::atomic<int>> calls(n + 1);
  auto work = [&](long i) {
    calls[static_cast<size_t>(i)]++;
    if (i == 7)
      throw std::runtime_error("image 7 failed");
  };
  REQUIRE_THROWS_WITH(eonc::forEachImage(n, work), "image 7 failed");
  for (long i = 1; i <= n; i++)
    REQUIRE(calls[static_cast<size_t>(i)] == 1);
}

TEST_CASE("forEachImage keeps its threads across bands", "[neb][threads]") {
  // Thread ids seen over many bands stay within the pool plus the caller:
  // no band creates threads of its own.
  const auto pool = eonc::detail::ImagePool::instance().threads();
  std::mutex m;
  std::vector<std::thread::id> ids;
  for (int band = 0; band < 50; band++) {
    eonc::forEachImage(8, [&](long) {
      std::lock_guard<std::mutex> lock(m);
      if (std::find(ids.begin(), ids.end(), std::this_thread::get_id()) ==
          ids.end())
        ids.push_back(std::this_thread::get_id());
      std::this_thread::sleep_for(std::chrono::microseconds(50));
    });
  }
  REQUIRE(ids.size() <= pool + 1);
}

TEST_CASE("forEachImage inside an image runs serially", "[neb][threads]") {
  // A nested band must not wait on the pool that runs its caller.
  std::atomic<long> inner{0};
  eonc::forEachImage(
      6, [&](long) { eonc::forEachImage(5, [&](long) { inner++; }); });
  REQUIRE(inner.load() == 30);
}

TEST_CASE("forEachImage uses at most the cores this process may run on",
          "[neb][threads]") {
  REQUIRE(eonc::detail::ImagePool::instance().threads() + 1 <=
          static_cast<size_t>(eonc::detail::usableCores()));
}
