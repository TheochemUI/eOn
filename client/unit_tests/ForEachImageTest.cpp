#include "../ForEachImage.h"
#include "catch2/catch_amalgamated.hpp"

#include <atomic>
#include <chrono>
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
