/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once

#include <atomic>
#include <memory>
#include <stdexcept>
#include <string>

namespace eonc {

/// Cooperative stop. Loops poll and throw; they do not abort the process.
class JobCancelled : public std::runtime_error {
public:
  explicit JobCancelled(const char *site = "cancelled")
      : std::runtime_error(std::string("job cancelled at ") +
                           (site != nullptr ? site : "cancelled")) {}
};

/// Shared flag. Copies observe the same request.
class CancelToken {
public:
  CancelToken() : flag_(std::make_shared<std::atomic<bool>>(false)) {}

  void request() noexcept {
    if (flag_) {
      flag_->store(true, std::memory_order_release);
    }
  }

  void reset() noexcept {
    if (flag_) {
      flag_->store(false, std::memory_order_release);
    }
  }

  [[nodiscard]] bool requested() const noexcept {
    return flag_ && flag_->load(std::memory_order_acquire);
  }

  void poll(const char *site = "step") const {
    if (requested()) {
      throw JobCancelled(site);
    }
  }

private:
  std::shared_ptr<std::atomic<bool>> flag_;
};

} // namespace eonc
