#pragma once

#include <stdexcept>

namespace eonc {

class DimerModeLostException : public std::runtime_error {
public:
  DimerModeLostException()
      : std::runtime_error("Dimer lost mode alignment") {}
};

class DimerModeRestoredException : public std::exception {
public:
  const char *what() const noexcept override {
    return "Dimer lost mode but restored to best negative curvature state";
  }
};

} // namespace eonc
