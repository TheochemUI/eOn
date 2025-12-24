#pragma once

#include <stdexcept>

class DimerModeLostException : public std::runtime_error {
public:
  DimerModeLostException()
      : std::runtime_error("Dimer lost mode alignment") {}
};
