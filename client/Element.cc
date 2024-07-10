#include "Element.hpp"
#include "magic_enum/magic_enum.hpp"

namespace eonc {

std::string mass2atom(double atomicmass) {
  for (const auto &[element, properties] : elementData) {
    if (std::abs(properties.atomicMass - atomicmass) < 0.5) {
      return std::string(magic_enum::enum_name(element));
    }
  }
  return "Unknown";
}

size_t symbol2atomicNumber(const std::string_view &symbol) {
  auto element = magic_enum::enum_cast<Element>(symbol);
  if (element.has_value()) {
    return static_cast<int>(element.value());
  }
  // invalid symbol
  return -1;
}

std::string atomicNumber2symbol(size_t n) {
  if (n >= 0 && n < static_cast<int>(Element::COUNTER)) {
    return std::string(magic_enum::enum_name(static_cast<Element>(n)));
  }
  return "Unknown";
}

} // namespace eonc
