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
#include "eon/ExactReduction.h"

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>

#if defined(_MSC_VER) && !defined(__SIZEOF_INT128__) && defined(_M_X64)
#include <intrin.h>
#endif

namespace eonc::exact {
namespace {

// A finite double is M * 2^e with M < 2^53 and e >= -1074, so a product of
// two is P * 2^E with P < 2^106 and E in [-2148, 1942]. Bit k of the
// accumulator weighs 2^(k - kBase). Chunks of 32 bits sit in signed 64-bit
// words: every addend piece is below 2^32, so a word absorbs 2^31 pieces
// before it can overflow, and signs need no borrow chains.
constexpr int kBase = 2148;
constexpr int kTopBit = 1942 + 106 + kBase; // one past the highest product bit
constexpr int kChunks = kTopBit / 32 + 4;   // headroom for carries and sign
constexpr std::size_t kNormalizeEvery = std::size_t{1} << 30;
constexpr uint64_t kMask32 = 0xffffffffULL;

struct Decoded {
  uint64_t m; // significand as an integer, < 2^53
  int e;      // value = m * 2^e
  bool neg;
};

inline bool decode(double v, Decoded &d) {
  uint64_t b;
  std::memcpy(&b, &v, sizeof b);
  const auto field = static_cast<int>((b >> 52) & 0x7ff);
  if (field == 0x7ff) {
    return false;
  }
  const uint64_t frac = b & ((uint64_t{1} << 52) - 1);
  d.neg = (b >> 63) != 0;
  if (field == 0) {
    d.m = frac;
    d.e = -1074;
  } else {
    d.m = frac | (uint64_t{1} << 52);
    d.e = field - 1075;
  }
  return true;
}

// Full 128-bit product of two significands, as (hi, lo).
inline void mul64(uint64_t a, uint64_t b, uint64_t &hi, uint64_t &lo) {
#if defined(__SIZEOF_INT128__)
  const unsigned __int128 p = static_cast<unsigned __int128>(a) * b;
  lo = static_cast<uint64_t>(p);
  hi = static_cast<uint64_t>(p >> 64);
#elif defined(_MSC_VER) && defined(_M_X64)
  lo = _umul128(a, b, &hi);
#else
  const uint64_t a0 = a & kMask32, a1 = a >> 32;
  const uint64_t b0 = b & kMask32, b1 = b >> 32;
  const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
  const uint64_t mid = (p00 >> 32) + (p01 & kMask32) + (p10 & kMask32);
  lo = (p00 & kMask32) | (mid << 32);
  hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
#endif
}

class Accumulator {
public:
  // Words outside [m_lo, m_hi) are never read before they are cleared, so
  // a short vector touches a few words, not all ~130.
  void add(uint64_t hi, uint64_t lo, int e, bool neg) {
    if ((hi | lo) == 0) {
      return;
    }
    const int pos = e + kBase;
    const int c = pos >> 5;
    const int s = pos & 31;
    // (hi, lo) * 2^s is below 2^138: five 32-bit pieces from word c up.
    const uint64_t t0 = lo << s;
    const uint64_t t1 = (hi << s) | (s != 0 ? lo >> (64 - s) : 0);
    const uint64_t t2 = s != 0 ? hi >> (64 - s) : 0;
    const std::array<int64_t, 5> piece{
        static_cast<int64_t>(t0 & kMask32), static_cast<int64_t>(t0 >> 32),
        static_cast<int64_t>(t1 & kMask32), static_cast<int64_t>(t1 >> 32),
        static_cast<int64_t>(t2 & kMask32)};
    cover(c, c + 5);
    if (neg) {
      for (int k = 0; k < 5; ++k) {
        m_w[c + k] -= piece[k];
      }
    } else {
      for (int k = 0; k < 5; ++k) {
        m_w[c + k] += piece[k];
      }
    }
  }

  // Carry every word into [0, 2^32); the carry out of the top word is the
  // sign (0 or -1).
  int64_t normalize() {
    if (m_lo >= m_hi) {
      return 0;
    }
    cover(m_lo, kChunks);
    int64_t carry = 0;
    for (int k = m_lo; k < kChunks; ++k) {
      const int64_t v = m_w[k] + carry;
      m_w[k] = v & static_cast<int64_t>(kMask32);
      carry = v >> 32; // arithmetic shift: floor division
    }
    m_hi = kChunks;
    return carry;
  }

  double round() {
    const int64_t sign = normalize();
    if (m_lo >= m_hi) {
      return 0.0;
    }
    const bool neg = sign < 0;
    if (neg) {
      // Two's complement over the covered words; the words below m_lo are
      // zero and stay zero.
      int64_t carry = 1;
      for (int k = m_lo; k < kChunks; ++k) {
        const int64_t v = ((~m_w[k]) & static_cast<int64_t>(kMask32)) + carry;
        m_w[k] = v & static_cast<int64_t>(kMask32);
        carry = v >> 32;
      }
    }
    int h = kChunks - 1;
    while (h >= m_lo && m_w[h] == 0) {
      --h;
    }
    if (h < m_lo) {
      return 0.0;
    }
    const int tb =
        32 * h + 31 - std::countl_zero(static_cast<uint32_t>(m_w[h]));
    // Lowest kept bit: 53 significant bits, but never below 2^-1074.
    int lsb = tb - 52;
    if (lsb < kBase - 1074) {
      lsb = kBase - 1074;
    }
    uint64_t m = bits(lsb, tb - lsb + 1);
    const bool half = bit(lsb - 1);
    const bool sticky = below(lsb - 1);
    if (half && (sticky || (m & 1U) != 0)) {
      ++m;
      if (m == (uint64_t{1} << 53)) {
        m >>= 1;
        ++lsb;
      }
    }
    // m < 2^53 times a power of two in range is exact; past 2^1024 it is
    // the correctly rounded infinity.
    const double r = std::ldexp(static_cast<double>(m), lsb - kBase);
    return neg ? -r : r;
  }

  void maybeNormalize(std::size_t count) {
    if (count % kNormalizeEvery == 0) {
      // Fold the carries back; the sign word stays in range because the
      // true sum is far below 2^(32 * kChunks - 1).
      const int64_t sign = normalize();
      m_w[kChunks - 1] += sign * (int64_t{1} << 32);
    }
  }

private:
  void cover(int lo, int hi) {
    if (m_lo >= m_hi) {
      for (int k = lo; k < hi; ++k) {
        m_w[k] = 0;
      }
      m_lo = lo;
      m_hi = hi;
      return;
    }
    for (int k = lo; k < m_lo; ++k) {
      m_w[k] = 0;
    }
    for (int k = m_hi; k < hi; ++k) {
      m_w[k] = 0;
    }
    m_lo = lo < m_lo ? lo : m_lo;
    m_hi = hi > m_hi ? hi : m_hi;
  }

  bool bit(int k) const {
    if (k < 0 || (k >> 5) < m_lo) {
      return false;
    }
    return ((m_w[k >> 5] >> (k & 31)) & 1) != 0;
  }

  // Any set bit strictly below bit k.
  bool below(int k) const {
    if (k <= 0) {
      return false;
    }
    const int c = k >> 5;
    for (int j = m_lo; j < c; ++j) {
      if (m_w[j] != 0) {
        return true;
      }
    }
    if (c < m_lo) {
      return false;
    }
    const uint64_t partMask = (uint64_t{1} << (k & 31)) - 1;
    return (static_cast<uint64_t>(m_w[c]) & partMask) != 0;
  }

  // Bits [k, k + width) for width <= 53; k >= 0.
  uint64_t bits(int k, int width) const {
    const int c = k >> 5;
    const int off = k & 31;
    const auto word = [&](int w) -> uint64_t {
      return (w >= m_lo && w < kChunks) ? static_cast<uint64_t>(m_w[w]) : 0;
    };
    const uint64_t low = word(c) | (word(c + 1) << 32);
    const uint64_t high = word(c + 2);
    const uint64_t v = (low >> off) | (off != 0 ? high << (64 - off) : 0);
    return v & ((uint64_t{1} << width) - 1);
  }

  std::array<int64_t, kChunks> m_w; // left uninitialised; see cover()
  int m_lo = 0;
  int m_hi = 0;
};

} // namespace

double dot(const double *x, const double *y, std::size_t n) {
  Accumulator acc;
  for (std::size_t i = 0; i < n; ++i) {
    Decoded a;
    Decoded b;
    if (!decode(x[i], a) || !decode(y[i], b)) [[unlikely]] {
      double plain = 0.0;
      for (std::size_t j = 0; j < n; ++j) {
        plain += x[j] * y[j];
      }
      return plain;
    }
    uint64_t hi;
    uint64_t lo;
    mul64(a.m, b.m, hi, lo);
    acc.add(hi, lo, a.e + b.e, a.neg != b.neg);
    acc.maybeNormalize(i + 1);
  }
  return acc.round();
}

double sum(const double *x, std::size_t n) {
  Accumulator acc;
  for (std::size_t i = 0; i < n; ++i) {
    Decoded a;
    if (!decode(x[i], a)) [[unlikely]] {
      double plain = 0.0;
      for (std::size_t j = 0; j < n; ++j) {
        plain += x[j];
      }
      return plain;
    }
    acc.add(0, a.m, a.e, a.neg);
    acc.maybeNormalize(i + 1);
  }
  return acc.round();
}

} // namespace eonc::exact
