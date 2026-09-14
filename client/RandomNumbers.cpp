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
#include "eon/RandomNumbers.h"

#include <cmath>

namespace {
struct Ran2State {
  long seed{-1};
  long seed2{123456789};
  long iy{0};
  long iv[eonc::NTAB]{};
};

// Parallel replica exchange (and any other std::thread MD) used to race on
// the process-wide ran2 tables. Each C++ thread owns a stream.
thread_local Ran2State tls;
} // namespace

double eonc::rng::random(long newSeed) {
  auto &st = tls;
  if (newSeed) {
    st.seed = -newSeed;
  }
  int j;
  long k;
  double temp;
  if (st.seed <= 0) {
    if (-(st.seed) < 1)
      st.seed = 3;
    else
      st.seed = -(st.seed);
    st.seed2 = (st.seed);
    for (j = NTAB + 7; j >= 0; j--) {
      k = (st.seed) / IQ1;
      st.seed = IA1 * (st.seed - k * IQ1) - k * IR1;
      if (st.seed < 0)
        st.seed += IM1;
      if (j < NTAB)
        st.iv[j] = st.seed;
    }
    st.iy = st.iv[0];
  }
  k = (st.seed) / IQ1;
  st.seed = IA1 * (st.seed - k * IQ1) - k * IR1;
  if (st.seed < 0)
    st.seed += IM1;
  k = st.seed2 / IQ2;
  st.seed2 = IA2 * (st.seed2 - k * IQ2) - k * IR2;
  if (st.seed2 < 0)
    st.seed2 += IM2;
  j = int(st.iy / NDIV);
  st.iy = st.iv[j] - st.seed2;
  st.iv[j] = st.seed;
  if (st.iy < 1)
    st.iy += IMM1;
  if ((temp = double(AM * st.iy)) > RNMX)
    return RNMX;
  else
    return temp;
}

double eonc::rng::randomDouble() { return (random()); }

double eonc::rng::randomDouble(int max) {
  double dmax = double(max);
  return (dmax * randomDouble());
}

double eonc::rng::randomDouble(long max) {
  double dmax = double(max);
  return (dmax * randomDouble());
}

double eonc::rng::randomDouble(double dmax) { return (dmax * randomDouble()); }

long eonc::rng::randomInt(int lower, int upper) {
  return lround((upper - lower) * randomDouble() + lower);
}

double eonc::rng::gaussRandom(double avg, double std) {
  double r = 2, v1, v2, l, result;
  while (r >= 1.0 || r < 1e-300) {
    v1 = 2.0 * randomDouble() - 1.0;
    v2 = 2.0 * randomDouble() - 1.0;
    r = v1 * v1 + v2 * v2;
  }
  l = v1 * sqrt(-2.0 * ::log(r) / r);
  result = avg + std * l;
  return (result);
}
