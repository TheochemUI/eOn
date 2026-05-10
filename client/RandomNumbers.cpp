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
#include "RandomNumbers.h"

#include <cmath>

double eonc::rng::random(long newSeed) {
  static long seed = -1;
  if (newSeed) {
    seed = -newSeed;
  }
  int j;
  long k;
  static long seed2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;
  if (seed <= 0) {
    if (-(seed) < 1)
      seed = 3;
    else
      seed = -(seed);
    seed2 = (seed);
    for (j = NTAB + 7; j >= 0; j--) {
      k = (seed) / IQ1;
      seed = IA1 * (seed - k * IQ1) - k * IR1;
      if (seed < 0)
        seed += IM1;
      if (j < NTAB)
        iv[j] = seed;
    }
    iy = iv[0];
  }
  k = (seed) / IQ1;
  seed = IA1 * (seed - k * IQ1) - k * IR1;
  if (seed < 0)
    seed += IM1;
  k = seed2 / IQ2;
  seed2 = IA2 * (seed2 - k * IQ2) - k * IR2;
  if (seed2 < 0)
    seed2 += IM2;
  j = int(iy / NDIV);
  iy = iv[j] - seed2;
  iv[j] = seed;
  if (iy < 1)
    iy += IMM1;
  // Compute outside the if-condition so the side-effecting
  // assignment-in-if pattern doesn't fire bugprone-assignment-in-if-
  // condition. Behaviour is identical -- temp gets the AM*iy value
  // either way.
  temp = static_cast<double>(AM * iy);
  if (temp > RNMX)
    return RNMX;
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
