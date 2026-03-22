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

namespace eonc {

// Random number generator constants
constexpr double IM = 2147483647.0;
constexpr double AM = 1.0 / IM;
constexpr int NTAB = 32;
constexpr int NDIV = 1 + (IM / NTAB);
constexpr double EPS = 1.2e-7;
constexpr double RNMX = 1.0 - EPS;
constexpr long IM1 = 2147483563;
constexpr long IM2 = 2147483399;
constexpr long IMM1 = IM1 - 1;
constexpr long IA1 = 40014;
constexpr long IA2 = 40692;
constexpr long IQ1 = 53668;
constexpr long IQ2 = 52774;
constexpr long IR1 = 12211;
constexpr long IR2 = 3791;

namespace rng {

double random(long newSeed = 0);
double randomDouble();
double randomDouble(int max);
double randomDouble(long max);
double randomDouble(double max);
long randomInt(int lower, int upper);
double gaussRandom(double avg, double std);

} // namespace rng

} // namespace eonc
