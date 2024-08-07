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
#include "client/matter/PeriodicBoundaryConditions.hpp"

namespace eonc::mat {

AtomMatrix PBC_A1::operator()(const AtomMatrix &diff) const {
  AtomMatrix ddiff = diff;
  double box =
      cell(0, 0); // Assuming a cubic box, all diagonal elements are equal
  double box2_r = 1.0 / (0.5 * box); // box2_r = 2.0 / box

  for (int i = 0; i < ddiff.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      int k = static_cast<int>(ddiff(i, j) * box2_r);
      ddiff(i, j) -= k * box;
      // Adjust to make sure values are within [-box/2, box/2]
      if (ddiff(i, j) >= box / 2.0)
        ddiff(i, j) -= box;
      if (ddiff(i, j) < -box / 2.0)
        ddiff(i, j) += box;
    }
  }

  return ddiff;
}

AtomMatrix PBC_B1::operator()(const AtomMatrix &diff) const {
  AtomMatrix ddiff = diff * cellInverse;

  for (int i = 0; i < ddiff.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      int k = static_cast<int>(ddiff(i, j) * 0.5);
      ddiff(i, j) -= k;
      k = static_cast<int>(ddiff(i, j) * 0.5);
      ddiff(i, j) -= k;
    }
  }

  return ddiff * cell;
}

AtomMatrix PBC_B2::operator()(const AtomMatrix &diff) const {
  AtomMatrix ddiff = diff * cellInverse;

  for (int i = 0; i < ddiff.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      if (ddiff(i, j) > 0.5) // box2 = 0.5 for normalized box
        ddiff(i, j) -= 1.0;
      else if (ddiff(i, j) < -0.5)
        ddiff(i, j) += 1.0;
    }
  }

  return ddiff * cell;
}

AtomMatrix PBC_B3::operator()(const AtomMatrix &diff) const {
  AtomMatrix ddiff = diff * cellInverse;

  for (int i = 0; i < ddiff.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      ddiff(i, j) = std::fabs(ddiff(i, j));
      if (ddiff(i, j) > 0.5) // box2 = 0.5 for normalized box
        ddiff(i, j) -= 1.0;
    }
  }

  return ddiff * cell;
}

AtomMatrix PBC_B4::operator()(const AtomMatrix &diff) const {
  AtomMatrix ddiff = diff * cellInverse;

  for (int i = 0; i < ddiff.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      ddiff(i, j) = ddiff(i, j) * 1.0 + 3.0;
      int k = static_cast<int>(ddiff(i, j) + 0.5);
      ddiff(i, j) = (ddiff(i, j) - k) * 1.0;
    }
  }

  return ddiff * cell;
}

AtomMatrix PBC_C1::operator()(const AtomMatrix &diff) const {
  AtomMatrix ddiff = diff * cellInverse;

  for (int i = 0; i < ddiff.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      ddiff(i, j) = std::remainder(ddiff(i, j), 1.0);
    }
  }

  return ddiff * cell;
}

AtomMatrix PBC_C5::operator()(const AtomMatrix &diff) const {
  AtomMatrix ddiff = diff * cellInverse;

  for (int i = 0; i < ddiff.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      int k = static_cast<int>(ddiff(i, j) * 1.0 +
                               ((ddiff(i, j) >= 0.0) ? 0.5 : -0.5));
      ddiff(i, j) -= k;
    }
  }

  return ddiff * cell;
}

AtomMatrix PBC_C6::operator()(const AtomMatrix &diff) const {
  AtomMatrix ddiff = diff * cellInverse;

  for (int i = 0; i < ddiff.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      int k = static_cast<int>(std::fabs(ddiff(i, j) * 1.0) + 0.5);
      ddiff(i, j) -= k;
    }
  }

  return ddiff * cell;
}

AtomMatrix PBC_ESVN::operator()(const AtomMatrix &diff) const {
  AtomMatrix ddiff = diff * cellInverse;

  for (int i = 0; i < ddiff.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      ddiff(i, j) = fmod(fmod(ddiff(i, j), 1.0) + 1.5, 1.0) - .5;
    }
  }

  return ddiff * cell;
}

} // namespace eonc::mat
