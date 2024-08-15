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
#include "MatterHelpers.hpp"
#include "client/BaseStructures.h"
#include "client/HelperFunctions.h"
#include "client/Parser.hpp"
#include "client/matter/MatterCreator.hpp"
#include <set>

// TODO(rg):: Cleanup, this shouldn't be required
namespace tmp {
struct atom {
  double r;
  int z;
};

struct by_atom {
  bool operator()(atom const &a, atom const &b) const {
    if (a.z != b.z) {
      return a.z < b.z;
    } else {
      return a.r < b.r;
    }
  }
};
} // namespace tmp

namespace eonc::mat {

bool rotationMatch(const Matter &m1, const Matter &m2, const double max_diff) {
  AtomMatrix r1 = m1.getPositions();
  AtomMatrix r2 = m2.getPositions();

  // Align centroids
  VectorType c1(3);
  VectorType c2(3);

  c1[0] = r1.col(0).sum();
  c1[1] = r1.col(1).sum();
  c1[2] = r1.col(2).sum();
  c2[0] = r2.col(0).sum();
  c2[1] = r2.col(1).sum();
  c2[2] = r2.col(2).sum();
  c1 /= r1.rows();
  c2 /= r2.rows();

  for (int i = 0; i < r1.rows(); i++) {
    r1(i, 0) -= c1[0];
    r1(i, 1) -= c1[1];
    r1(i, 2) -= c1[2];

    r2(i, 0) -= c2[0];
    r2(i, 1) -= c2[1];
    r2(i, 2) -= c2[2];
  }

  RotationMatrix R = eonc::helper_functions::rotationExtract(r1, r2);

  // Eigen is transposed relative to numpy
  r2 = r2 * R;

  for (int i = 0; i < r1.rows(); i++) {
    double diff = (r2.row(i) - r1.row(i)).norm();
    if (diff > max_diff) {
      return false;
    }
  }
  return true;
}

void rotationRemove(const AtomMatrix r1_passed, std::shared_ptr<Matter> m2) {
  AtomMatrix r1 = r1_passed;
  AtomMatrix r2 = m2->getPositions();

  // Align centroids
  VectorType c1(3);
  VectorType c2(3);

  c1[0] = r1.col(0).sum();
  c1[1] = r1.col(1).sum();
  c1[2] = r1.col(2).sum();
  c2[0] = r2.col(0).sum();
  c2[1] = r2.col(1).sum();
  c2[2] = r2.col(2).sum();
  c1 /= r1.rows();
  c2 /= r2.rows();

  for (int i = 0; i < r1.rows(); i++) {
    r1(i, 0) -= c1[0];
    r1(i, 1) -= c1[1];
    r1(i, 2) -= c1[2];

    r2(i, 0) -= c2[0];
    r2(i, 1) -= c2[1];
    r2(i, 2) -= c2[2];
  }

  RotationMatrix R = eonc::helper_functions::rotationExtract(r1, r2);

  // Eigen is transposed relative to numpy
  r2 = r2 * R;

  // Move centroid back to initial position
  for (int i = 0; i < r2.rows(); i++) {
    r2(i, 0) += c2[0];
    r2(i, 1) += c2[1];
    r2(i, 2) += c2[2];
  }

  m2->setPositions(r2);
  return;
}

void rotationRemove(const std::shared_ptr<Matter> m1,
                    std::shared_ptr<Matter> m2) {
  AtomMatrix r1 = m1->getPositions();
  rotationRemove(r1, m2);
  return;
}

void translationRemove(Matter &m1, const AtomMatrix r2_passed) {
  AtomMatrix r1 = m1.getPositions();
  AtomMatrix r2 = r2_passed;

  // net displacement
  VectorType disp(3);
  AtomMatrix r12 = m1.pbc(r2 - r1);

  disp[0] = r12.col(0).sum();
  disp[1] = r12.col(1).sum();
  disp[2] = r12.col(2).sum();
  disp /= r1.rows();

  for (int i = 0; i < r1.rows(); i++) {
    r1(i, 0) += disp[0];
    r1(i, 1) += disp[1];
    r1(i, 2) += disp[2];
  }

  m1.setPositions(r1);
  return;
}

void translationRemove(Matter &m1, const Matter &m2) {
  AtomMatrix r2 = m2.getPositions();
  translationRemove(m1, r2);
  return;
}

bool identical(const Matter &m1, const Matter &m2,
               const double distanceDifference) {

  AtomMatrix r1 = m1.getPositions();
  AtomMatrix r2 = m2.getPositions();

  std::set<int> matched;
  double tolerance = distanceDifference;

  if (r1.rows() != r2.rows()) {
    return false;
  }
  int N = r1.rows();

  for (int i = 0; i <= N; i++) {
    if (fabs((m1.pbc(r1.row(i) - r2.row(i))).norm()) < tolerance &&
        m1.getAtomicNr(i) == m2.getAtomicNr(i)) {
      matched.insert(i);
    }
  }

  for (int j = 0; j < N; j++) {

    if (matched.count(j) == 1)
      continue;

    for (int k = 0; k < N; k++) {
      if (matched.count(j) == 1)
        break;

      if (fabs((m1.pbc(r1.row(j) - r2.row(k))).norm()) < tolerance &&
          m1.getAtomicNr(j) == m2.getAtomicNr(k)) {
        matched.insert(j);
      }
    }

    // XXX: can we abort early if no match was found?
  }

  if (matched.size() == (unsigned)N) {
    return true;
  } else {
    return false;
  }
}

bool sortedR(const Matter &m1, const Matter &m2,
             const double distanceDifference) {
  SPDLOG_INFO("In sortedR");
  AtomMatrix r1 = m1.getPositions();
  AtomMatrix r2 = m2.getPositions();
  double tolerance = distanceDifference;
  int matches = 0;
  if (r1.rows() != r2.rows()) {
    return false;
  }

  // Allocate memory for rdf1 and rdf2
  std::vector<std::set<tmp::atom, tmp::by_atom>> rdf1(r1.rows());
  std::vector<std::set<tmp::atom, tmp::by_atom>> rdf2(r2.rows());

  for (int i2 = 0; i2 < r2.rows(); i2++) {
    rdf2[i2].clear();
    for (int j2 = 0; j2 < r2.rows(); j2++) {
      if (j2 == i2)
        continue;
      tmp::atom a2;
      a2.r = m2.distance(i2, j2);
      a2.z = m2.getAtomicNr(j2);
      rdf2[i2].insert(a2);
      rdf2[j2].insert(a2);
    }
  }
  for (int i1 = 0; i1 < r1.rows(); i1++) {
    if (matches == i1 - 2) {
      return false;
    }
    for (int j1 = 0; j1 < r1.rows(); j1++) {
      if (j1 == i1)
        continue;
      tmp::atom a;
      a.r = m1.distance(i1, j1);
      a.z = m1.getAtomicNr(j1);
      rdf1[i1].insert(a);
      rdf1[j1].insert(a);
    }
    for (int x = 0; x < r2.rows(); x++) {
      auto it2 = rdf2[x].begin();
      auto it = rdf1[i1].begin();
      int c = 0;
      int counter = 0;
      for (; c < r1.rows(); c++) {
        if (it == rdf1[i1].end() || it2 == rdf2[x].end())
          break;
        tmp::atom k1 = *it;
        tmp::atom k2 = *it2;
        if (fabs(k1.r - k2.r) < tolerance && k1.z == k2.z) {
          counter++;
        } else {
          SPDLOG_INFO("No match");
          break;
        }
        ++it;
        ++it2;
      }
      if (counter == r1.rows()) {
        matches++;
      } else {
        SPDLOG_INFO("No match");
      }
    }
  }

  return matches >= r1.rows();
}

void pushApart(std::shared_ptr<Matter> m1, double minDistance) {
  if (minDistance <= 0)
    return;

  AtomMatrix r1 = m1->getPositions();
  MatrixType Force(r1.rows(), 3);
  double f = 0.025;
  double cut = minDistance;
  double pushAparts = 500;
  for (int p = 0; p < r1.rows(); p++) {
    for (int axis = 0; axis <= 2; axis++) {
      Force(p, axis) = 0;
    }
  }
  for (int count = 0; count < pushAparts; count++) {
    int moved = 0;
    for (int i = 0; i < r1.rows(); i++) {
      for (int j = i + 1; j < r1.rows(); j++) {
        double d = m1->distance(i, j);
        if (d < cut) {
          moved++;
          for (int axis = 0; axis <= 2; axis++) {
            double componant = f * (r1(i, axis) - r1(j, axis)) / d;
            Force(i, axis) += componant;
            Force(j, axis) -= componant;
          }
        }
      }
    }
    if (moved == 0)
      break;
    for (int k = 0; k < r1.rows(); k++) {
      for (int axis = 0; axis <= 2; axis++) {
        r1(k, axis) += Force(k, axis);
        Force(k, axis) = 0;
      }
    }
    m1->setPositions(r1);
    // m1->matter2con("movie.con", true);
  }
}

// TODO(rg):: Doesn't belong here, should be in the jobs which need it
void saveMode(FILE *modeFile, std::shared_ptr<Matter> matter, AtomMatrix mode) {
  long const nAtoms = matter->numberOfAtoms();
  for (long i = 0; i < nAtoms; ++i) {
    if (matter->getFixed(i)) {
      fprintf(modeFile, "0 0 0\n");
    } else {
      fprintf(modeFile, "%lf\t%lf \t%lf\n", mode(i, 0), mode(i, 1), mode(i, 2));
    }
  }
  return;
}

std::vector<Matter> make_matter(const toml::table &config,
                                const std::shared_ptr<PotBase> &pot) {
  std::vector<Matter> mats;

  std::vector<std::string> input_files;
  if (auto inputs = config["Main"]["inputs"].as_array()) {
    for (const auto &input : *inputs) {
      if (input.is_string()) {
        input_files.push_back(input.as_string()->get());
      }
    }
  }

  // Process each input file
  eonc::mat::ConFileParser cfp;
  auto jtype = get_enum_toml<JobType>(config["Main"]["job"]);
  switch (jtype) {
  // All the single ended ones
  case JobType::Point:
    [[fallthrough]];
  case JobType::Minimization: {
    for (const auto &confile : input_files) {
      mats.emplace_back(Matter(pot));
      cfp.parse(mats.back(), confile);
    }
    return mats;
  }
  default: {
    throw std::runtime_error("No known job could be constructed");
  }
  }
}
} // namespace eonc::mat
