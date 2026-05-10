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
#include "Eigen.h"
#include "GeometryAnalysis.h"
#include "Matter.h"
#include "RandomNumbers.h"
#include <string>
#include <vector>

namespace eonc {

/* Collection of supporting functions that handle arrays of doubles as vectors
 * and different random number generators */
namespace helpers {

inline constexpr double pi = 3.14159265358979323846;

// Backward-compatible wrappers delegating to eonc::rng
using eonc::rng::gaussRandom;
using eonc::rng::random;
using eonc::rng::randomDouble;
using eonc::rng::randomInt;

// Backward-compatible wrappers delegating to eonc::geometry
using eonc::geometry::identical;
using eonc::geometry::maxAtomMotion;
using eonc::geometry::maxAtomMotionApplied;
using eonc::geometry::maxAtomMotionAppliedV;
using eonc::geometry::maxAtomMotionV;
using eonc::geometry::maxMotionApplied;
using eonc::geometry::maxMotionAppliedV;
using eonc::geometry::numAtomsMoved;
using eonc::geometry::projectOutRotTrans;
using eonc::geometry::pushApart;
using eonc::geometry::rotationExtract;
using eonc::geometry::rotationMatch;
using eonc::geometry::rotationRemove;
using eonc::geometry::sortedR;
using eonc::geometry::translationRemove;

AtomMatrix makeOrthogonal(
    const AtomMatrix &v1,
    const AtomMatrix &v2); // return orthogonal component of v1 from v2
bool relaxMatter(Matter &matter, const Parameters &params, bool quiet = false,
                 bool writeMovie = false, bool checkpoint = false,
                 const std::string &prefixMovie = std::string(),
                 const std::string &prefixCheckpoint = std::string());
void getTime(double *real, double *user, double *sys);
bool existsFile(const std::string &filename); // does filename exist
std::string
getRelevantFile(std::string filename); // return filename containing _checkpoint
                                       // or _passed if such a file exists
VectorXd loadMasses(std::string filename, int nAtoms);
AtomMatrix loadMode(FILE *modeFile, int nAtoms);
AtomMatrix loadMode(std::string filename, int nAtoms);
void saveMode(FILE *modeFile, const std::shared_ptr<Matter> &matter,
              const AtomMatrix &mode);
void saveMode(const std::string &filename,
              const std::shared_ptr<Matter> &matter, const AtomMatrix &mode);
std::vector<int> split_string_int(const std::string &s,
                                  const std::string &delim);

} // namespace helpers

} // namespace eonc
