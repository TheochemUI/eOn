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
#include <optional>
#include <string>
#include <string_view>
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
    const AtomMatrix v1,
    const AtomMatrix v2); // return orthogonal component of v1 from v2
bool relaxMatter(Matter &matter, const Parameters &params, bool quiet = false,
                 bool writeMovie = false, bool checkpoint = false,
                 std::string prefixMovie = std::string(),
                 std::string prefixCheckpoint = std::string(),
                 std::vector<readcon::ConFrame> *outFrames = nullptr);
void getTime(double *real, double *user, double *sys);
bool existsFile(std::string filename); // does filename exist
std::string
getRelevantFile(std::string filename); // return filename containing _checkpoint
                                       // or _passed if such a file exists
VectorXd loadMasses(std::string filename, int nAtoms);
AtomMatrix loadMode(FILE *modeFile, int nAtoms);
AtomMatrix loadMode(std::string filename, int nAtoms);
// Load displacement.con, or synthesize initial + scale * unit(direction.dat)
// (issue #79). scale is saddle_search.displace_magnitude (default 0.1).
// Fixed-atom rows are always restored from initial. Returns false if neither
// displacement nor mode file is usable, or mode has zero norm when
// synthesizing.
bool loadOrSynthesizeDisplacement(Matter &target, const Matter &initial,
                                  const std::string &displacementPath,
                                  const std::string &modePath, double scale);
void saveMode(FILE *modeFile, std::shared_ptr<Matter> matter, AtomMatrix mode);
void saveMode(const std::string &filename, std::shared_ptr<Matter> matter,
              AtomMatrix mode);
std::vector<int> split_string_int(std::string s, std::string delim);

/// Display label for a force-convergence metric, or nullopt when the name is
/// none of the three the optimizer-driven searches understand. Single source
/// of the accepted spellings, shared by the INI loader and every search that
/// dispatches on the string.
std::optional<std::string_view> convergenceMetricLabel(std::string_view metric);
/// Throws std::invalid_argument naming context when metric is unrecognized.
/// Call at construction: the metric is compared once per optimizer step, and
/// a typo should not surface on iteration 4000.
void requireKnownConvergenceMetric(std::string_view metric,
                                   std::string_view context);

} // namespace helpers

} // namespace eonc
