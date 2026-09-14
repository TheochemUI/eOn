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
#include <memory>

namespace eonc {
class Matter;

namespace geometry {

RotationMatrix rotationExtract(const AtomMatrix r1, const AtomMatrix r2);
bool rotationMatch(const Matter &m1, const Matter &m2, const double max_diff);
void projectOutRotTrans(Eigen::VectorXd &step, const AtomMatrix &positions);
void rotationRemove(const AtomMatrix r1, std::shared_ptr<Matter> m2);
void rotationRemove(const std::shared_ptr<Matter> m1,
                    std::shared_ptr<Matter> m2);
void translationRemove(Matter &m1, const AtomMatrix r1);
void translationRemove(Matter &m1, const Matter &m2);
double maxAtomMotion(const AtomMatrix v1);
double maxAtomMotionV(const VectorXd v1);
long numAtomsMoved(const AtomMatrix v1, double cutoff);
AtomMatrix maxAtomMotionApplied(const AtomMatrix v1, double maxMotion);
VectorXd maxAtomMotionAppliedV(const VectorXd v1, double maxMotion);
AtomMatrix maxMotionApplied(const AtomMatrix v1, double maxMotion);
VectorXd maxMotionAppliedV(const VectorXd v1, double maxMotion);

bool identical(const Matter &m1, const Matter &m2,
               const double distanceDifference);
bool sortedR(const Matter &m1, const Matter &m2,
             const double distanceDifference);
void pushApart(std::shared_ptr<Matter> m1, double minDistance);

} // namespace geometry
} // namespace eonc
