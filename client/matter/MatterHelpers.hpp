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
#include "client/matter/Matter.h"

namespace eonc::mat {

bool rotationMatch(const Matter &m1, const Matter &m2, const double max_diff);
void rotationRemove(const AtomMatrix r1, std::shared_ptr<Matter> m2);
void rotationRemove(const std::shared_ptr<Matter> m1,
                    std::shared_ptr<Matter> m2);
void translationRemove(Matter &m1, const AtomMatrix r1);
void translationRemove(Matter &m1, const Matter &m2);

bool identical(const Matter &m1, const Matter &m2,
               const double distanceDifference);
bool sortedR(const Matter &m1, const Matter &m2,
             const double distanceDifference);
void pushApart(std::shared_ptr<Matter> m1, double minDistance);
void saveMode(FILE *modeFile, std::shared_ptr<Matter> matter, AtomMatrix mode);

} // namespace eonc::mat
