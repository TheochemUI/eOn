/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#pragma once

#include "BaseStructures.h"
#include "Job.h"

#include <functional>
#include <memory>

namespace eonc {

using JobFactory =
    std::function<std::unique_ptr<Job>(std::unique_ptr<Parameters>, Runtime)>;

void registerJob(JobType type, JobFactory factory);
// Pulls Job.cpp into static links so the registrars run (MSVC).
void forceJobRegistration();

} // namespace eonc
