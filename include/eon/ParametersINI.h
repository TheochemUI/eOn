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

class INIReader;

namespace eonc {

class Parameters;

namespace config {

/// Parse all INI sections into the Parameters struct.
/// Returns 0 on success, nonzero on error.
int load_ini(INIReader &ini, Parameters &params);

/// Resolve cross-group defaults and time unit conversions.
/// Called after INI (or JSON) loading to finalize computed fields.
void validate_and_link(Parameters &params);

/// Copy an explicit [Optimizer] opt_method onto job sections that did not
/// set their own. NEB keeps LBFGS, and the dimer keeps CG, when the
/// optimizer key was omitted.
void inherit_job_opt_methods(Parameters &params);

} // namespace config
} // namespace eonc
