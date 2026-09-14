/*
** Cap'n Proto L0 field-graph defaults for covered Parameters groups.
** Generated defaults live in generated/ParametersSSOTDefaults.h (from
** schema/eon_params.capnp). INI/JSON remain adapters into Parameters.
*/
#pragma once

#include "Parameters.h"

namespace eonc::config {

/// Apply SSoT defaults for covered groups (Main, Potential, Optimizer,
/// Structure Comparison, Process Search). Safe to call after construction;
/// does not clear uncovered groups.
void apply_ssot_defaults(Parameters &params);

/// Return true if *section*/*key* (snake_case INI name) is declared in the
/// Cap'n Proto L0 catalog for covered groups.
bool ssot_has_field(const char *section, const char *key);

} // namespace eonc::config
