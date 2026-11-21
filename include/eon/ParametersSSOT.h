/*
** Cap'n Proto L0 field-graph defaults for covered Parameters groups.
** Generated defaults live in generated/ParametersSSOTDefaults.h (from
** schema/eon_params.capnp). INI/JSON remain adapters into Parameters.
*/
#pragma once

#include "Parameters.h"

namespace eonc::config {

/// Apply SSoT defaults for covered groups (Main, Potential, Optimizer,
/// Structure Comparison, Process Search, RgpotPot, Lineb). Safe to call
/// after construction; does not clear uncovered groups.
void apply_ssot_defaults(Parameters &params);

/// Return true if *section*/*key* (snake_case INI name) is declared in the
/// Cap'n Proto L0 catalog for covered groups.
bool ssot_has_field(const char *section, const char *key);

/// Validate projected [Lineb] values after INI or JSON load.
/// *energy_present* is true when the adapter saw an explicit energy key.
/// When job is lineb, energy must be present. beads must be 0 (path only)
/// or at least 2 (harmonic QTST); beads = 1 is refused.
void validate_lineb_options(const Parameters &params, bool energy_present);

} // namespace eonc::config
