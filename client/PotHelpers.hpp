#pragma once

#include "C_Structs.h"

namespace eonc::pot {
// Typically this is done by the caller, however it is here as a sanity check
void zeroForceOut(const size_t &nAtoms, ForceOut *efvd);
} // namespace eonc::pot
