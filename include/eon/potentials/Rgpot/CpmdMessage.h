#pragma once

#include "eon/potentials/Rgpot/RGPotEngine.h"

#include <capnp/message.h>

#include "rgpot/rpc/Potentials.capnp.h"

// The CPMDParams root the engine passes to CPMDPot.
//
// A non-empty params_path is that message: functional, cutOffRy, charge,
// multiplicity, title, memory, inputSections, and inputBlocks. Scalar
// fields on the options are not written over a message that loaded.
// An empty path uses those scalars as the message.
//
// engine_path, engine_library, engine_root, scratch_dir, and permanent_dir
// are written after either source. They place the process.
//
// input_block, or RGPOT_CPMD_INPUT_BLOCK when that string is empty, is
// appended to inputBlocks. inputSections stay as the file stored them.
namespace eon {
::CPMDParams::Builder fillCpmdParams(::capnp::MallocMessageBuilder &msg,
                                     const RGPotEngineOptions &opt);
}
