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
// Precompiled header for the metatomic_pot translation unit.
//
// The torch + metatensor + metatomic include chain dominates the
// MetatomicPotential.cpp compile time -- ~2.5 GB of preprocessor
// output per build, ~30 s wall-clock on a warm cache, and several
// minutes from cold. PCH'ing the headers cuts the per-rebuild cost
// to a few seconds whenever the .cpp is touched without any of the
// torch / metatomic versions changing.
//
// Wrap pragmas mirror what MetatomicPotential.h does so the PCH and
// the consuming TU agree on which warnings the third-party headers
// silence; meson rebuilds the PCH whenever this header or any
// transitive include changes.

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"

#include <torch/cuda.h>
#include <torch/mps.h>
#include <torch/script.h>
#include <torch/version.h>
#include <torch/csrc/jit/runtime/graph_executor.h>

#include "metatensor/torch.hpp"
#include "metatensor/torch/module.hpp"
#include "metatomic/torch.hpp"

#pragma GCC diagnostic pop
