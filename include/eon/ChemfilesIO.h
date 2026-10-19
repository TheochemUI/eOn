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
#include "eon/ConFileIO.h"

#include <string>
#include <string_view>
#include <vector>

namespace eonc {
class Matter;

namespace io {

/// True when this build links chemfiles (-Dwith_chemfiles found it).
[[nodiscard]] bool chemfiles_available() noexcept;

/// True for paths readcon owns: .con, .convel and their .gz/.zst forms.
[[nodiscard]] bool is_readcon_path(std::string_view filename);

/// True for structure paths routed to chemfiles: a known extension
/// (xyz, extxyz, pdb, cif, mmcif, gro, data, lmp, mol2, sdf, vasp, xtc, trr,
/// dcd, nc, arc, cml, cssr, mmtf, tng) or a POSCAR/CONTCAR name. Other dotted
/// names keep the classic meaning and get .con appended.
[[nodiscard]] bool is_chemfiles_path(std::string_view filename);

/**
 * [Main] input_format: empty or "con" keeps readcon only. Any other value
 * (xyz, pdb, cif, ...) makes a request for "pos.con" read "pos.<fmt>" when
 * that file exists, so jobs that name .con inputs start from any format.
 * Output stays .con only.
 */
void set_structure_input_format(std::string fmt);
[[nodiscard]] const std::string &structure_input_format() noexcept;

/**
 * Fallback cell for foreign files that carry none ([Main] cell / pbc).
 * `values` holds 3 numbers (orthorhombic a b c) or 9 (row-major lattice,
 * one vector per row). Empty clears it. Throws std::invalid_argument on
 * any other count.
 */
void set_config_cell(const std::vector<double> &values, bool pbc);
/// Parse "a b c" or nine numbers; throws std::invalid_argument.
[[nodiscard]] std::vector<double> parse_cell_values(const std::string &text);
/// Vacuum padding per side (Angstrom) for the bounding-box default cell.
void set_default_vacuum(double angstrom);
[[nodiscard]] double default_vacuum() noexcept;

/**
 * Input entry point: .con/.convel (and .gz/.zst) go to readcon, everything
 * else goes to chemfiles. An unknown extension is treated as .con, as before.
 */
[[nodiscard]] IoStatus inp2matter(Matter &m, const std::string &filename);

/**
 * Read the first frame of any chemfiles-supported file into Matter.
 *
 * Maps positions, atomic numbers, masses and velocities when present.
 * Fixed atoms come from per-atom extxyz properties: `fixed` (true = frozen)
 * or ASE's `move_mask` (true = free), as a bool, a number, or a 3-vector
 * for per-axis masks.
 *
 * Cell precedence: the file's own cell (periodic); else [Main] cell with
 * [Main] pbc; else a loud WARNING naming the file and both keys, and an
 * orthorhombic, non-periodic box spanning the atoms plus
 * [Main] default_vacuum (10 Angstrom) on every side.
 *
 * `format` overrides chemfiles' extension guess ("" = guess).
 */
/// Text of the last default-cell warning on this thread ("" when the last
/// read needed none). Lets tests and callers surface it.
[[nodiscard]] const std::string &last_cell_warning() noexcept;

[[nodiscard]] IoStatus chemfiles2matter(Matter &m, const std::string &filename,
                                        const std::string &format = "");

} // namespace io
} // namespace eonc
