#!/usr/bin/env python3

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "ase",
#   "click",
#   "numpy",
#   "ovito",
#   "rich",
# ]
# ///
"""
Identifies atoms in a structure file that do not match a specified crystal
structure (e.g., FCC) and prints their 0-based indices to standard output.

This script is designed to be called by eOn's displacement atom list feature
(``displace_atom_kmc_state_script``). It uses OVITO's Polyhedral Template
Matching (PTM) to classify each atom's local environment and returns the
indices of atoms that are *not* in the bulk crystal structure. These are
typically defect atoms (vacancies, interstitials, surfaces) that should be
preferentially displaced during saddle searches.

Usage::

    uvx ptmdisp.py pos.con
    uvx ptmdisp.py pos.con --structure-type FCC --verbose
"""

# 1. WARNING SUPPRESSION (before other imports)
import warnings

warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
warnings.filterwarnings(
    "ignore", category=FutureWarning, message=".*Please use atoms.calc.*"
)

# 2. IMPORTS
import logging
import sys
from enum import StrEnum

import ase.io as aseio
import click
import numpy as np
from ovito.io.ase import ase_to_ovito
from ovito.modifiers import PolyhedralTemplateMatchingModifier, SelectTypeModifier
from ovito.pipeline import Pipeline, StaticSource
from rich.logging import RichHandler

# 3. CONSTANTS and ENUMERATIONS
logging.basicConfig(
    level="WARNING",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
)
log = logging.getLogger(__name__)


class CrystalStructure(StrEnum):
    OTHER = "Other"
    FCC = "FCC"
    HCP = "HCP"
    BCC = "BCC"
    ICO = "Icosahedral"


STRUCTURE_TYPE_MAP = {
    CrystalStructure.FCC: PolyhedralTemplateMatchingModifier.Type.FCC,
    CrystalStructure.HCP: PolyhedralTemplateMatchingModifier.Type.HCP,
    CrystalStructure.BCC: PolyhedralTemplateMatchingModifier.Type.BCC,
    CrystalStructure.ICO: PolyhedralTemplateMatchingModifier.Type.ICO,
    CrystalStructure.OTHER: PolyhedralTemplateMatchingModifier.Type.OTHER,
}
STRUCTURE_PROPERTY_NAME = "Structure Type"


def find_mismatch_indices(
    filename: str, target_structure: CrystalStructure
) -> np.ndarray:
    """
    Analyzes a structure file with PTM and returns indices of atoms that
    do NOT match the target crystal structure.
    """
    try:
        log.info(f"Reading structure from '{filename}'...")
        atoms = aseio.read(filename)
        # NOTE: con readers in ase may lose PBC information
        atoms.set_pbc([True] * 3)
    except FileNotFoundError:
        log.critical(f"Error: The file '{filename}' was not found.")
        sys.exit(1)
    except Exception as e:
        log.critical(f"Failed to read or parse file '{filename}'. Error: {e}")
        sys.exit(1)

    pipeline = Pipeline(source=StaticSource(data=ase_to_ovito(atoms)))

    ptm = PolyhedralTemplateMatchingModifier()
    pipeline.modifiers.append(ptm)

    ovito_type = STRUCTURE_TYPE_MAP[target_structure]
    select_modifier = SelectTypeModifier(
        operate_on="particles",
        property=STRUCTURE_PROPERTY_NAME,
        types={ovito_type},
    )
    pipeline.modifiers.append(select_modifier)

    log.info(f"Running PTM analysis to find non-{target_structure.value} atoms...")
    data = pipeline.compute()

    # selection == 1 for matching atoms, 0 for non-matching
    mismatch_indices = np.where(data.particles.selection.array == 0)[0]
    log.info(f"Found {len(mismatch_indices)} non-{target_structure.value} atoms.")
    return mismatch_indices


# 4. CLI
@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument(
    "filename",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
)
@click.option(
    "-s",
    "--structure-type",
    "structure",
    type=click.Choice(CrystalStructure),
    default=CrystalStructure.FCC,
    show_default=True,
    help="The crystal structure to identify and exclude.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Enable verbose informational output to stderr.",
)
def main(filename: str, structure: CrystalStructure, verbose: bool):
    """
    Analyzes FILENAME to find all atoms that are NOT the specified
    crystal structure type and prints their 0-based indices as a
    comma-separated list, suitable for use in other programs.
    """
    if verbose:
        log.setLevel(logging.INFO)

    indices = find_mismatch_indices(filename, structure)
    print(", ".join(map(str, indices)))


if __name__ == "__main__":
    main()
