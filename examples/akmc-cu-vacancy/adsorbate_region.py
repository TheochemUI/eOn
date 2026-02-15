#!/usr/bin/env python3

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "ase",
#   "click",
#   "numpy",
# ]
# ///
"""
Identifies adsorbate atoms and nearby surface atoms in a structure file and
prints their 0-based indices to standard output.

This script is designed to be called by eOn's displacement atom list feature
(``displace_atom_kmc_state_script``). It selects atoms based on element type
and/or z-coordinate, then expands the selection to include all atoms within a
cutoff radius of the initial selection. This is useful for targeting saddle
search displacements at an adsorbate and its immediate surface environment.

Usage::

    uvx adsorbate_region.py pos.con --adsorbate-elements C O --cutoff 4.0
    uvx adsorbate_region.py pos.con --z-above 12.5 --cutoff 3.5
    uvx adsorbate_region.py pos.con --adsorbate-elements C O --z-above 10.0 --cutoff 4.0
"""

import logging
import sys

import ase.io as aseio
import click
import numpy as np

logging.basicConfig(
    level=logging.WARNING,
    format="%(levelname)s: %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)


def find_adsorbate_region(
    filename: str,
    adsorbate_elements: tuple[str, ...],
    z_above: float | None,
    cutoff: float,
) -> np.ndarray:
    """
    Identifies adsorbate atoms and their surface neighbours.

    Parameters
    ----------
    filename
        Path to a .con or other ASE-readable structure file.
    adsorbate_elements
        Chemical symbols of adsorbate atoms (e.g. ("C", "O")).
    z_above
        If set, also select atoms with z-coordinate above this value.
    cutoff
        Distance cutoff for expanding the selection to nearby atoms.

    Returns
    -------
    np.ndarray
        Sorted array of 0-based atom indices.
    """
    atoms = aseio.read(filename)
    atoms.set_pbc([True] * 3)
    n = len(atoms)
    symbols = np.array(atoms.get_chemical_symbols())

    # Build initial selection mask
    mask = np.zeros(n, dtype=bool)

    if adsorbate_elements:
        for elem in adsorbate_elements:
            mask |= symbols == elem
        log.info(
            f"Selected {mask.sum()} atoms by element: {adsorbate_elements}"
        )

    if z_above is not None:
        z_mask = atoms.positions[:, 2] > z_above
        mask |= z_mask
        log.info(
            f"Selected {z_mask.sum()} atoms with z > {z_above}"
        )

    if not mask.any():
        log.warning("No adsorbate atoms found with the given criteria.")
        print("")
        return np.array([], dtype=int)

    seed_indices = np.where(mask)[0]
    log.info(f"Seed atoms: {seed_indices.tolist()}")

    # Expand selection by distance cutoff
    if cutoff > 0:
        all_distances = atoms.get_all_distances(mic=True)
        expanded_mask = mask.copy()
        for idx in seed_indices:
            expanded_mask |= all_distances[idx] <= cutoff
        result = np.where(expanded_mask)[0]
    else:
        result = seed_indices

    log.info(f"Total selected atoms (with neighbours): {len(result)}")
    return result


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.argument(
    "filename",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
)
@click.option(
    "-e",
    "--adsorbate-elements",
    multiple=True,
    help="Chemical symbol(s) of adsorbate atoms (can be repeated).",
)
@click.option(
    "-z",
    "--z-above",
    type=float,
    default=None,
    help="Select atoms with z-coordinate above this value (Angstroms).",
)
@click.option(
    "-c",
    "--cutoff",
    type=float,
    default=4.0,
    show_default=True,
    help="Distance cutoff for expanding selection to nearby atoms (Angstroms).",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Enable verbose informational output to stderr.",
)
def main(
    filename: str,
    adsorbate_elements: tuple[str, ...],
    z_above: float | None,
    cutoff: float,
    verbose: bool,
):
    """
    Identifies adsorbate atoms and nearby surface atoms in FILENAME and prints
    their 0-based indices as a comma-separated list to stdout.

    Select adsorbate atoms by element (--adsorbate-elements) and/or by
    z-coordinate (--z-above). The selection is then expanded to include all
    atoms within --cutoff distance of any selected atom.
    """
    if verbose:
        log.setLevel(logging.INFO)

    if not adsorbate_elements and z_above is None:
        click.echo(
            "Error: specify at least one of --adsorbate-elements or --z-above.",
            err=True,
        )
        sys.exit(1)

    indices = find_adsorbate_region(filename, adsorbate_elements, z_above, cutoff)
    print(", ".join(map(str, indices)))


if __name__ == "__main__":
    main()
