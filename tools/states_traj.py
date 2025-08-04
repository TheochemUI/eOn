#!/usr/bin/env python3
"""
This script constructs an ASE trajectory file from EON state files.

It can generate trajectories in three ways:
1. From a continuous range of states.
2. From an explicit list of state numbers.
3. From an AKMC dynamics path defined in dynamics.txt.
"""

# /// script
# requires-python = ">=3.9"
# dependencies = [
#   "ase==3.23.0",
#   "pandas",
#   "rich==13.7.1",
#   "click",
# ]
# ///

import logging
import sys
from enum import Enum
from pathlib import Path

import click
import pandas as pd
from ase.io import read
from ase.io.trajectory import Trajectory
from rich.console import Console
from rich.logging import RichHandler

# --- Constants and Setup ---
# By convention, constants are in UPPER_SNAKE_CASE.
# Using pathlib.Path makes paths object-oriented and cross-platform compatible.
CWD = Path.cwd()
STATE_MAIN_DIR = CWD / "states"
STATE_TABLE_PATH = STATE_MAIN_DIR / "state_table"
DYNAMICS_FILE_PATH = CWD / "dynamics.txt"

# --- Setup for Rich library and Logging ---
# Use Rich for better terminal output and logging.
console = Console()
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=console, rich_tracebacks=True)],
)

# --- Data Structures ---
# Enums make the different modes of operation explicit and readable.
class TrajectoryMode(Enum):
    """Defines the different modes for trajectory generation."""
    CONTINUOUS = "continuous"
    EXPLICIT = "explicit"
    AKMC = "akmc"

# --- Helper Functions ---

def read_state_energies(path: Path) -> dict[int, float]:
    """Reads the state_table file and returns a dictionary of state energies."""
    logging.info(f"Reading state energies from '{path}'.")
    try:
        with open(path) as f:
            # This is a more robust way to parse the file
            return {
                int(pair[0]): float(pair[1])
                for line in f
                if (pair := line.strip().split(None, 1)) and len(pair) == 2
            }
    except FileNotFoundError:
        logging.error(f"Error: 'state_table' not found in '{path.parent}'.")
        sys.exit(1)
    except (ValueError, IndexError) as e:
        logging.error(f"Error parsing 'state_table': {e}")
        sys.exit(1)

def write_trajectory(output_filename: str, states_to_process: list[int]):
    """Writes a trajectory file from a list of state numbers."""
    logging.info(f"Generating trajectory for states -> '{output_filename}'")
    written_count = 0
    with Trajectory(output_filename, "w") as traj_writer:
        for state_num in states_to_process:
            reactant_file = STATE_MAIN_DIR / str(state_num) / "reactant.con"
            try:
                atoms = read(reactant_file)
                traj_writer.write(atoms)
                written_count += 1
            except FileNotFoundError:
                logging.warning(f"Skipping: reactant.con not found for state {state_num}.")
    if written_count > 0:
        console.print(f"\n[green]✔[/green] Trajectory with {written_count} states successfully written to '{output_filename}'.")
    else:
        console.print("\n[yellow]⚠[/yellow] No states were found to create a trajectory.")

def process_akmc_path(start_step: int, end_state: int, states_e: dict[int, float]) -> tuple[list[int], str]:
    """Processes dynamics.txt to find the path and write an energy profile."""
    logging.info(f"Processing AKMC path from step {start_step} to state {end_state}.")
    try:
        dynamics = pd.read_table(
            DYNAMICS_FILE_PATH,
            delimiter=r"\s+",
            skiprows=[0, 1],
            names=[
                "step-number", "reactant-id", "process-id", "product-id",
                "step-time", "total-time", "barrier", "rate", "energy",
            ],
        )
    except FileNotFoundError:
        logging.error(f"Error: '{DYNAMICS_FILE_PATH.name}' not found. Cannot use --akmc_step.")
        sys.exit(1)

    path_states, barriers = [], []
    selected_dynamics = dynamics[dynamics["step-number"] >= start_step]

    for _, row in selected_dynamics.iterrows():
        reactant_id = int(row["reactant-id"])
        product_id = int(row["product-id"])
        path_states.append(reactant_id)
        barriers.append(row["barrier"])
        if product_id in path_states:
            # Remove loops
            # Keep the state being looped back,
            # remove only the excursion
            index = path_states.index(product_id)
            del path_states[index + 1:]
            del barriers[index + 1:]
        if product_id == end_state:
            path_states.append(product_id)
            break
    else:
        logging.warning("Could not find a complete path to the end state.")

    # Write energy profile
    energy_profile_filename = f"energy_profile_{start_step}_{end_state}.dat"
    logging.info(f"Writing energy profile to '{energy_profile_filename}'")
    breakpoint()
    with open(energy_profile_filename, "w") as f:
        f.write("# State_ID    Energy(eV)    Barrier_Energy(eV)\n")
        for i, state_n in enumerate(path_states):
            energy = states_e.get(state_n, float('nan'))
            f.write(f"{int(state_n):<12d}  {energy:<15.4f}\n")
            if i < len(barriers):
                ts_energy = energy + barriers[i]
                f.write(f"{state_n + 0.5:<12.2f}  {ts_energy:<15.4f}\n")

    output_filename = f"akmc_path_{start_step}_to_{end_state}.traj"
    return path_states, output_filename

# --- Main Command-Line Interface ---
# Click provides a much cleaner and more powerful way to create CLIs.
@click.command()
@click.option('--states', type=int, multiple=True, help='A space-separated list of states to include.')
@click.option('--start', type=int, help='The starting state for a continuous range.')
@click.option('--end', type=int, help='The ending state for a continuous range.')
@click.option('--akmc_step', type=int, help='Starting step number from dynamics.txt for an AKMC path.')
@click.option('-o', '--output', type=str, default='', help='Name of the output trajectory file (default is auto-generated).')
def main(states, start, end, akmc_step, output):
    """Constructs an ASE binary trajectory (.traj) from EON states."""
    # --- Pre-computation Checks ---
    if not STATE_MAIN_DIR.is_dir():
        logging.error(f"Error: 'states' directory not found in '{CWD}'.")
        sys.exit(1)

    states_e = read_state_energies(STATE_TABLE_PATH)

    # --- Determine Mode and Process ---
    mode = None
    if start is not None and end is not None:
        mode = TrajectoryMode.CONTINUOUS
        states_to_process = list(range(start, end))
        output_filename = output or f"path_{start}_to_{end - 1}.traj"
    elif states:
        mode = TrajectoryMode.EXPLICIT
        states_to_process = list(states)
        output_filename = output or f"path_{states[0]}_{states[-1]}.traj"
    elif akmc_step is not None and end is not None:
        mode = TrajectoryMode.AKMC
        states_to_process, output_filename = process_akmc_path(akmc_step, end, states_e)
        if output:
            output_filename = output
    else:
        # Show help if no valid combination of arguments is given.
        click.echo(click.get_current_context().get_help())
        sys.exit(0)

    # --- Write Trajectory ---
    if states_to_process:
        write_trajectory(output_filename, states_to_process)
    else:
        logging.warning("No state path was generated.")

if __name__ == "__main__":
    main()
