#!/usr/bin/env python
"""
This code is used to construct trajectory file with states
for help:
states_traj.py --h
"""

import sys
import os
import argparse
import pandas as pd
from ase.io import read
from ase.io.trajectory import Trajectory

# --- Script Setup ---
current = os.getcwd()
state_main_dir = os.path.join(current, "states")

# Check if states directory exists
if not os.path.isdir(state_main_dir):
    print(f"Error: 'states' directory not found in {current}")
    sys.exit(1)

# Read state energies into a dictionary
try:
    with open(os.path.join(state_main_dir, 'state_table')) as f:
        states_e = dict([int(pair[0]), float(pair[1])] for pair in [line.strip().split(None, 1) for line in f])
except FileNotFoundError:
    print(f"Error: 'state_table' not found in {state_main_dir}")
    sys.exit(1)

parser = argparse.ArgumentParser(description="Construct an ASE binary trajectory (.traj) from EON states.")
parser.add_argument('--states', type=int, nargs='+', metavar='StateNumber',
                    help='A space-separated list of states to include in the trajectory.')
parser.add_argument('--start', type=int, metavar='StateNumber',
                    help='The starting state number for a continuous trajectory.')
parser.add_argument('--end', type=int, metavar='StateNumber',
                    help='The ending state number for a continuous trajectory.')
parser.add_argument('--akmc_step', type=int, metavar='AKMCStepNumber',
                    help='Generate a trajectory based on dynamics.txt, starting from a specific AKMC step.')
parser.add_argument('-o', '--output', type=str, metavar='FILENAME', default='',
                    help='Name of the output trajectory file (default is auto-generated).')

args = parser.parse_args()

# --- Main Logic ---
output_filename = args.output
traj_writer = None

# Mode 1: A continuous range of states (--start and --end)
if args.start is not None and args.end is not None:
    if not output_filename:
        output_filename = f"path_{args.start}_to_{args.end-1}.traj"
    print(f"Generating trajectory from state {args.start} to {args.end-1} -> {output_filename}")

    traj_writer = Trajectory(output_filename, 'w')
    for state_num in range(args.start, args.end):
        try:
            atoms = read(os.path.join(state_main_dir, str(state_num), 'reactant.con'))
            traj_writer.write(atoms)
        except FileNotFoundError:
            print(f"Warning: reactant.con not found for state {state_num}. Skipping.")

# Mode 2: An explicit list of states (--states)
elif args.states:
    if not output_filename:
        output_filename = f"path_{args.states[0]}_{args.states[-1]}.traj"
    print(f"Generating trajectory from specified states: {args.states} -> {output_filename}")

    traj_writer = Trajectory(output_filename, 'w')
    for state_num in args.states:
        try:
            atoms = read(os.path.join(state_main_dir, str(state_num), 'reactant.con'))
            traj_writer.write(atoms)
        except FileNotFoundError:
            print(f"Warning: reactant.con not found for state {state_num}. Skipping.")

# Mode 3: A path from the AKMC dynamics log (--akmc_step and --end)
elif args.akmc_step is not None and args.end is not None:
    if not output_filename:
        output_filename = f"akmc_path_{args.akmc_step}_to_{args.end}.traj"
    print(f"Generating trajectory from AKMC step {args.akmc_step} to product state {args.end} -> {output_filename}")

    try:
        dynamics = pd.read_table('dynamics.txt', delimiter=r'\s+', skiprows=[0,1],
                                 names=['step-number', 'reactant-id', 'process-id', 'product-id', 'step-time', 'total-time', 'barrier', 'rate', 'energy'])
    except FileNotFoundError:
        print("Error: dynamics.txt not found. Cannot use --akmc_step.")
        sys.exit(1)

    rs = []
    barrier = []
    selected_dynamics = dynamics[dynamics['step-number'] >= args.akmc_step]

    for i in range(len(selected_dynamics)):
        reactant_id = selected_dynamics['reactant-id'].iloc[i]
        product_id = selected_dynamics['product-id'].iloc[i]
        rs.append(reactant_id)
        barrier.append(selected_dynamics['barrier'].iloc[i])
        try:
            index = rs.index(product_id)
            del rs[index:]
            del barrier[index:]
        except ValueError:
            pass
        if product_id == args.end:
            rs.append(product_id)
            break

    # Write energy profile data and build atoms list for trajectory
    energy_profile_filename = f"energy_profile_{args.akmc_step}_{args.end}.dat"
    print(f"Writing energy profile to {energy_profile_filename}")
    traj_writer = Trajectory(output_filename, 'w')
    with open(energy_profile_filename, 'w') as output:
        output.write("# State_ID    Energy(eV)    Barrier_Energy(eV)\n")
        for i, state_n in enumerate(rs):
            try:
                # Write state energy and barrier energy to file
                output.write("%-12d  %-15.4f\n" % (state_n, states_e[state_n]))
                output.write("%-12.2f  %-15.4f\n" % (state_n + 0.5, barrier[i] + states_e[state_n]))
            except IndexError:
                # Last state has no barrier, write its energy only
                output.write("%-12d  %-15.4f\n" % (state_n, states_e[state_n]))

            try:
                atoms = read(os.path.join(state_main_dir, str(state_n), 'reactant.con'))
                traj_writer.write(atoms)
            except FileNotFoundError:
                print(f"Warning: reactant.con not found for state {state_n}. Skipping.")

else:
    parser.print_help()
    sys.exit(0)

# --- Finalize ---
if traj_writer:
    traj_writer.close()
    print(f"\nTrajectory successfully written to '{output_filename}'.")
else:
    print("\nNo states were found to create a trajectory.")
