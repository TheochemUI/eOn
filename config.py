##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

import ConfigParser
import numpy
import os.path
import sys
import string
import config

config.init_done = False

config.format = []


## New config section; does not affect current code, ignore for now. ===========

class ConfigSection:
    def __init__(self, name, description):
        self.name = name
        self.description = description
        self.keys = []
class ConfigKey:
    def __init__(self, name, kind, description, default):
        self.name = name
        self.kind = kind
        self.description = description
        self.default = default
        self.values = []
class ConfigValue:
    def __init__(self, name, description):
        self.name = name
        self.description = description

        
def fget(section, key = None):
    for f in format:
        if f.name == section:
            if key == None:
                return f
            for k in f.keys:
                if k.name == key:
                    return k
            return None
    return None

def fadd(section, key = None, value = None, description = "", kind = None, default = ""):
    if fget(section) is None:
        config.format.append(ConfigSection(section, description))
        return
    if fget(section, key) is None:
        s = fget(section)
        s.keys.append(ConfigKey(key, kind, description, default))
        return
    k = fget(section, key)
    k.values.append(ConfigValue(value, description))


# Main
fadd("Main", description = "These are the options that go in the 'Main' section of config.ini")
fadd("Main", "job", kind = "string", description = "The type of job to execute.", default = None)
fadd("Main", "job", "akmc", description = "Run an adaptive kinetic monte carlo simulation.")
fadd("Main", "job", "parallel_replica", description = "Calculate the rare-event dynamics of the system by combining transitions observed from multiple trajectories run in parallel.")
fadd("Main", "job", "process_search", description = "Combined saddle search, minimizations, and prefactor calculations. Used by the aKMC method.")
fadd("Main", "job", "saddle_search", description = "Do a saddle point search using a minimum mode method.")
fadd("Main", "job", "minimization", description = "Find the minimum from an initial configuration.")
fadd("Main", "job", "hessian", description = "Calculate the Hessian matrix for the specified configuration in a process.")
fadd("Main", "job", "dimer_dr", description = "Rye is changing this.")
fadd("Main", "job", "dimer_rotation", description = "Rye is changing this.")
fadd("Main", "job", "displacement_sampling", description = "Job to sample different displacement methods and parameters to see which are the most efficient.")
fadd("Main", "job", "basin_hopping", description = "Search for global minimum using basin hopping method.")
fadd("Main", "temperature", kind = "float", description = "The temperature that the job will run at.")
fadd("Main", "random_seed", kind = "int", description = "Takes an integer number for the random seed. If this number is less than zero the current time is used as the random seed.")
fadd("Main", "potential", kind = "string", description = "the type of potential to execute")
fadd("Main", "potential", "lj", description = "Lennard-Jones potential in reduced units.")
fadd("Main", "potential", "morse_pt", description = "Morse potential for platinum.")
fadd("Main", "potential", "emt", description = "Effective medium theory, for metals.")
fadd("Main", "potential", "edip", description = "Environment-Dependent Interatomic Potential, for carbon.")
fadd("Main", "potential", "vasp", description = "Vienna Ab-Initio Simulation Program (VASP) interface.")
fadd("Main", "potential", "tersoff_si", description = "Tersoff pair potential with angular terms, for silicon.")
fadd("Main", "potential", "sw_si", description = "Stillinger-Weber potential, for silicon.")
fadd("Main", "potential", "lenosky_Si", description = "Lenosky potential, for silicon.")
fadd("Main", "potential", "eam_al", description = "Embedded atom method parameterized for aluminum.")
fadd("Main", "potential", "qsc", description = "Quantum Sutton-Chen potential, for FCC metals.")
fadd("Main", "potential", "zpice", description = "Water on platinum.")
fadd("Main", "potential", "tip4p", description = "Point charge model for water.")
fadd("Main", "potential", "bopfox", description = "Bond order potential, for metals.")


# AKMC
fadd("AKMC", description = "Parameters for the AKMC section of config.ini.")
fadd("AKMC", "confidence", kind = "float", description = "The confidence (out of 1.0) criterion for moving to the next state.")
fadd("AKMC", "max_kmc_steps", kind = "int", description = "The maximum number of transitions per execution of the server.")
fadd("AKMC", "thermally_accessible_window", kind = "float", description = "Processes with barriers within this number of kT above the lowest barrier will be used in the rate table and for confidence calculations.")
fadd("AKMC", "thermally_accessible_buffer", kind = "float", description = "Processes with barriers of thermally_accessible_window + thermally_accessible_buffer will be stored, in the event that they are thermally accessible later, but are not used in the rate table or for the confidence calculations. Processes with barriers higher than the sum of these two values will be discarded.")
fadd("AKMC", "confidence_scheme", kind = "string", description = "")# add description
fadd("AKMC", "confidence_correction", kind = "boolean", description = "")# add description


# Structure Comparison
fadd("Structure Comparison", description = "")
fadd("Structure Comparison", "energy_difference", kind = "float", description = "How close in energy two configurations must be to be considered energetically equivalent.")
fadd("Structure Comparison", "distance_difference", kind = "float", description = "The maximum distance two mapped atoms may be for two configurations to be considered equivalent.")
fadd("Structure Comparison", "indistinguishable_atoms", kind = "boolean", description = "Use an algorithm to compare structures that does not distinguish between atoms of the same element. That is to say the numbering of the atoms does not affect the structural comparison.")
fadd("Structure Comparison", "check_rotation", kind = "boolean", description = "Finds optimal overlap of structures via rotation before comparing them. Use this option in systems where structures can become rotated, such as nanoparticles.")
fadd("Structure Comparison", "neighbor_cutoff", kind = "float", description = "Atoms within this distance of each other are considered neighbors.")
fadd("Structure Comparison", "use_covalent", kind = "boolean", description = "Use the covalent radii of atoms to determine neighbors.")
fadd("Structure Comparison", "covalent_scale", kind = "float", description = "Multiply covalent radii by this amount before determining neighbors.")
fadd("Structure Comparison", "brute_neighbors", kind = "boolean", description = "Determine neighbors by brute force (use this with nonorthogonal boxes).")


# Paths
fadd("Paths", description = "Location of files related to the sending and receiving data on the server.")
fadd("Paths", "main_directory", kind = "string", description = "This is the root directory of the simulation. Configuration files and the initial reactant are here and by default all of the simulation data will be stored under this directory.")
fadd("Paths", "searches_in", kind = "string", description = "")
fadd("Paths", "states", kind = "string", description = "Where all of the information about individual states is located.")
fadd("Paths", "scratch", kind = "string", description = "")
fadd("Paths", "potential_files", kind = "string", description = "For extra files needed by the client for the potential.")
fadd("Paths", "jobs_out", kind = "string", description = "")#add description
fadd("Paths", "jobs_in", kind = "string", description = "")# add description
fadd("Paths", "incomplete", kind = "string", description = "")# add description
fadd("Paths", "results", kind = "string", description = "")#add description
fadd("Paths", "bh_minima", kind = "string", description = "")#add descripiton
fadd("Paths", "kdb_scratch", kind = "string", description = "")#add description
fadd("Paths", "kdb", kind = "string", description = "")#add description
fadd("Paths", "superbasin_recycling", kind = "string", description = "")#add description



# Process Search
fadd("Process Search", description = "The AKMC method can ask clients to do a saddle search, find connecting minima, and calculate prefactors all within the Process Search job type.")
fadd("Process Search", "minimize_first", kind = "boolean", description = "Every time a process search is run by a client the reactant will be minimized first before doing any saddle searches.")
fadd("Process Search", "prefactor_min", kind = "float", description = "Minimum value for a resonable prefactor.")
fadd("Process Search", "prefactor_max", kind = "float", description = "Maximum value for a resonable prefactor.")
fadd("Process Search", "default_prefactor", kind = "int", description = "Calculate prefactor if zero, otherwise use given value instead of doing a full prefactor calculation.")
fadd("Process Search", "minimization_offset", kind = "float", description = "After a saddle is found, images are placed on either side of the saddle along the mode and minimized to ensure that the saddle is connected to the original minimum and to locate the product state. This is the distance those images are displaced from the saddle.")
fadd("Saddle Search", description = "A saddle search is initiated by making a local displacement of atoms from their position at the minimum of the current state. This displacement can be done using the different strategies indicated by the displace_type option, and the following parameters. If the user knows something about the local environment where reactions are likely to take place in the system, this information can be used to make saddle searches more efficient by getting them started in the right part of configuration space.")


# Saddle Search
fadd("Saddle Search", "displace_type", kind = "string", description = "Type of displace to use")
fadd("Saddle Search", "displace_type", "random", description = "Select an atom at random from the free atoms in the configuration.")
fadd("Saddle Search", "displace_type", "least_coordinated", description = "Determine the lowest coordination number of all atoms in the configuration and select one atom at random with that coordination number.")
fadd("Saddle Search", "displace_type", "under_coordinated", description = "Select a random atom with coordination less than displace_max_coordination.")
fadd("Saddle Search", "displace_type", "listed_atoms", description = "Select an atom from the list passed as displace_atomlist")
fadd("Saddle Search", "displace_type", "not_FCC_HCP", description = "Select an atom that is not in a local HCP or FCC coordination.")
fadd("Saddle Search", "displace_type", "client_least_coordinated", description = "Displacement is made on the client centered on an atom with the minimal coordination.")
fadd("Saddle Search", "displace_type", "cleint_not_FCC_HCP_coordinated", description = "Displacement is made on the client centered on an atom that is not in a local HCP or FCC coordination.")
fadd("Saddle Search", "displace_type", "client_last_atom", description = "Displacement is made on the client centered on the last atom in the configuration.")
fadd("Saddle Search", "max_iterations", kind = "int", description = "The maximum number of translation steps to be taken.")
fadd("Saddle Search", "displace_radius", kind = "float", description = "Atoms within this distance of the epicenter will be displaced.")
fadd("Saddle Search", "displace_magnitude", kind = "float", description = "The standard deviation of the magnitude of the displacement.")
fadd("Saddle Search", "displace_min_norm", kind = "float", description = "The total length of the displacement vector is ensured to exceeds this value. Is useful when only a few degrees of freedoms are displaced to guarantee that the starting point for the saddle point search is significantly different from the initial minimum.")
fadd("Saddle Search", "displace_max_coordination", kind = "int", description = "When using under_coordinated as the displacement type, choose only atoms with a coordination equal to or less than this.")
fadd("Saddle Search", "min_mode_method", kind = "string", description = "")
fadd("Saddle Search", "min_mode_method", "dimer", description = "Use the dimer min-mode method.")
fadd("Saddle Search", "min_mode_method", "lanczos", description = "Use the Lanczos min-mode method.")
fadd("Saddle Search", "max_energy", kind = "float", description = "The energy at which a saddle search is considered bad and terminated.")
fadd("Saddle Search", "displace_atomlist", kind = "string", description = "The individual index should be seperated by a comma 10, 20,-1 would be the 10, 20 and the last atom.")
fadd("Saddle Search", "client_max_single_displace", kind = "int", description = "Only functional when displacement is done on the client. Defines the maximal allowed absolute value for a single component in the displacement vector.")
fadd("Saddle Search", "stdev_translation", kind = "float", description = "")#add description
fadd("Saddle Search", "stdev_rotation", kind = "float", description = "")#add description
fadd("Saddle Search", "molecule_list", kind = "string", description = "")#add description
fadd("Saddle Search", "disp_at_random", kind = "int", description = "")#add description


# Dimer
fadd("Dimer", description = "Options for controlling the dimer minimum mode finding method on the client side.")
fadd("Dimer", "seperation", kind = "float", description = "Seperation between dimer images.")
fadd("Dimer", "finite_angles", kind = "float", description = "Finite difference angle over which the dimer is rotated to find the lowest curvature.")
fadd("Dimer", "rotations_min", kind = "int", description = "Minimum number of rotations allowed for the dimer in each step.")
fadd("Dimer", "rotations_max", kind = "int", description = "Maximum number of rotations allowed for the dimer in each step.")
fadd("Dimer", "torque_min", kind = "float", description = "Minimum torque above which the dimer rotates only once and below which is does not rotate.")
fadd("Dimer", "torque_max", kind = "float", description = "Maximum torque above which the dimer rotates up to rotations_max times.")
fadd("Dimer", "improved", kind = "boolean", description = "Improvements to the dimer method from Kastner.")
fadd("Dimer", "coverged_rotatoin", kind = "float", description = "Dimer is considered converged if it rotates fewer degrees than this.")
fadd("Dimer", "opt_method", kind = "string", description = "Optimization algorithm to choose the dimer rotation direction")
fadd("Dimer", "opt_method", "sd", description = "steepest descent, rotate along the rotational force.")
fadd("Dimer", "opt_method", "cg", description = "conjudate gradient, rotate along conjugate directions.")
fadd("Dimer", "opt_method", "lbfgs", description = "quasi-Newton method [not implemeted yet].")


# Lanczos
fadd("Lanczos", description = "Options for controlling the Lanczos minimum mode finding method on the client side.")
fadd("Lanczos", "max_iterations", kind = "int", description = "The maximum number of refinement iterations when calculating the minimum eigenvalue.")
fadd("Lanczos", "tolerance", kind = "float", description = "This is the convergence critera for relative error of the lowest eigenvalue.")
fadd("Lanczos", "finite_dist", kind = "float", description = "Finite difference step size when computing the second derivative of the potential along the Lanczos vectors.")


# Hessian
fadd("Hessian", description = "Options for controlling how Hessian matricies are calculated.")
fadd("Hessian", "type", kind = "string", description = "The Hessian to be calculated has to be one of reactant, saddle, or product.")
fadd("Hessian", "finite_dist", kind = "float", description = "Finite difference distance between forces used to construct the Hessian.")
fadd("Hessian", "min_displacement", kind = "float", description = "Minimum amount that an atom has to move to be included in the Hessian calculation.")
fadd("Hessian", "within_radius", kind = "float", description = "Atoms within this radius of moving atoms are included in the Hessian.")


# Communicator
fadd("Communicator", description = "Options that apply to all of the different communicator types.")
fadd("Communicator", "type", kind = "string", description = "")
fadd("Communicator", "type", "local", description = "The local communicator runs the calculations on the same computer that the server is run on.")
fadd("Communicator", "type", "cluster", description = "A job scheduler can be used to run jobs through user supplied shell scripts. Examples are given for SGE.")
fadd("Communicator", "type", "boinc", description = "Jobs can be submitted to a BOINC project.")
fadd("Communicator", "type", "arc", description = "Jobs can be submitted to the grid computing software ARC.")
fadd("Communicator", "num_jobs", kind = "int", description = "Local( The number of jobs that will be run every time the program is invoked) Cluster( The desired sum of the queued and running jobs.) Boinc( The number of jobs to keep in the queue.")
fadd("Communicator", "jobs_per_bundle", kind = "int", description = "In eon a job is defined as task that the eon client executes, such as a process search or a parallel replica run. Sometimes it makes sense to run more than one of the same type of job at a time.")
fadd("Communicator", "client_path", kind = "string", description = "Either the name or path to the eon client binary. If only a name and not a path is given then eon looks for the binary in same directory as config.ini failing to find it there it will search though the directories in the $PATH environment variable.")
fadd("Communicator", "number_of_CPUs", kind = "int", description = "The number of jobs that will run simultaneously.")
fadd("Communicator", "script_path", kind = "string", description = "The path to the user defined scripts for submitting jobs to the communicator.")
fadd("Communicator", "name_prefix", kind = "string", description = "When jobs are submitted to the scheduler they are given a unique internally used named. In order to make the jobs identifiable by the user the name_prefix can be set to a meaningful string that will always be prepended to the job names.")
fadd("Communicator", "queued_jobs", kind = "string", description = "This is the name of the script that returns the job ids of all the running and queued jobs. It does not have to return the job ids of only eon related jobs.")
fadd("Communicator", "submit_job", kind = "string", description = "This is the name of the script that submits a single job to the queuing system. It takes two command line arguments. The first is the name of the job. This is not required for eon use, but is highly recommended so that users can identify which job is which. The second argument is the working directory. This is the path where the eon client should be executed. All of the needed client files will be placed in this directory. The script must return the job id of the submitted job. This is how eon internally keeps track of jobs.")
fadd("Communicator", "cancel_job", kind = "string", description = "This is the name of the script that cancels a job. It takes a single argument the job id.")
fadd("Communicator", "boinc_project_dir", kind = "string", description = "This is the full path to the root of the BOINC project directory.")
fadd("Communicator", "boinc_wu_template_path", kind = "string", description = "This is the path, relative from the boinc_project_dir, to the boinc workunit template.")
fadd("Communicator", "boinc_re_template_path", kind = "string", description = "This is the path, relative from the boinc_project_dir, to the boinc result template.")
fadd("Communicator", "boinc_appname", kind = "string", description = "This is the name of the application in BOINC.")
fadd("Communicator", "boinc_results_path", kind = "string", description = "This is the path where BOINC puts the final results. If you are using the sample_assimilator the results are stored in the project directory in a folder named sample_results.")
fadd("Communicator", "blacklist", kind = "string", description = "")#add description
fadd("Communicator", "boinc_priority", kind = "int", description = "The priority of the BOINC workunits that will be submitted.")

# Parallel Replica
fadd("Parallel Replica", description = "Parallel Replica Dynamics (PRD) is the simplest and the most accurate way to do accelerated-MD simulation. The only assumption made in this method is that the reactions satisfy first order kinetics.")
fadd("Parallel Replica", "time_step", kind = "float", description = "The length of each MD step in femtoseconds.")
fadd("Parallel Replica", "auto_stop", kind = "boolean", description = "Whether or not stop the job when a new state is found. For boinc communicator this value should be set to false.")
fadd("Parallel Replica", "steps", kind = "int", description = "The number of MD steps to run.")
fadd("Parallel Replica", "dephase_steps", kind = "int", description = "Number of steps used to decorrelate the replica trajectories. The momenta will be inversed when reaching the dividing surface to prevent transitions occurring during this period.")
fadd("Parallel Replica", "check_period", kind = "int", description = "How frequently the state of system is checked. Every check_period steps, the current structure and the initial one will be compared to tell whether a newstate has been reached. Also note when you set refine as true, the code will keep a buffer array consisting of check_period/record_resolution+1 atomic configurations, which may increase the usage of memory.")
fadd("Parallel Replica", "refine_transition_time", kind = "boolean", description = "Whether or not the transition time is refined. When this option is turned on, the code will keep an array consisted by check_period/record_resolution+1 atomic configurations. A Binary search algorithm is employed to determine the transition step. Otherwise the transition step would be the first in which a new state was found. This function reduces the need for a smaller check_period. And the accuracy of transition time is record_resolution*timestep.")
fadd("Parallel Replica", "record_resolution", kind = "int", description = "How often the system is recorded to the buffer array when the refine_transition_time option is activated. Increasing the value of record_resolution lowers the accuracy of the transition time estimate but also reduces memory usage and speeds up refinement of the transition step.")
fadd("Parallel Replica", "post_transition_steps", kind = "int", description = "Number of MD steps which will be performed after a new state has been found. A state check will be employed after these post_transition_steps to confirm that the state is stable. This additional check helps avoid meta-stable states. A value similar to dephase_steps is recommended.")
fadd("Parallel Replica", "thermo_type", kind = "string", description = "")
fadd("Parallel Replica", "thermo_type", "andersen", description = "Andersen thermostat with Verlet algorithm")
fadd("Parallel Replica", "thermo_type", "nose_hoober", description = "Nose-Hover thermostat with Verlet algorithm")
fadd("Parallel Replica", "thermo_type", "langevin", description = "Langevin thermostat with Verlet algorithm")
fadd("Parallel Replica", "andersen_alpha", kind = "float", description = "The collision strength in the Andersen thermostat")
fadd("Parallel Replica", "andersen_collision_steps", kind = "float", description = "The collision period (in MD steps) for the Andersen thermostat.")
fadd("Parallel Replica", "nose_mass", kind = "float", description = "The effective mass of the additional degree of freedom in the Nose-Hover thermostat, which determines the rate of heat transfer.")
fadd("Parallel Replica", "langevin_friction", kind = "float", description = "The damping coefficient for langevin dynamics.")
fadd("Parallel Replica", "bias_potential", kind = "string", description = "")
fadd("Parallel Replica", "bias_potential", "none", description = "with no bias potential, run regular MD")
fadd("Parallel Replica", "bias_potential", "bond_boost", description = "bond boost method from Miron and Fichthorn")
fadd("Parallel Replica", "bb_dvmax", kind = "float", description = "The magnitude of the bond-boost bias potential. It should be smaller than the barrier of any transition.")
fadd("Parallel Replica", "bb_rmd_steps", kind = "int", description = "Number of MD steps used to determine the equilibrium bond length before the bias potential is added.")
fadd("Parallel Replica", "bb_stretch_threshold", kind = "float", description = "Defines the bond-boost dividing surface. It should be smaller than the maximum fractional nearest-neighbor bond stretch or compression at any transition state.")
fadd("Parallel Replica", "bb_ds_curvature", kind = "float", description = "The curvature near the bond-boost dividing surface, it should has a value <= 1. We recommend the value to be 0.9-0.98.")
fadd("Parallel Replica", "bb_rcut", kind = "float", description = "All bonds which belong to the tagged atoms and are shorter than a cutoff of rcut will be included in the bond-boost potential.")


# Basin Hopping
fadd("Basin Hopping", description = "Basin hopping is a Monte Carlo method in which the energy of each configuration is taken to be the energy of a local minimization.")
fadd("Basin Hopping", "steps", kind = "int", description = "Number at steps to take at the temeprature assigned.")
fadd("Basin Hopping", "max_displacement", kind = "float", description = "Max displacement in each degree of freedom.")
fadd("Basin Hopping", "max_displacement_algorithm", kind = "string", description = "The algorithm used to assign max displacement of each atom.")
fadd("Basin Hopping", "max_displacement_algorithm", "standard", description = "The max displacement of all the atoms will be the value assigned in max_displacement.")
fadd("Basin Hopping", "max_displacement_algorithm", "linear", description = "The max displacement of each atom will be linearly correlated to its distance from the geometric center, with the overal maximum displacement being the value assigned in max_displacement.")
fadd("Basin Hopping", "max_displacement_algorithm", "quadratic", description = "The max displacement of each atom will be quadratically correlated to its distance from the geometric center, with the overal maximum displacement being the value assigned in max_displacement.")
fadd("Basin Hopping", "displacement_distribution", kind ="string", description = "The distribution used for the displacement of each atom.")
fadd("Basin Hopping", "displacement_distribution", "uniform", description = "A random number is selected between 0 and max_displacement")
fadd("Basin Hopping", "displacement_distribution", "gaussian", description = "max_displacement serves as the width of a gaussian distribution used to select displacements")
fadd("Basin Hopping", "quenching_steps", kind = "int", description = "Number at steps at 0 temperature.")
fadd("Basin Hopping", "single_atom_displace", kind = "boolean", description = "Displace only one atom per step.")
fadd("Basin Hopping", "stay_minimized", kind = "boolean", description = "Displace minimized structures.")
fadd("Basin Hopping", "swap_probability", kind = "float", description = "The probability in range [0,1.0] that a swapping step takes place instead of a displacement step. The swap step selects two atoms of different elements and switches them.")
fadd("Basin Hopping", "jump_max", kind = "int", description = "The number of consecutive rejected steps after which jump steps should be taken. This serves to provide a more global search when the structure is stuck in a certain basin. The number of jump steps is assigned in jump_steps. See paper on the Basin Hopping with Occasional Jumping algorithm by Iwamatsu and Okabe.")
fadd("Basin Hopping", "jump_steps", kind = "int", description = "The number of jump steps to take after the jump_max number of consecutive rejections have taken place.")
fadd("Basin Hopping", "md_probability", kind = "float", description = "The probability that molecular dynamics is run before the first step is taken. The parameters for molecular dynamics can be set in the [Dynamics] section. As default, the temperature specified in [Main] will be used, or another temperature can be specificed in MD_temp.")
fadd("Basin Hopping", "md_temp", kind = "float", description = "The temperature to use for the molecular dynamics steps that occur before basin hopping steps.")
fadd("Basin Hopping", "md_probability", kind = "boolean", description = "The probability that a basin hopping job will run Molecular Dynamics before starting BH.")


# Optimizers
fadd("Optimizers", description = "The eon client provides two optimizers quick-min and conjugate gradients. These options are in the [Optimizers] section.")
fadd("Optimizers", "opt_method", kind = "string", description = "The optimization method to use.")
fadd("Optimizers", "opt_method", "cg", description = "Conjugate gradient")
fadd("Optimizers", "opt_method", "qm", description = "Quickmin")
fadd("Optimizers", "opt_method", "box", description = "Optimizes the atom positions and box using quickmin")
fadd("Optimizers", "max_iterations", kind = "int", description = "The maximum number of optimization iterations that the will be performed.")
fadd("Optimizers", "converged_force", kind = "float", description = "When the maximum force (in eV/A) on any one atom is smaller than this value, the structure is considered minimized.")
fadd("Optimizers", "max_move", kind = "float", description = "Maximum distance that an atom may be moved in a single optimization step.")
fadd("Optimizers", "time_step", kind = "float", description = "The dynamical timestep for the quickmin algorithm.")
fadd("Optimizers", "finite_dist", kind = "float", description = "The finite difference step size.")


#Corase Graining
fadd("Coarse Graining", description = "In AKMC simulations where there are vastly different rates, the simulation can get stuck in a group of states connected by relatively fast rates. In order to explore slower transitions, a prohibitively large number of KMC steps may be needed. In order to circumvent this problem, eOn implements two methods. The first method, projective dynamics [1], groups states that are joined by fast rates into superbasins. Information about transitions between states in a superbasin is lost, but the rates for transitions across a superbasin are correct. The second method, accelerated superbasining kinetic Monte Carlo (AS-KMC) [2], artificially raises low barriers. The dynamics between states connected by fast rates are simulated, but an error is introduced in the dynamics direction and time. Both methods cannot be used simultaneously.")
fadd("Coarse Graining", "state_file", kind = "string", description = "File name for the state specific data stored within each of the state directories.")
fadd("Coarse Graining", "use_projective_dynamics", kind = "boolean", description = "This option determines whether the projective dynamics coarse graining method will be used. This mutually excludes the use_askmc option.")
fadd("Coarse Graining", "superbasin_scheme", kind = "string", description = "Projective dynamics provides a method for calculating transition rates across superbasins. An additional method is needed in order to decide when to combine states into a superbasin. eOn provides two methods. The first method, called transition counting, counts the number of times that the simulation has transitioned between a given pair of states. After a critical number of transitions have occured, the pair of states are merged to form a superbasin. (If one is already in a superbasin, the other is added to that superbasin. If both are already in superbasins, the two superbasins are merged). This method can be selected by setting scheme equal to transition_counting. Jean-Claude: document your scheme here This method can be elected by setting scheme equal to energy_level.")
fadd("Coarse Graining", "number_of_transitions", kind = "int", description = "If the transition counting scheme is being used (scheme=transition_counting), this is the number of transitions that must occur between two states before they are merged into a superbasin.")
fadd("Coarse Graining", "energy_increment", kind = "float", description = "If the energy level scheme is being used (scheme=energy_level). Each state, the first time it is visited, is assigned an energy level first equal to the energy of the minimum. Every time the state is visited again by the Monte Carlo simulation, the energy level is increased by this amount")
fadd("Coarse Graining", "use_askmc", kind = "boolean", description = "This option determines whether the AS-KMC coarse graining method will be used. This mutually excludes the use_projective_dynamics option.")
fadd("Coarse Graining", "askmc_condifence", kind = "float", description = "The confidence for AS-KMC. This value determines the accuracy of the direction of the dynamics trajectory.")
fadd("Coarse Graining", "askmc_barrier_raise_param", kind = "float", description = "This parameter sets how much the barriers are raised during AS-KMC. ( in the reference.)")
fadd("Coarse Graining", "askmc_high_barrier_def", kind = "int", description = "This parameter sets how high a barrier must be to be considered high in AS-KMC.")
fadd("Coarse Graining", "askmc_barrier_test_on", kind = "boolean", description = "")
fadd("Coarse Graining", "askmc_connections_test_on", kind = "boolean", description = "This parameter determines whether to ensure that there are no processes which connect states in the defined superbasin which have not been visited yet and which have a low-barrier. This check is somewhat more computationally expensive than the previous because structure comparisons have to be made when finding product states of unvisited processes.")


# KDB
fadd("Kdb", description = "One of the bottlenecks in an aKMC simulation is performing the saddle point searches. The kinetic database is used to ameliorate this cost by storing information about processes as they are found and using it to predict future saddle points.")
fadd("Kdb", "use_kdb", kind = "boolean", description = "Turn Kdb on/off.")
fadd("Kdb", "wait", kind = "boolean", description = "Wait for the query to finish before submitting jobs (for debugging purposes).")
fadd("Kdb", "keep", kind = "boolean", description = "Keep the saddle suggestions (for debugging purposes).")
fadd("Kdb", "Kdb_only", kind = "boolean", description = "")#add description
fadd("Kdb", "addpath", kind = "string", description = "")#add description
fadd("Kdb", "querypath", kind = "string", description = "")#add description

# Recycling         add descriptions
fadd("Recycling", description = "")
fadd("Recycling", "use_recycling", kind = "boolean", description = "")
fadd("Recycling", "save_suggestions", kind = "boolean" , description = "")
fadd("Recycling", "displace_moved_only", kind = "boolean",  description = "")
fadd("Recycling", "move_distance", kind = "float",  description = "")
fadd("Recycling", "use_sb_recycling", kind = "boolean", description = "")

# Debug
fadd("Debug", description = "Parameters that are generally used to help debug calculations")
fadd("Debug", "keep_bad_saddles", kind = "boolean", description = "Keep data about bad saddles. If true, the result files for failed saddle searches are kept in the badprocdata directory within the state directory for that search.")
fadd("Debug", "keep_all_results_files", kind = "boolean", description = "Stores all result files in main_directory/results")
fadd("Debug", "result_files_path", kind="string", description="Where to store all result files. Defaults to 'debug_results'.")
fadd("Debug", "register_extra_results", kind = "boolean", description = "Register processes found for a state after leaving that state.")
fadd("Debug", "use_mean_time", kind = "boolean", description = "Select transition times from the mean of the exponential distribution of escape times.")
fadd("Debug", "target_trajectory", kind = "boolean", description = "Follow the state-to-state trajectory of another akmc simulation.")
fadd("Debug", "save_stdout", kind = "boolean", description = "Save the standard output from the client to a file named stdout_0.dat")
fadd("Debug", "interactive_shell", kind = "boolean", description = "")#add description

## End new config section. =====================================================

def init(config_file = ""):
    if config.init_done:
        return None
    config.init_done = True

    parser = ConfigParser.SafeConfigParser()

    parser.read(os.path.join(sys.path[0], 'default_config.ini'))

    gave_config = True
    if config_file != "":
        if os.path.isfile(config_file):
            parser.read(config_file)
            config.config_path = os.path.abspath(config_file)
        else:
            print >> sys.stderr, "specified configuration file %s does not exist" % config_file
            sys.exit(2)
    elif os.path.isfile('config.ini'):
        parser.read('config.ini')
        config.config_path = os.path.abspath('config.ini')
        gave_config = False
    else:
        print >> sys.stderr, "You must provide a configuration file either by providing it as a command line argument or by placing a config.ini in the current directory."
        sys.exit(2)        

    #Main options
    config.main_job = parser.get('Main', 'job')
    config.main_temperature = parser.getfloat('Main', 'temperature')

    try:
        config.main_random_seed = parser.getint('Main', 'random_seed')
        numpy.random.seed(config.main_random_seed)
    except:
        config.main_random_seed = None

    #Structure Comparison options
    config.comp_eps_e = parser.getfloat('Structure Comparison', 'energy_difference')
    config.comp_eps_r = parser.getfloat('Structure Comparison', 'distance_difference')
    config.comp_use_identical = parser.getboolean('Structure Comparison', 'indistinguishable_atoms')
    config.comp_brute_neighbors = parser.getboolean('Structure Comparison', 'brute_neighbors')
    config.comp_neighbor_cutoff = parser.getfloat('Structure Comparison', 'neighbor_cutoff')
    config.comp_use_covalent = parser.getboolean('Structure Comparison', 'use_covalent')
    config.comp_covalent_scale = parser.getfloat('Structure Comparison', 'covalent_scale')

    #AKMC options
    config.akmc_confidence              = parser.getfloat('AKMC', 'confidence')
    config.akmc_thermal_window          = parser.getfloat('AKMC', 'thermally_accessible_window')
    config.akmc_max_thermal_window      = parser.getfloat('AKMC', 'thermally_accessible_buffer')
    config.akmc_max_kmc_steps           = parser.getint('AKMC', 'max_kmc_steps')
    config.akmc_confidence_scheme       = parser.get('AKMC', 'confidence_scheme')
    config.akmc_confidence_correction   = parser.getboolean('AKMC', "confidence_correction")

    #Basin Hopping options
    config.bh_md_probability = parser.getfloat('Basin Hopping', 'md_probability')
    
    #path options
    config.path_root         = parser.get('Paths', 'main_directory')
    config.path_jobs_out     = parser.get('Paths', 'jobs_out')
    config.path_jobs_in      = parser.get('Paths', 'jobs_in')
    config.path_incomplete   = parser.get('Paths', 'incomplete')
    config.path_states       = parser.get('Paths', 'states')
    config.path_results      = parser.get('Paths', 'results')
    config.path_pot          = parser.get('Paths', 'potential_files')
    config.path_bh_minima    = parser.get('Paths', 'bh_minima')

    #Rye-requested check
    if not gave_config and not os.path.samefile(config.path_root, os.getcwd()):
        res = raw_input("The config.ini file in the current directory does not point to the current directory. Are you sure you want to continue? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
            pass
        else:
            sys.exit(3)
        
    #communicator options
    config.comm_type = parser.get('Communicator', 'type')
    #print comm_type
    config.comm_job_bundle_size = parser.getint('Communicator', 'jobs_per_bundle')
    config.comm_job_buffer_size = parser.getint('Communicator', 'num_jobs')
    config.path_scratch = parser.get('Paths', 'scratch')
    #print path_scratch
    if config.comm_type == 'local':
        config.comm_local_client = parser.get('Communicator', 'client_path')
        config.comm_local_client = parser.get('Communicator', 'client_path')
        config.comm_local_ncpus = parser.getint('Communicator', 'number_of_CPUs')
    if config.comm_type == 'cluster':
        config.comm_script_path = parser.get('Communicator', 'script_path')
        config.comm_script_name_prefix = parser.get('Communicator', 'name_prefix')
        config.comm_script_queued_jobs_cmd = parser.get('Communicator', 'queued_jobs')
        config.comm_script_cancel_job_cmd = parser.get('Communicator', 'cancel_job')
        config.comm_script_submit_job_cmd = parser.get('Communicator', 'submit_job')
    if config.comm_type == 'mpi':
        pass
    if config.comm_type == 'boinc':
        config.comm_boinc_project_dir = parser.get('Communicator', 'boinc_project_dir')
        config.comm_boinc_wu_template_path = parser.get('Communicator', 'boinc_wu_template_path')
        config.comm_boinc_re_template_path = parser.get('Communicator', 'boinc_re_template_path')
        config.comm_boinc_appname = parser.get('Communicator', 'boinc_appname')
        config.comm_boinc_results_path = parser.get('Communicator', 'boinc_results_path')
        config.comm_boinc_priority = parser.getint('Communicator', 'boinc_priority')
        #print comm_boinc_project_dir
        #print comm_boinc_wu_template_path
        #print comm_boinc_re_template_path
        #print comm_boinc_appname
        #print comm_boinc_results_path
    if config.comm_type == 'arc':
        if parser.has_option('Communicator', 'client_path'):
            config.comm_client_path = parser.get('Communicator', 'client_path')
        else:
            config.comm_client_path = ""

        if parser.has_option('Communicator', 'blacklist'):
            config.comm_blacklist = [ string.strip(c) for c in parser.get('Communicator', 'blacklist').split(',') ]
        else:
            config.comm_blacklist = []

    #Saddle Search options
    config.disp_type = parser.get('Saddle Search', 'displace_type')
    if config.disp_type == 'water':
        config.stdev_translation = parser.getfloat('Saddle Search', 'stdev_translation') # undocumented
        config.stdev_rotation = parser.getfloat('Saddle Search', 'stdev_rotation') # undocumented
        config.molecule_list = eval(parser.get('Saddle Search', 'molecule_list')) # undocumented
        config.disp_at_random = parser.getint('Saddle Search', 'disp_at_random') # undocumented
    else:
        config.disp_magnitude= parser.getfloat('Saddle Search', 'displace_magnitude')
        config.disp_radius = parser.getfloat('Saddle Search', 'displace_radius')
        config.disp_min_norm = parser.getfloat('Saddle Search', 'displace_min_norm')
    if config.disp_type == 'under_coordinated':
        config.disp_max_coord = parser.getint('Saddle Search', 'displace_max_coordination')
    if config.disp_type == 'listed_atoms':
        config.disp_listed_atoms = [ int(string.strip(c)) for c in parser.get('Saddle Search', 'displace_atomlist').split(',') ]

    #KDB
    config.kdb_on = parser.getboolean('KDB', 'use_kdb')
    config.kdb_only = parser.getboolean('KDB', 'kdb_only')
    config.kdb_scratch_path = parser.get('Paths', 'kdb_scratch')
    config.kdb_path = parser.get('Paths', 'kdb')
    config.kdb_wait = parser.get('KDB', 'wait')
    config.kdb_addpath = parser.get('KDB', 'addpath')
    if config.kdb_addpath == "False":
        config.kdb_addpath = os.path.join(os.path.dirname(__file__), "kdb", "kdbinsert.py")
    config.kdb_querypath = parser.get('KDB', 'querypath')
    if config.kdb_querypath == "False":
        config.kdb_querypath = os.path.join(os.path.dirname(__file__), "kdb", "kdbquery.py")

    #Recycling
    config.recycling_on = parser.getboolean('Recycling', 'use_recycling')
    config.recycling_save_sugg = parser.getboolean('Recycling', 'save_suggestions')
    if not config.recycling_on:
        config.disp_moved_only = False
    else:
        config.disp_moved_only = parser.getboolean('Recycling', 'displace_moved_only')
    config.recycling_move_distance = parser.getfloat('Recycling', 'move_distance')
    config.sb_recycling_on = parser.getboolean('Recycling','use_sb_recycling')
    config.sb_recycling_path = None
    if config.sb_recycling_on:
        config.sb_recycling_path = parser.get('Paths', 'superbasin_recycling')

    #Coarse Graining
    config.sb_on = parser.getboolean('Coarse Graining', 'use_projective_dynamics')
    config.sb_state_file = parser.get('Coarse Graining', 'state_file') 
    config.sb_path = None
    if config.sb_on:
        config.sb_path = parser.get('Paths', 'superbasins')
        config.sb_scheme = parser.get('Coarse Graining', 'superbasin_scheme')
        if config.sb_scheme == 'transition_counting':
            config.sb_tc_ntrans = parser.getint('Coarse Graining', 'number_of_transitions')
        elif config.sb_scheme == 'energy_level':
            config.sb_el_energy_increment = parser.getfloat('Coarse Graining', 'energy_increment')

    config.askmc_on = parser.getboolean('Coarse Graining','use_askmc')
    if config.askmc_on:
        config.askmc_confidence = parser.getfloat('Coarse Graining','askmc_confidence')
        config.askmc_alpha = parser.getfloat('Coarse Graining','askmc_barrier_raise_param')
        config.askmc_gamma = parser.getfloat('Coarse Graining','askmc_high_barrier_def')
        config.askmc_barrier_test_on = parser.getboolean('Coarse Graining','askmc_barrier_test_on')
        config.askmc_connections_test_on = parser.getboolean('Coarse Graining','askmc_connections_test_on')

    #Debug options
    config.debug_interactive_shell = parser.getboolean('Debug', 'interactive_shell')
    if config.debug_interactive_shell:
        import signal, code
        signal.signal(signal.SIGQUIT, lambda signum, frame: code.interact(local=locals()))
    config.debug_keep_bad_saddles  = parser.getboolean('Debug', 'keep_bad_saddles')
    config.debug_keep_all_results  = parser.getboolean('Debug', 'keep_all_result_files')
    config.debug_results_path = parser.get('Debug', 'result_files_path')
    config.debug_register_extra_results = parser.getboolean('Debug', 'register_extra_results')
    config.debug_use_mean_time = parser.getboolean('Debug', 'use_mean_time')
    config.debug_target_trajectory = parser.get('Debug', 'target_trajectory')

    del parser
