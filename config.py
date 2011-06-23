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
    def __init__(self, name, kind, description):
        self.name = name
        self.description = description
        self.kind = kind
        self.values = []

class ConfigValue:
    def __init__(self, name, description):
        self.name = name
        self.description = description
        
# Main

tempSection = ConfigSection("Main", "These are the options that go in the 'Main' section of config.ini")

tempKey = ConfigKey("job", "string", "The type of job to execute.")
tempKey.values.append(ConfigValue("akmc", "Run an adaptive kinetic monte carlo simulation."))
tempKey.values.append(ConfigValue("parallel_replica", "Calculate the rare-event dynamics of the system by combining transitions observed from multiple trajectories run in parallel."))
tempKey.values.append(ConfigValue("process_search", "Combined saddle search, minimizations, and prefactor calculations. Used by the aKMC method."))
tempKey.values.append(ConfigValue("saddle_search", "Do a saddle point search using a minimum mode method."))
tempKey.values.append(ConfigValue("minimization", "Find the minimum from an initial configuration."))
tempKey.values.append(ConfigValue("hessian", "Calculate the Hessian matrix for the specified configuration in a process."))
tempKey.values.append(ConfigValue("dimer_dr", "Rye is changing this."))
tempKey.values.append(ConfigValue("dimer_rotation", "Rye is changing this."))
tempKey.values.append(ConfigValue("displacement_sampling", "Job to sample different displacement methods and parameters to see which are the most efficient."))

tempSection.keys.append(tempKey)

tempKey = ConfigKey("temperature", "float", "The temperature that the job will run at.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("random_seed", "int", "Takes an integer number for the random seed. If this number is less than zero the current time is used as the random seed.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("potential", "string", "the type of potential to execute")
tempKey.values.append(ConfigValue("lj", "Lennard-Jones potential in reduced units."))
tempKey.values.append(ConfigValue("morse_pt", "Morse potential for platinum."))
tempKey.values.append(ConfigValue("emt", "Effective medium theory, for metals."))
tempKey.values.append(ConfigValue("edip", "Environment-Dependent Interatomic Potential, for carbon."))
tempKey.values.append(ConfigValue("vasp", "Vienna Ab-Initio Simulation Program (VASP) interface."))
tempKey.values.append(ConfigValue("tersoff_si", "Tersoff pair potential with angular terms, for silicon."))
tempKey.values.append(ConfigValue("sw_si", "Stillinger-Weber potential, for silicon."))
tempKey.values.append(ConfigValue("lenosky_Si", "Lenosky potential, for silicon."))
tempKey.values.append(ConfigValue("eam_al", "Embedded atom method parameterized for aluminum."))
tempKey.values.append(ConfigValue("qsc", "Quantum Sutton-Chen potential, for FCC metals."))
tempKey.values.append(ConfigValue("zpice", "Water on platinum."))
tempKey.values.append(ConfigValue("tip4p", "Point charge model for water."))
tempKey.values.append(ConfigValue("bopfox", "Bond order potential, for metals."))

tempSection.keys.append(tempKey)

config.format.append(tempSection)

# AKMC

tempSection = ConfigSection("AKMC", "Parameters for the AKMC section of config.ini.")

tempKey = ConfigKey("confidence", "float", "The confidence (out of 1.0) criterion for moving to the next state.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("max_kmc_steps", "int", "The maximum number of transitions per execution of the server.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("thermally_acessible_window", "float", "Processes with barriers within this number of kT above the lowest barrier will be used in the rate table and for confidence calculations.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("thermally_accessible_buffer", "float", "Processes with barriers of thermally_accessible_window + thermally_accessible_buffer will be stored, in the event that they are thermally accessible later, but are not used in the rate table or for the confidence calculations. Processes with barriers higher than the sum of these two values will be discarded.")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

# Structure Comparison

tempSection = ConfigSection("Structure Comparison", "") #####edit descrition#####

tempKey = ConfigKey("energy_difference", "float", "How close in energy two configurations must be to be considered energetically equivalent.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("distance_difference", "float", "The maximum distance two mapped atoms may be for two configurations to be considered equivalent.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("indistinguishable_atoms", "boolean", "Use an algorithm to compare structures that does not distinguish between atoms of the same element. That is to say the numbering of the atoms does not affect the structural comparison.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("check_rotation", "boolean", "Finds optimal overlap of structures via rotation before comparing them. Use this option in systems where structures can become rotated, such as nanoparticles.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("neighbor_cutoff", "float", "Atoms within this distance of each other are considered neighbors.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("use_covalent", "boolean", "Use the covalent radii of atoms to determine neighbors.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("covalent_scale", "float", "Multiply covalent radii by this amount before determining neighbors.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("brute_neighbors", "boolean", "Determine neighbors by brute force (use this with nonorthogonal boxes).")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

# PathsMaximum number of rotations allowed for the dimer in each step.Improvements to the dimer method from Kastner.Dimer is considered converged if it rotates fewer degrees than this.

tempSection = ConfigSection("Paths", "Location of files related to the sending and receiving data on the server.")

tempKey = ConfigKey("main_directory", "string" , "This is the root directory of the simulation. Configuration files and the initial reactant are here and by default all of the simulation data will be stored under this directory.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("searches_in", "string", "") ####edit description later####

tempSection.keys.append(tempKey)

tempKey = ConfigKey("states", "string", "Where all of the information about individual states is located.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("scratch", "string", "") ####edit description later####

tempSection.keys.append(tempKey)

tempKey = ConfigKey("potential_files", "string", "For extra files needed by the client for the potential.")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

#Process Search

tempSection=ConfigSection("Process Search", "The AKMC method can ask clients to do a saddle search, find connecting minima, and calculate prefactors all within the Process Search job type.")

tempKey = ConfigKey("minimize_first", "boolean", "Every time a process search is run by a client the reactant will be minimized first before doing any saddle searches.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("prefactor_min", "float", "Minimum value for a resonable prefactor.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("prefactor_max", "float", "Maximum value for a resonable prefactor.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("default_prefactor", "int", "Calculate prefactor if zero, otherwise use given value instead of doing a full prefactor calculation.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("minimization_offset", "float", "After a saddle is found, images are placed on either side of the saddle along the mode and minimized to ensure that the saddle is connected to the original minimum and to locate the product state. This is the distance those images are displaced from the saddle.")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

#Saddle Search

tempSection = ConfigSection("Saddle Search", "A saddle search is initiated by making a local displacement of atoms from their position at the minimum of the current state. This displacement can be done using the different strategies indicated by the displace_type option, and the following parameters. If the user knows something about the local environment where reactions are likely to take place in the system, this information can be used to make saddle searches more efficient by getting them started in the right part of configuration space.")

tempKey = ConfigKey("displace_type", "string", "Type of displace to use")
tempKey.values.append(ConfigValue("random", "Select an atom at random from the free atoms in the configuration."))
tempKey.values.append(ConfigValue("least_coordinated", "Determine the lowest coordination number of all atoms in the configuration and select one atom at random with that coordination number."))
tempKey.values.append(ConfigValue("under_coordinated", "Select a random atom with coordination less than displace_max_coordination."))
tempKey.values.append(ConfigValue("listed_atoms", "Select an atom from the list passed as displace_atomlist"))
tempKey.values.append(ConfigValue("not_FCC_HCP", "Select an atom that is not in a local HCP or FCC coordination."))
tempKey.values.append(ConfigValue("client_least_coordinated", "Displacement is made on the client centered on an atom with the minimal coordination."))
tempKey.values.append(ConfigValue("cleint_not_FCC_HCP_coordinated", "Displacement is made on the client centered on an atom that is not in a local HCP or FCC coordination."))
tempKey.values.append(ConfigValue("client_last_atom", "Displacement is made on the client centered on the last atom in the configuration."))

tempSection.keys.append(tempKey)

tempKey = ConfigKey("max_iterations", "int", "The maximum number of translation steps to be taken.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("displace_radius", "float", "Atoms within this distance of the epicenter will be displaced.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("displace_magnitude", "float", "The standard deviation of the magnitude of the displacement.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("displace_min_norm", "float" , "The total length of the displacement vector is ensured to exceeds this value. Is useful when only a few degrees of freedoms are displaced to guarantee that the starting point for the saddle point search is significantly different from the initial minimum.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("displace_max_coordination", "int", "When using under_coordinated as the displacement type, choose only atoms with a coordination equal to or less than this.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("min_mode_method", "string", "") ####edit description later####
tempKey.values.append(ConfigValue("dimer", "Use the dimer min-mode method."))
tempKey.values.append(ConfigValue("lanczos", "Use the Lanczos min-mode method."))

tempSection.keys.append(tempKey)

tempKey = ConfigKey("max_energy", "float", "The energy at which a saddle search is considered bad and terminated.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("displace_atomlist","string", "The individual index should be seperated by a comma 10, 20,-1 would be the 10, 20 and the last atom.") ####possibily not a string####

tempSection.keys.append(tempKey)

tempKey = ConfigKey("client_max_single_displace","int", "Only functional when displacement is done on the client. Defines the maximal allowed absolute value for a single component in the displacement vector.")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

# Dimer

tempSection = ConfigSection("Dimer", "Options for controlling the dimer minimum mode finding method on the client side.")

tempKey = ConfigKey("seperation", "float", "Seperation between dimer images.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("finite_angles", "float", "Finite difference angle over which the dimer is rotated to find the lowest curvature.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("rotations_min", "int", "Minimum number of rotations allowed for the dimer in each step.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("rotations_max", "int", "Maximum number of rotations allowed for the dimer in each step.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("torque_min", "float", "Minimum torque above which the dimer rotates only once and below which is does not rotate.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("torque_max", "float", "Maximum torque above which the dimer rotates up to rotations_max times.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("improved", "boolean", "Improvements to the dimer method from Kastner.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("coverged_rotatoin", "float", "Dimer is considered converged if it rotates fewer degrees than this.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("opt_method", "string", "Optimization algorithm to choose the dimer rotation direction")
tempKey.values.append(ConfigValue("sd", "steepest descent, rotate along the rotational force."))
tempKey.values.append(ConfigValue("cg", "conjudate gradient, rotate along conjugate directions."))
tempKey.values.append(ConfigValue("lbfgs", "quasi-Newton method [not implemeted yet].")) ####not implemented yet####

tempSection.keys.append(tempKey)

config.format.append(tempSection)

# Lanczos

tempSection = ConfigSection("Lanczos", "Options for controlling the Lanczos minimum mode finding method on the client side.")

tempKey = ConfigKey("max_iterations", "int", "The maximum number of refinement iterations when calculating the minimum eigenvalue.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("tolerance", "float", "This is the convergence critera for relative error of the lowest eigenvalue.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("finite_dist", "float", "Finite difference step size when computing the second derivative of the potential along the Lanczos vectors.")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

# Hessian

tempSection = ConfigSection("Hessian", "Options for controlling how Hessian matricies are calculated.")

tempKey = ConfigKey("type", "string", "The Hessian to be calculated has to be one of reactant, saddle, or product.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("finite_dist", "float", "Finite difference distance between forces used to construct the Hessian.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("min_displacement", "float", "Minimum amount that an atom has to move to be included in the Hessian calculation.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("within_radius", "float", "Atoms within this radius of moving atoms are included in the Hessian.")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

# Nudged Elastic Band

tempSection = ConfigSection("Nudged Elastic Band", "Options for controlling the dimer minimum mode finding method on the client side.")

tempKey = ConfigKey("images", "int", "Number of NEB images between the fixed endpoints.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("spring", "float", "Spring constant between images, in eV/Ang.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("climbing_image_method", "boolean", "Use the climbing image method to move the highest energy image to the saddle.")

tempSection.keys.append(tempKey)

# Communicator

tempSection = ConfigSection("Communicator", "Options that apply to all of the different communicator types.")

tempKey = ConfigKey("type", "string", "") #no description
tempKey.values.append(ConfigValue("local", "The local communicator runs the calculations on the same computer that the server is run on."))
tempKey.values.append(ConfigValue("cluster", "A job scheduler can be used to run jobs through user supplied shell scripts. Examples are given for SGE."))
tempKey.values.append(ConfigValue("boinc", "Jobs can be submitted to a BOINC project."))
tempKey.values.append(ConfigValue("arc", "Jobs can be submitted to the grid computing software ARC."))

tempSection.keys.append(tempKey)

tempKey = ConfigKey("num_jobs", "int", "Local( The number of jobs that will be run every time the program is invoked) Cluster( The desired sum of the queued and running jobs.) Boinc( The number of jobs to keep in the queue.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("jobs_per_bundle", "int", "In eon a job is defined as task that the eon client executes, such as a process search or a parallel replica run. Sometimes it makes sense to run more than one of the same type of job at a time.")

tempSection.keys.append(tempKey)
    # Local options
tempKey = ConfigKey("client_path", "string", "Either the name or path to the eon client binary. If only a name and not a path is given then eon looks for the binary in same directory as config.ini failing to find it there it will search though the directories in the $PATH environment variable.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("number_of_cpus", "int", "The number of jobs that will run simultaneously.")

tempSection.keys.append(tempKey)
    #Cluster options
tempKey = ConfigKey("script_path", "string", "The path to the user defined scripts for submitting jobs to the communicator.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("name_prefix", "string", "When jobs are submitted to the scheduler they are given a unique internally used named. In order to make the jobs identifiable by the user the name_prefix can be set to a meaningful string that will always be prepended to the job names.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("queued_jobs", "string", "This is the name of the script that returns the job ids of all the running and queued jobs. It does not have to return the job ids of only eon related jobs.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("submit_job", "string", "This is the name of the script that submits a single job to the queuing system. It takes two command line arguments. The first is the name of the job. This is not required for eon use, but is highly recommended so that users can identify which job is which. The second argument is the working directory. This is the path where the eon client should be executed. All of the needed client files will be placed in this directory. The script must return the job id of the submitted job. This is how eon internally keeps track of jobs.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("cancel_job", "string", "This is the name of the script that cancels a job. It takes a single argument the job id.")

tempSection.keys.append(tempKey)
    #Boinc options
tempKey = ConfigKey("boinc_project_dir", "string", "This is the full path to the root of the BOINC project directory.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("boinc_wu_template_path", "string", "This is the path, relative from the boinc_project_dir, to the boinc workunit template.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("boinc_re_template_path", "string", "This is the path, relative from the boinc_project_dir, to the boinc result template.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("boinc_appname", "string", "This is the name of the application in BOINC.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("boinc_results_path", "string", "This is the path where BOINC puts the final results. If you are using the sample_assimilator the results are stored in the project directory in a folder named sample_results.")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

#Parallel Replica

tempSection = ConfigSection("Parallel Replica", "Parallel Replica Dynamics (PRD) is the simplest and the most accurate way to do accelerated-MD simulation. The only assumption made in this method is that the reactions satisfy first order kinetics.")
    #Dynamics options
tempKey = ConfigKey("time_step", "float", "The length of each MD step in femtoseconds.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("auto_stop", "boolean", "Whether or not stop the job when a new state is found. For boinc communicator this value should be set to false.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("steps", "int", "The number of MD steps to run.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("dephase_steps", "int", "Number of steps used to decorrelate the replica trajectories. The momenta will be inversed when reaching the dividing surface to prevent transitions occurring during this period.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("check_period", "int" , "How frequently the state of system is checked. Every check_period steps, the current structure and the initial one will be compared to tell whether a newstate has been reached. Also note when you set refine as true, the code will keep a buffer array consisting of check_period/record_resolution+1 atomic configurations, which may increase the usage of memory.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("refine_transition_time", "boolean", "Whether or not the transition time is refined. When this option is turned on, the code will keep an array consisted by check_period/record_resolution+1 atomic configurations. A Binary search algorithm is employed to determine the transition step. Otherwise the transition step would be the first in which a new state was found. This function reduces the need for a smaller check_period. And the accuracy of transition time is record_resolution*timestep.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("record_resolution", "int", "How often the system is recorded to the buffer array when the refine_transition_time option is activated. Increasing the value of record_resolution lowers the accuracy of the transition time estimate but also reduces memory usage and speeds up refinement of the transition step.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("post_transition_steps", "int", "Number of MD steps which will be performed after a new state has been found. A state check will be employed after these post_transition_steps to confirm that the state is stable. This additional check helps avoid meta-stable states. A value similar to dephase_steps is recommended.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("thermo_type", "string", "") #no description
tempKey.values.append(ConfigValue("andersen", "Andersen thermostat with Verlet algorithm"))
tempKey.values.append(ConfigValue("nose_hoober", "Nose-Hover thermostat with Verlet algorithm"))
tempKey.values.append(ConfigValue("langevin", "Langevin thermostat with Verlet algorithm"))

tempSection.keys.append(tempKey)

tempKey = ConfigKey("andersen_alpha", "float", "The collision strength in the Andersen thermostat")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("andersen_collision_steps", "float", "The collision period (in MD steps) for the Andersen thermostat.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("nose_mass", "float", "The effective mass of the additional degree of freedom in the Nose-Hover thermostat, which determines the rate of heat transfer.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("nose_mass", "float", "The effective mass of the additional degree of freedom in the Nose-Hoover thermostat, which determines the rate of heat transfer.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("langevin_friction", "float", "The damping coefficient for langevin dynamics.")

tempSection.keys.append(tempKey)
    #Hyperdynamics options
tempKey = ConfigKey("bias_potential", "string", "") #no description
tempKey.values.append(ConfigValue("none", "with no bias potential, run regular MD"))
tempKey.values.append(ConfigValue("bond_boost","bond boost method from Miron and Fichthorn"))

tempSection.keys.append(tempKey)

tempKey = ConfigKey("bb_dvmax", "float", "The magnitude of the bond-boost bias potential. It should be smaller than the barrier of any transition.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("bb_rmd_steps", "int", "Number of MD steps used to determine the equilibrium bond length before the bias potential is added.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("bb_stretch_threshold", "float", "Defines the bond-boost dividing surface. It should be smaller than the maximum fractional nearest-neighbor bond stretch or compression at any transition state.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("bb_ds_curvature", "float", "The curvature near the bond-boost dividing surface, it should has a value <= 1. We recommend the value to be 0.9-0.98.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("bb_rcut", "float", "All bonds which belong to the tagged atoms and are shorter than a cutoff of rcut will be included in the bond-boost potential.")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

# Basin Hopping

tempSection = ConfigSection("Basin Hopping", "Basin hopping is a Monte Carlo method in which the energy of each configuration is taken to be the energy of a local minimization.")

tempKey = ConfigKey("step_size", "float", "Mean displacement in each degree of freedom")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("steps", "int", "Total number of steps.")

tempSection.keys.append(tempKey)

#config.format.append(tempSection)

# Optimizers

tempsection = ConfigSection("Optimizers", "The eon client provides two optimizers quick-min and conjugate gradients. These options are in the [Optimizers] section.")

tempKey = ConfigKey("opt_method", "string", "The optimization method to use.")
tempKey.values.append(ConfigValue("cg", "Conjugate gradient"))
tempKey.values.append(ConfigValue("qm", "Quickmin"))
tempKey.values.append(ConfigValue("box", "Optimizes the atom positions and box using quickmin"))

tempSection.keys.append(tempKey)

tempKey = ConfigKey("max_iterations", "int", "The maximum number of optimization iterations that the will be performed.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("converged_force", "float", "When the maximum force (in eV/A) on any one atom is smaller than this value, the structure is considered minimized.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("max_move", "float", "Maximum distance that an atom may be moved in a single optimization step.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("time_step", "float", "The dynamical timestep for the quickmin algorithm.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("finite_dist", "float", "The finite difference step size.")

config.format.append(tempSection)

# Coarse Graining

tempSection = ConfigSection("Coarse Graining", "In AKMC simulations where there are vastly different rates, the simulation can get stuck in a group of states connected by relatively fast rates. In order to explore slower transitions, a prohibitively large number of KMC steps may be needed. In order to circumvent this problem, eOn implements two methods. The first method, projective dynamics [1], groups states that are joined by fast rates into superbasins. Information about transitions between states in a superbasin is lost, but the rates for transitions across a superbasin are correct. The second method, accelerated superbasining kinetic Monte Carlo (AS-KMC) [2], artificially raises low barriers. The dynamics between states connected by fast rates are simulated, but an error is introduced in the dynamics direction and time. Both methods cannot be used simultaneously.")

tempKey = ConfigKey("state_file", "string", "File name for the state specific data stored within each of the state directories.")

tempSection.keys.append(tempKey)
    #projective dynamics
tempKey = ConfigKey("use_projective_dynamics", "boolean", "This option determines whether the projective dynamics coarse graining method will be used. This mutually excludes the use_askmc option.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("superbasin_scheme", "string", "Projective dynamics provides a method for calculating transition rates across superbasins. An additional method is needed in order to decide when to combine states into a superbasin. eOn provides two methods. The first method, called transition counting, counts the number of times that the simulation has transitioned between a given pair of states. After a critical number of transitions have occured, the pair of states are merged to form a superbasin. (If one is already in a superbasin, the other is added to that superbasin. If both are already in superbasins, the two superbasins are merged). This method can be selected by setting scheme equal to transition_counting. Jean-Claude: document your scheme here This method can be elected by setting scheme equal to energy_level.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("number_of_transitions", "int", "If the transition counting scheme is being used (scheme=transition_counting), this is the number of transitions that must occur between two states before they are merged into a superbasin.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("energy_increment", "float", "If the energy level scheme is being used (scheme=energy_level). Each state, the first time it is visited, is assigned an energy level first equal to the energy of the minimum. Every time the state is visited again by the Monte Carlo simulation, the energy level is increased by this amount")

tempSection.keys.append(tempKey)
    # accelerated superbasin kinetic monte carlo
tempKey = ConfigKey("use_askmc", "boolean", "This option determines whether the AS-KMC coarse graining method will be used. This mutually excludes the use_projective_dynamics option.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("askmc_condifence", "float", "The confidence for AS-KMC. This value determines the accuracy of the direction of the dynamics trajectory.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("askmc_barrier_raise_param", "float", "This parameter sets how much the barriers are raised during AS-KMC. ( in the reference.)")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("askmc_high_barrier_def", "int", "This parameter sets how high a barrier must be to be considered high in AS-KMC.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("askmc_barrier_test_on", "boolean", "") # no description

tempSection.keys.append(tempKey)

tempKey = ConfigKey("askmc_connections_test_on", "boolean", "This parameter determines whether to ensure that there are no processes which connect states in the defined superbasin which have not been visited yet and which have a low-barrier. This check is somewhat more computationally expensive than the previous because structure comparisons have to be made when finding product states of unvisited processes.")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

# Kdb

tempSection = ConfigSection("Kdb", "One of the bottlenecks in an aKMC simulation is performing the saddle point searches. The kinetic database is used to ameliorate this cost by storing information about processes as they are found and using it to predict future saddle points.")

tempKey = ConfigKey("use_kdb", "boolean", "Turn Kdb on/off.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("wait", "boolean", "Wait for the query to finish before submitting jobs (for debugging purposes).")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("keep", "boolean", "Keep the saddle suggestions (for debugging purposes).")

tempSection.keys.append(tempKey)

config.format.append(tempSection)

# Debug

tempSection = ConfigSection ("Debug", "Parameters that are generally used to help debug calculations")

tempKey = ConfigKey("keep_bad_saddles", "boolean", "Keep data about bad saddles. If true, the result files for failed saddle searches are kept in the badprocdata directory within the state directory for that search.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("keep_all_results_files", "boolean", "Stores all result files in main_directory/results")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("register_extra_results", "boolean", "Register processes found for a state after leaving that state.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("use_mean_time", "boolean", "Select transition times from the mean of the exponential distribution of escape times.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("target_trajectory", "boolean", "Follow the state-to-state trajectory of another akmc simulation.")

tempSection.keys.append(tempKey)

tempKey = ConfigKey("write_movies", "boolean", "Causes the client to output movies of minimizations and saddle searches.")

tempKey = ConfigKey("save_stdout", "boolean", "Save the standard output from the client to a file named stdout_0.dat")

tempSection.keys.append(tempKey)

config.format.append(tempSection)



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
    config.comp_check_rotation = parser.getboolean('Structure Comparison', 'check_rotation')
    config.comp_brute_neighbors = parser.getboolean('Structure Comparison', 'brute_neighbors')
    config.comp_neighbor_cutoff = parser.getfloat('Structure Comparison', 'neighbor_cutoff')
    config.comp_use_covalent = parser.getboolean('Structure Comparison', 'use_covalent')
    config.comp_covalent_scale = parser.getfloat('Structure Comparison', 'covalent_scale')

    #aKMC options
    config.akmc_confidence              = parser.getfloat('AKMC', 'confidence')
    config.akmc_thermal_window          = parser.getfloat('AKMC', 'thermally_accessible_window')
    config.akmc_max_thermal_window      = parser.getfloat('AKMC', 'thermally_accessible_buffer')
    config.akmc_max_kmc_steps           = parser.getint('AKMC', 'max_kmc_steps')
    config.akmc_confidence_scheme       = parser.get('AKMC', 'confidence_scheme')
    config.akmc_confidence_correction   = parser.getboolean('AKMC', "confidence_correction")
    
    #path options
    config.path_root         = parser.get('Paths', 'main_directory')
    config.path_jobs_out     = parser.get('Paths', 'jobs_out')
    config.path_jobs_in      = parser.get('Paths', 'jobs_in')
    config.path_incomplete   = parser.get('Paths', 'incomplete')
    config.path_states       = parser.get('Paths', 'states')
    config.path_results      = parser.get('Paths', 'results')
    config.path_pot          = parser.get('Paths', 'potential_files')
    

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
    if config.sb_recycling_on:
        config.sb_recycling_path = parser.get('Paths', 'superbasin_recycling')

    #Coarse Graining
    config.sb_on = parser.getboolean('Coarse Graining', 'use_projective_dynamics')
    config.sb_state_file = parser.get('Coarse Graining', 'state_file') 
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
    config.debug_register_extra_results = parser.getboolean('Debug', 'register_extra_results')
    config.debug_use_mean_time = parser.getboolean('Debug', 'use_mean_time')
    config.debug_target_trajectory = parser.get('Debug', 'target_trajectory')

    del parser
