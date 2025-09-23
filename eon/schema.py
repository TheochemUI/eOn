from pydantic import (
    BaseModel,
    Field,
    validator,
    model_validator,
    ConfigDict,
)
from typing import Optional, Any, Union
from typing_extensions import Literal
import math
import random


class MainConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    job: Literal[
        "akmc",
        "basin_hopping",
        "displacement_sampling",
        "dynamics",
        "escaperate",  # TODO(rg): Document
        "finite_differences",
        "global_optimization",
        "gp_surrogate",
        "hessian",
        "minimization",
        "monte_carlo",
        "nudged_elastic_band",
        "parallel_replica",  # Alias for unbiased_parallel_replica
        "point",
        "prefactor",
        "process_search",
        "replica_exchange",
        "saddle_search",
        "safe_hyperdynamics",
        "structure_comparison",
        "tad",
    ] = Field(
        default="akmc",
        description="The type of job to execute.",
    )
    """
    Options:
     - ``molecular_dynamics``: Molecular dynamics.
     - ``parallel_replica``: Calculate the rare-event dynamics of the system by combining transitions observed from multiple trajectories run in parallel.
     - ``saddle_search``: Do a saddle point search using a minimum mode method.
     - ``displacement_sampling``: Job to sample different displacement methods and parameters to see which are the most efficient.
     - ``process_search``: Combined saddle search, minimizations, and prefactor calculations. Used by the aKMC method.
     - ``basin_hopping``: Search for global minimum using basin hopping method.
     - ``minimization``: Find the minimum from an initial configuration.
     - ``akmc``: Run an adaptive kinetic monte carlo simulation.
     - ``hessian``: Calculate the Hessian matrix for the specified configuration in a process.
     - ``finite_differences``: Finite difference methods.
     - ``global_optimization``: Global optimization methods.
     - ``gp_surrogate``: Gaussian process surrogate methods.
     - ``monte_carlo``: Monte Carlo methods.
     - ``nudged_elastic_band``: Nudged elastic band methods.
     - ``point``: Single point calculation.
     - ``prefactor``: Prefactor calculations.
     - ``replica_exchange``: Replica exchange methods.
     - ``safe_hyperdynamics``: Safe hyperdynamics methods.
     - ``structure_comparison``: Structure comparison methods.
     - ``tad``: Temperature-accelerated dynamics.
    """
    random_seed: int = Field(
        # See https://click.rgoswami.me/pyrandint
        # for a few caveats on getrandbits
        default_factory=lambda: random.getrandbits(32),
        description="Takes an integer for the random seed. If this number is less than zero the current time is used as the random seed. If it is not defined, a 32 bit sized integer is used.",
    )
    temperature: float = Field(
        default=300.0, description="The temperature that the job will run at."
    )
    finite_difference: float = Field(
        default=0.01,
        description="The finite difference distance to use for dimer, hessian, lanczos, and optimization methods.",
    )
    checkpoint: bool = Field(
        default=False, description="Resuming search from checkpoint files."
    )
    quiet: bool = Field(default=False, description="Disable logging to stdout.")
    write_log: bool = Field(
        default=True, description="Enable writing log to client.log file."
    )
    max_force_calls: int = Field(
        default=0,
        description="The maximum number of total force calls per job. The default, 0, means unlimited force calls. If this limit is reached, an error code 1017 is thrown and shows up in the client log.",
    )
    remove_net_force: bool = Field(
        default=True,
        description="If True, ensures that the net force on the system of atoms is zero by adjusting the force on each free atom.",
    )


class StructureComparisonConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    energy_difference: float = Field(
        default=0.01,
        description="How different in energy two configurations must be to be considered different structures.",
    )
    distance_difference: float = Field(
        default=0.1,
        description="The maximum distance two mapped atoms may be for two configurations to be considered equivalent.",
    )
    indistinguishable_atoms: bool = Field(
        default=True,
        description="Use an algorithm to compare structures that does not distinguish between atoms of the same element. The numbering of the atoms does not affect the structural comparison.",
    )
    check_rotation: bool = Field(
        default=False,
        description="Finds optimal overlap of structures via rotation before comparing them. Use this option in systems where structures can become rotated, such as nanoparticles.",
    )
    brute_neighbors: bool = Field(
        default=False,
        description="Determine neighbors by brute force (use this with nonorthogonal boxes).",
    )
    neighbor_cutoff: float = Field(
        default=3.3,
        description="Atoms within this distance of each other are considered neighbors.",
    )
    use_covalent: bool = Field(
        default=False,
        description="Use the covalent radii of atoms to determine neighbors.",
    )
    covalent_scale: float = Field(
        default=1.3,
        description="Scale factor for covalent radii when determining neighbors.",
    )
    remove_translation: bool = Field(
        default=True,
        description="Remove translational components when comparing structures.",
    )


class AKMCConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    confidence: float = Field(
        default=0.99,
        description="The confidence (out of 1.0) criterion for moving to the next state.",
    )
    server_side_process_search: bool = Field(
        default=False,
        description='If true, the server does not send out "process_search" jobs. Instead, it manages individual dimer, minimization, and Hessian calculations. This option is usually used when a single "process_search" job will take a long time (hours or longer).',
    )
    thermally_accessible_window: float = Field(
        default=20.0,
        description="Processes with barriers within this number of kT above the lowest barrier will be used in the rate table and for confidence calculations.",
    )
    thermally_accessible_buffer: float = Field(
        default=0.0,
        description="Processes with barriers of thermally_accessible_window + thermally_accessible_buffer will be stored, in the event that they are thermally accessible later, but are not used in the rate table or for the confidence calculations. Processes with barriers higher than the sum of these two values will be discarded.",
    )
    max_kmc_steps: int = Field(
        default=0,
        description="The maximum number of KMC transitions in a row. In MPI or continuous mode, EON will exit after performing this many KMC steps. If this is set to 0, EON will run forever.",
    )
    confidence_scheme: Literal["old", "new", "sampling", "dynamics"] = Field(
        default="new", description="The scheme used for confidence calculation."
    )
    """
    Options:
    - ``old``: Traditional confidence scheme based on the number of unique processes (:math:`Nf`) and the number of searches (:math:`Ns`).
    - ``new``: Enhanced confidence scheme using a mathematical function to account for the ratio of :math:`Nf` to :math:`Ns`.
    - ``sampling``: Confidence based on sampling probabilities derived from the repeat counts of processes.
    - ``dynamics``: Confidence derived from molecular dynamics, considering time extrapolation and rates at different temperatures.
    """
    confidence_correction: bool = Field(
        default=False,
        description="If true, applies a correction to the confidence calculation.",
    )
    max_rate: float = Field(
        default=0.0,
        description="Sets the maximum allowable rate for a process. If a calculated rate exceeds this value, it is capped to prevent exceedingly high rates from dominating the simulation. This helps in maintaining numerical stability and realistic dynamics in the simulation.",
    )
    eq_rate: float = Field(
        default=0.0,
        description="Sets the equilibrium rate. If the forward and reverse rates exceed this value, they are adjusted to maintain equilibrium conditions. This prevents exceedingly high rates from dominating the simulation and ensures that both forward and reverse processes are considered symmetrically.",
    )


class BasinHoppingConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    steps: int = Field(
        default=10000,
        description="Number of steps to take at the assigned temperature.",
    )
    displacement: float = Field(
        default=0.5, description="Displacement in each degree of freedom."
    )
    write_unique: bool = Field(
        default=False,
        description="If true, the client will write out all the unique geometries and their energies",
    )
    """
    The files are named ``min_xxxxx.con`` and ``energy_xxxxx.con``, where
    'xxxxx' corresponds to the Monte Carlo step where the structure was first
    seen. For clusters it is important to set
    :any:`eon.schema.StructureComparisonConfig.check_rotation`.
    """
    initial_state_pool_size: int = Field(
        default=1,
        description="Number of lowest energy states in the pool from which the initial structures are drawn.",
    )
    """
    A value of 1 ensures simulation always starts from the lowest energy found so far.
    """
    initial_random_structure_probability: float = Field(
        default=0.0,
        description="Probability that the initial structure is generated randomly.",
    )
    push_apart_distance: float = Field(
        default=0.4,
        description="Push atoms apart until no atoms are closer than this distance.",
    )
    """
    This criterion is enforced for **both** the initial structure and all those
    generated by random displacement.
    """
    adjust_displacement: bool = Field(
        default=True,
        description="Flag to automatically adjust the displacement to meet the target acceptance ratio.",
    )
    adjust_fraction: float = Field(
        default=0.05,
        description="Fraction by which to change the step size in order to meet the target acceptance ratio.",
    )
    adjust_period: int = Field(
        default=10,
        description="Number of Monte Carlo steps between adjustments of the step size.",
    )
    target_ratio: float = Field(
        default=0.5,
        description="Target acceptance ratio used to determine whether to increase or decrease the step size.",
    )
    displacement_distribution: Literal["gaussian", "uniform"] = Field(
        default="gaussian",
        description="Distribution used for the displacement of each atom.",
    )
    """
    :any:`eon.schema.BasinHoppingConfig.displacement` for the **gaussian** will
    serve as the standard deviation of the normal distribution used to select
    displacements.

    **uniform** will select a random number between the positive and negative
    values of :any:`eon.schema.BasinHoppingConfig.displacement`.
    """
    swap_probability: float = Field(
        default=0.0,
        description="Probability (in range [0,1]) that a swapping step takes place instead of a displacement step.",
    )
    """
    The swap step selects two atoms of different elements and swaps their position.
    """
    single_atom_displace: bool = Field(
        default=False, description="Displace only one atom per step."
    )
    jump_max: int = Field(
        default=0,
        description="Number of consecutive rejected steps after which jump steps should be taken.",
    )
    """
    This provides a more global search when the structure is stuck in a certain
    basin. The number of jump steps is assigned in
    :any:`eon.schema.BasinHoppingConfig.jump_steps`. See
    :cite:t:`bh-iwamatsuBasinHoppingOccasional2004` for more details.
    """
    significant_structure: bool = Field(
        default=True, description="Displace from minimized structures."
    )
    """Implements the method of :cite:t:`bh-whiteInvestigationTwoApproaches1998`."""
    jump_steps: int = Field(
        default=0,
        description="Number of jump steps to take after jump_max consecutive rejections.",
    )


class PathsConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    main_directory: str = Field(
        default="./",
        description="This is the root directory of the simulation. Configuration files and the initial reactant are here and by default all of the simulation data will be stored under this directory.",
    )
    jobs_out: str = Field(default=None, description="Directory for job outputs.")
    jobs_in: str = Field(default=None, description="Directory for job inputs.")
    incomplete: str = Field(default=None, description="Directory for incomplete jobs.")
    states: str = Field(
        default=None,
        description="Directory where all of the information about individual states is located.",
    )
    results: str = Field(default=None, description="Directory for storing results.")
    potential_files: str = Field(
        default=None,
        description="Directory for extra files needed by the client for the potential.",
    )
    bh_minima: str = Field(default=None, description="Directory for BH minima.")
    kdb_scratch: str = Field(
        default=None, description="Directory for KDB scratch files."
    )
    kdb: str = Field(default=None, description="Directory for KDB files.")
    superbasins: str = Field(default=None, description="Directory for superbasins.")
    superbasin_recycling: str = Field(
        default=None, description="Directory for superbasin recycling."
    )
    scratch: str = Field(default=None, description="Directory for scratch files.")

    @model_validator(mode="before")
    def set_default_paths(cls, values: dict[str, Any]) -> dict[str, Any]:
        main_directory = values.get("main_directory", "./")
        defaults = {
            "jobs_out": f"{main_directory}/jobs/out/",
            "jobs_in": f"{main_directory}/jobs/in/",
            "incomplete": f"{main_directory}/jobs/incomplete/",
            "states": f"{main_directory}/states/",
            "results": main_directory,
            "potential_files": f"{main_directory}/potfiles",
            "bh_minima": f"{main_directory}/minima",
            "kdb_scratch": f"{main_directory}/kdbscratch/",
            "kdb": f"{main_directory}/kdb/",
            "superbasins": f"{main_directory}/superbasins/",
            "superbasin_recycling": f"{main_directory}/SB_recycling",
            "scratch": f"{main_directory}/jobs/scratch/",
        }
        for field, default in defaults.items():
            if values.get(field) is None:
                values[field] = default
        return values


class CommunicatorConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    type: Literal["local", "cluster", "mpi"] = Field(
        default="local",
        description="Communicator type",
    )
    """
    Options:
     - 'local': The local communicator runs the calculations on the same computer that the server is run on.
     - 'cluster': A job scheduler can be used to run jobs through user supplied shell scripts.
     - 'mpi': Allows for the server and clients to run as a MPI job.
    """
    jobs_per_bundle: int = Field(
        default=1,
        description="Number of jobs per bundle.",
    )
    """
    In ``eON``, a job is defined as a task that the ``eonclient`` executes, such
    as a process search or a parallel replica run. Sometimes it makes sense to
    run more than one of the same type of job at a time.

    For example, when using empirical potentials to do saddle searches a single
    search might only take several seconds on modern CPUs. In order to improve
    performance more than one client job (*e.g.*, process search, dimer,
    minimization) can be run at the same time.
    """
    num_jobs: int = Field(
        default=1,
        description="Number of jobs.",
    )
    """
    The meaning of this variable changes depending on the communicator type. For
    ``local``, it is the number of jobs run every time the program is invoked. For
    ``cluster``, it is the desired sum of the queued and running jobs.
    """
    max_jobs: int = Field(
        default=0,
        description="Maximum number of akmc jobs that can be running at once for the current state.",
    )
    """
    For communicators with queues (``cluster``), no more jobs will be queued if
    the number of jobs queued and in progress equals or exceeds this number. A
    default of 0 means unlimited.
    """
    client_path: str = Field(
        default="eonclient",
        description="Path to the eon client binary.",
    )
    """
    If only a name and not a path is given, ``eON`` looks for the binary in the
    same directory as the configuration file. If not found there, it searches
    through the directories in the ``$PATH`` environment variable.
    """
    number_of_CPUs: int = Field(
        default=1,
        description="Number of jobs that will run simultaneously for the local communicator.",
    )
    script_path: str = Field(
        default="./",
        description="Path to the user-defined scripts for submitting jobs to the communicator for the cluster communicator.",
    )
    name_prefix: str = Field(
        default="eon",
        description="Prefix added to job names to make them identifiable by the user for the cluster communicator.",
    )
    queued_jobs: str = Field(
        default="queued_jobs.sh",
        description="Name of the script that returns the job IDs of all the running and queued jobs for the cluster communicator.",
    )
    """
    This may return more than just ``eON`` jobs.
    """
    cancel_job: str = Field(
        default="cancel_job.sh",
        description="Name of the script that cancels a job for the cluster communicator.",
    )
    """
    Takes a single argument, the JobID.
    """
    submit_job: str = Field(
        default="submit_job.sh",
        description="Name of the script that submits a single job to the queuing system for the cluster communicator.",
    )
    """
    It takes two command line arguments. The first is the name of the job. This
    is not required for ``eON`` use, but is highly recommended so that users can
    identify which job is which. The second argument is the working directory.
    This is the path where the ``eON`` client should be executed. All of the
    needed client files will be placed in this directory. The script must return
    the job id of the submitted job. This is how ``eON`` internally keeps track
    of jobs.
    """


class ProcessSearchConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    minimization_offset: float = Field(
        default=None,
        description="This is the distance images bracketing the saddle are displaced.",
    )
    """
    After a saddle is found, images are placed on either side of the saddle
    along the mode and minimized to ensure that the saddle is connected to the
    original minimum and to locate the product state.
    This defaults to being the same as
    :any:`eon.schema.OptimizerConfig.max_move`
    """
    minimize_first: bool = Field(
        default=True,
        description="Every time a process search is run by a client the reactant will be minimized first before doing any saddle searches.",
    )


class PrefactorConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    all_free_atoms: bool = Field(
        default=False,
        description="Account for all free atoms when determining the prefactor.",
    )
    filter_scheme: Literal["cutoff", "fraction"] = Field(
        default="fraction",
        description="Determines how to filter the atoms for use in the prefactor calculation.",
    )
    """
    Options:
    - ``cutoff``: includes atoms that move more than :any:`eon.schema.PrefactorConfig.min_displacement`
    - ``fraction``: includes atoms that make up :any:`eon.schema.PrefactorConfig.filter_fraction` of the total motion, prioritizing the atoms that move the most
    """
    filter_fraction: float = Field(
        default=0.9,
        description="When using filter_scheme ``fraction``, includes the atoms that move the most, limited to the number that make up ``filter_fraction`` of the total motion.",
    )
    min_displacement: float = Field(
        default=0.25,
        description="Minimum displacement for an atom to be included in the Hessian calculation. Used only with filter_scheme 'cutoff'.",
    )
    within_radius: float = Field(
        default=3.3,
        description="Atoms within this radius of moving atoms are included in the Hessian. Used only with filter_scheme 'cutoff'.",
    )
    default_value: float = Field(
        default=0.0,
        description="Calculate prefactor if zero, otherwise use given value instead of doing a full prefactor calculation.",
    )
    min_value: float = Field(
        default=1e9, description="Minimum value for a reasonable prefactor."
    )
    max_value: float = Field(
        default=1e21, description="Maximum value for a reasonable prefactor."
    )
    configuration: Literal["reactant", "saddle", "product"] = Field(
        default="reactant",
        description="Configuration for which the eigenfrequencies will be determined in a prefactor job.",
    )
    """
    Options:
     - ``reactant``
     - ``saddle``
     - ``product``.
    """
    rate_estimation: str = Field(
        default="htst",
        description="Rate estimation method used for the prefactor calculation.",
    )


class PotentialConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    mpi_poll_period: float = Field(
        default=0.25, description="Polling period for MPI potential."
    )
    potential: Literal[
        "ams",
        "ams_io",
        "ase_orca",
        "bop",
        "bopfox",
        "cuh2",
        "eam_al",
        "edip",
        "emt",
        "ext",
        "fehe",
        "gpr",
        "imd",
        "lammps",
        "lenosky_si",
        "lj",
        "ljcluster",
        "morse_pt",
        "mpi",
        "mpi",
        "pyamff",
        "python",
        "qsc",
        "qsc",
        "spce",
        "sw_si",
        "tersoff_si",
        "tip4p",
        "tip4p_h",
        "tip4p_pt",
        "unknown",
        "vasp",
        "xtb",
        "zbl",
        "socket_nwchem",
    ] = Field(
        default="lj",
        description="Type of potential to execute.",
    )
    """
    Options:
     - ``ams``: Amsterdam Modeling Suite potential.
     - ``ams_io``: Amsterdam Modeling Suite via the I/O.
     - ``ase_orca``: ASE interface for ORCA quantum chemistry package.
     - ``bop``: Bond order potential for metals. [unused]
     - ``bopfox``: Bond order potential, for metals. [unused]
     - ``cuh2``: Potential for copper hydride systems.
     - ``eam_al``: Embedded atom method parameterized for aluminum.
     - ``edip``: Environment-Dependent Interatomic Potential, for carbon.
     - ``emt``: Effective medium theory, for metals.
     - ``ext``: External potential with system call interface.
     - ``fehe``: Potential for iron-hydrogen systems.
     - ``gpr``: Gaussian process regression potential.
     - ``imd``: IMD simulation package interface.
     - ``lammps``: The LAMMPS potentials.
     - ``lenosky_si``: Lenosky potential, for silicon.
     - ``lj``: Lennard-Jones potential in reduced units.
     - ``ljcluster``: Lennard-Jones cluster potential.
     - ``morse_pt``: Morse potential for platinum.
     - ``mpi``: Communicate with an MPI process to calculate energy and forces.
     - ``pyamff``: Python implementation of the AMFF potential.
     - ``python``: Custom python potential.
     - ``qsc``: Quantum Sutton-Chen potential, for FCC metals.
     - ``spce``: Simple Point Charge model for water.
     - ``sw_si``: Stillinger-Weber potential, for silicon.
     - ``tersoff_si``: Tersoff pair potential with angular terms, for silicon.
     - ``tip4p``: Point charge model for water.
     - ``tip4p_h``: TIP4P model for water with hydrogen.
     - ``tip4p_pt``: TIP4P model for water on platinum.
     - ``unknown``: Placeholder for unknown potential type.
     - ``vasp``: Vienna Ab-Initio Simulation Program (VASP) interface.
     - ``xtb``: Extended Tight Binding model.
    """
    log_potential: Optional[bool] = Field(
        default=None,
        description="If true, write timing information about each force call to client.log.",
    )

    @validator("log_potential", always=True)
    def set_log_potential(cls, v, values):
        if v is None:
            potential = values.get("potential")
            if potential in {"mpi", "vasp", "bop", "bopfox"}:
                return True
            else:
                return False
        return v


class XTBPot(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    paramset: str = Field(
        default="GFNFF", description="Parameter set for XTB potential."
    )
    accuracy: float = Field(default=1.0, description="Accuracy of the XTB calculation.")
    electronic_temperature: float = Field(
        default=0.0, description="Electronic temperature for XTB."
    )
    max_iterations: int = Field(
        default=250, description="Maximum number of XTB iterations."
    )


class ZBLPot(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    cut_inner: float = Field(
        default=2.0, description="Distance where switching function begins."
    )
    cut_global: float = Field(
        default=2.5, description="Global cutoff for the ZBL interaction."
    )


class SocketNWChemPot(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    host: str = Field(
        default="127.0.0.1", description="Host where the NWChem client will connect."
    )
    port: int = Field(
        default=9999, description="Port for the NWChem client to send data to EON."
    )
    nwchem_settings: str = Field(
        default="nwchem_settings.nwi",
        description="Path to the user-provided file containing NWChem scientific settings (basis, theory, task, etc.).",
    )
    unix_socket_path: str = Field(
        default="eon_nwchem",
        description="The basename for the UNIX socket file. The full path will be /tmp/ipi_<basename>.",
    )
    unix_socket_mode: bool = Field(
        default=False,
        description="If true, use a UNIX domain socket for communication instead of TCP/IP.",
    )
    mem_in_gb: int = Field(
        default=2,
        description="Memory (in GB) to be allocated for the NWChem calculation.",
    )
    make_template_input: bool = Field(
        default=True,
        description="If false, the socket input for the right geometry is generated by the user.",
    )


class Metatomic(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    model_path: str = Field(default="", description="Path to the Metatomic model.")
    device: str = Field(
        default="auto",
        description="Device to use for Metatomic calculations (e.g., 'cpu', 'cuda', 'auto').",
    )
    length_unit: str = Field(
        default="Angstrom", description="Length unit used by the Metatomic model."
    )
    extensions_directory: str = Field(
        default="", description="Directory for Metatomic extensions."
    )
    check_consistency: bool = Field(
        default=False,
        description="Whether to check consistency of the Metatomic model.",
    )


class ASE_NWCHEM(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    nwchem_path: str = Field(
        default="NONE", description="Path to the NWChem executable."
    )
    nproc: Union[int, Literal["auto"]] = Field(
        default="auto",
        description="Number of processors to use for NWChem. Can be 'auto' or an integer string.",
    )
    multiplicity: str = Field(
        default="1", description="Spin multiplicity for the NWChem calculation."
    )
    scf_thresh: float = Field(
        default=1e-5, description="SCF convergence threshold for NWChem."
    )
    scf_maxiter: int = Field(
        default=200, description="Maximum number of SCF iterations for NWChem."
    )


class ASE_ORCA(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    orca_path: str = Field(default="NONE", description="Path to the ORCA executable.")
    nproc: Union[int, Literal["auto"]] = Field(
        default="auto",
        description="Number of processors to use for ORCA. Can be 'auto' or an integer string.",
    )
    simpleinput: str = Field(
        default="ENGRAD HF-3c",
        description="Simple input string for ORCA, specifying method and basis set.",
    )


class AMSConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    engine: Literal["", "REAXFF", "MOPAC"] = Field(
        default="", description="Engine for AMS calculation. One of REAXFF, MOPAC."
    )
    forcefield: str = Field(
        default="", description="Force field to use (e.g., 'OPt.ff')."
    )
    model: str = Field(default="", description="Model to use (e.g., 'PM7', 'PM3').")
    xc: str = Field(default="", description="Exchange-correlation functional.")
    basis: str = Field(
        default="",
        description="Basis set to use with the exchange-correlation functional.",
    )
    resources: str = Field(default="", description="Resources path for DFTB.")


class AMSIOConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    engine: Literal["", "REAXFF", "MOPAC"] = Field(
        default="", description="Engine for AMS calculation. One of REAXFF, MOPAC."
    )
    forcefield: str = Field(
        default="", description="Force field to use (e.g., 'OPt.ff')."
    )
    model: str = Field(default="", description="Model to use (e.g., 'PM7', 'PM3').")
    xc: str = Field(default="", description="Exchange-correlation functional.")


class AMSEnvConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    amshome: str = Field(
        default="",
        description="Path to AMS home directory (e.g., '/some/path/to/amshome/').",
    )
    scm_tmpdir: str = Field(
        default="", description="Temporary directory for SCM (e.g., '/tmp')."
    )
    scm_pythondir: str = Field(
        default="", description="Python directory for SCM (e.g., '/.scm/python')."
    )
    amsbin: str = Field(
        default="",
        description="Path to AMS binary directory (will be appended to amshome).",
    )
    scmlicense: str = Field(
        default="",
        description="Path to SCM license file (will be appended to amshome).",
    )
    amsresources: str = Field(
        default="",
        description="Path to AMS atomic data resources (will be appended to amshome).",
    )


class SaddleSearchConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    method: Literal["min_mode", "dynamics"] = Field(
        default="min_mode", description="Method to locate the saddle point."
    )
    """
     Options:
      - ``min_mode``: Use a min-mode following scheme to locate the saddle
        point.
      - ``dynamics``: Experimental method that uses molecular dynamics to find
        new states and then runs a climbing image NEB calculation to find the
        saddle and a dimer calculation to estimate the eigenmode at the saddle.
    """
    min_mode_method: Literal["dimer", "lanczos", "gprdimer"] = Field(
        default="dimer", description="Min-mode method to use."
    )
    """
    Options:
     - ``dimer``: Use the dimer min-mode method from :cite:t:`ss-henkelmanDimerMethodFinding1999`
     - ``lanczos``: Use the Lanczos min-mode method from :cite:t:`ss-malekDynamicsLennardJonesClusters2000`
     - ``gprdimer``: Use the GP accelerated dimer method.
     """
    max_energy: float = Field(
        default=20.0,
        description="The energy at which a saddle search is considered bad and terminated.",
    )
    displace_radius: float = Field(
        default=5.0,
        description="Atoms within this distance of the epicenter will be displaced.",
    )
    displace_magnitude: float = Field(
        default=0.1,
        description="The standard deviation of the Gaussian displacement in each degree of freedom for the selected atoms.",
    )
    displace_random_weight: float = Field(
        default=1.0,
        description="Relative probability to displace with a random epicenter.",
    )
    displace_not_FCC_HCP_weight: float = Field(
        default=0.0,
        description="Relative probability to displace with an epicenter that is not FCC or HCP coordinated.",
    )
    displace_least_coordinated_weight: float = Field(
        default=0.0,
        description="Relative probability to displace with an epicenter that has a coordination number equal to the least-coordinated atom in the configuration.",
    )
    displace_under_coordinated_weight: float = Field(
        default=0.0,
        description="Relative probability to displace with an epicenter with a coordination equal to or less than displace_max_coordination.",
    )
    displace_softest_mode_weight: float = Field(
        default=0.0,
        description="Relative probability to displace the system along the softest vibrational mode.",
    )
    displace_listed_atom_weight: float = Field(
        default=0.0,
        description="Relative probability to displace with an epicenter listed in displace_atom_list.",
    )
    displace_atom_list: list[int] = Field(
        default=[-1],
        description="The individual index should be separated by a comma. Example: 10, 20, -1 would be the 10, 20, and the last atom.",
    )
    displace_listed_type_weight: float = Field(
        default=0.0,
        description="Relative probability to displace with an epicenter listed in displace_type_list.",
    )
    displace_type_list: list[str] = Field(
        default=[], description="The atom types should be separated by a comma."
    )
    displace_all_listed: bool = Field(
        default=False,
    )
    """
    This can be disabled by setting `displace_radius` to 0. Otherwise:

    - If true, each displacement will include all of the degrees of freedom of all of the listed atoms in `displace_atom_list` or `displace_type_list`.
    - If false, one of the atoms in `displace_atom_list` or `displace_type_list` will be selected at random for each displacement.
    - In either case, all atoms up to `displace_radius` distance away from any displaced atom will be included in the displacement.
    """
    displace_max_coordination: int = Field(
        default=11,
        description="When using under_coordinated as the displacement type, choose only atoms with a coordination equal to or less than this.",
    )
    converged_force: Optional[float] = Field(
        default=None,
        description="When the maximum force (in eV/A) on any one atom is smaller than this value, the structure is considered converged onto a saddle point.",
    )
    max_iterations: Optional[int] = Field(
        default=None, description="The maximum number of translation steps to be taken."
    )
    nonlocal_count_abort: int = Field(
        default=0,
        description="If this is not zero, the saddle search will abort when this many atoms have moved more than nonlocal_distance_abort from the initial displacement.",
    )
    nonlocal_distance_abort: float = Field(
        default=0.0,
        description="If nonlocal_count_abort is not zero, the saddle search will abort when nonlocal_count_abort atoms have moved more than this distance.",
    )
    client_displace_type: Literal[
        "load",
        "random",
        "last_atom",
        "min_coordinated",
        "not_fcc_or_hcp",
        "softest_mode",
    ] = Field(default="load", description="Type of displacement method used.")
    zero_mode_abort_curvature: float = Field(
        default=0.0,
        description="The saddle search will abort when the magnitude of the minmode curvature is less than this value.",
    )
    confine_positive: bool = Field(
        default=False,
        description="Activates a confinement scheme when the search is within a positive region of the PES.",
    )
    bowl_breakout: bool = Field(
        default=False,
        description="When activated, the search within positive regions of PES is confined to a subset of atoms.",
    )
    """
    Determines :any:`bowl_active_atoms` that are subject to the largest forces.
    To activate, :any:`confine_positive` must also be true. Method of
    :cite:t:`ss-pedersenBowlBreakoutEscaping2014`.
    """
    bowl_active_atoms: int = Field(
        default=20,
        description="Size of the applied confinement in the bowl breakout scheme.",
    )
    dynamics_temperature: Optional[float] = Field(
        default=None,
        description="The temperature, in Kelvin, for the molecular dynamics run.",
    )
    """
    A good initial choice might be near the melting temperature of the material.
    """
    dynamics_state_check_interval: float = Field(
        default=100.0,
        description="The time interval, in femtoseconds, to minimize the geometry and check if the system has left the initial state.",
    )
    dynamics_record_interval: float = Field(
        default=10.0,
        description="The time interval, in femtoseconds, between snapshots of the molecular dynamics trajectory.",
    )
    """
    Snapshots of MD trajectories are used to locate when the system first left
    the initial state. A binary search is used to locate the first snapshot that
    minimizes to a new geometry.
    """
    dynamics_linear_interpolation: bool = Field(
        default=True,
    )
    """
     - If true, then the band connecting the initial and final states will be
       initialized using a linear interpolation.
     - If false, then the band is interpolated through the first snapshot that
       minimizes to the final state.
    """
    dynamics_max_init_curvature: float = Field(
        default=0.0,
        description="The maximum initial curvature for the dynamics method in eV/Å^2.",
    )
    zero_mode_abort_curvature: float = Field(
        default=0.0,
        description="The saddle search will abort when the magnitude of the minmode curvature is less than this value.",
    )
    nonnegative_displacement_abort: bool = Field(
        default=False,
        description="Abort the search if a non-negative displacement is detected.",
    )
    max_single_displace: float = Field(
        default=10.0, description="The maximum single displacement value."
    )
    remove_rotation: bool = Field(
        default=False, description="Remove rotational components from the displacement."
    )
    perp_force_ratio: float = Field(
        default=0.0, description="The ratio of perpendicular force, undocumented."
    )
    confine_positive_min_force: float = Field(
        default=0.5,
        description="The minimum force for confining the positive region of the PES, undocumented.",
    )
    confine_positive_scale_ratio: float = Field(
        default=0.9,
        description="The scaling ratio for confining the positive region of the PES, undocumented.",
    )
    confine_positive_boost: float = Field(
        default=10.0,
        description="The boost factor for confining the positive region of the PES, undocumented.",
    )
    confine_positive_min_active: int = Field(
        default=30,
        description="The minimum number of active atoms for confining the positive region of the PES, undocumented.",
    )

    @validator("converged_force", always=True)
    def set_converged_force(cls, v, values):
        if v is None:
            optimizer_config = values.get("optimizer")
            if optimizer_config:
                return optimizer_config.converged_force
        return v

    @validator("max_iterations", always=True)
    def set_max_iterations(cls, v, values):
        if v is None:
            optimizer_config = values.get("optimizer")
            if optimizer_config:
                return optimizer_config.max_iterations
        return v


class KDBConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    use_kdb: bool = False
    kdb_only: bool = Field(
        default=False,
        description="Only use KDB in an AKMC simulation, no random searches.",
    )
    remove_duplicates: bool = Field(
        default=False,
        description="KDB will not make duplicate suggestions. This can slow KDB querying, so it may be best to use this only for slow potentials (DFT, etc.).",
    )
    kdb_name: str = "kdb.db"
    # TODO(rg): These are in config.yaml, not sure what they do..
    kdb_nf: float = 0.2
    kdb_dc: float = 0.3
    kdb_mac: float = 0.7


class RecyclingConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    use_recycling: bool = Field(default=True, description="Turn recycling on and off.")
    move_distance: float = Field(
        default=0.2,
        description="The distance an atom must move to be considered in the 'hole'.",
    )
    displace_moved_only: bool = Field(
        default=False,
        description="When performing saddle search displacements, only use atoms in the vicinity of the hole as displacement epicenters.",
    )
    save_suggestions: bool = Field(
        default=False,
        description="If True, the suggestions made by saddle recycling are saved to the 'saddle_suggestions' directory of the state directory for which the saddles are being suggested.",
    )
    active_region: float = Field(
        default=1.0,
        description="Defines a region around the atoms that moved in the previous KMC step to target new saddle searches in the current state.",
    )
    """
    For dynamics-based saddle searches, atoms in the active region will have
    mass weighting applied; for min-mode following saddle searches, only these
    atoms will be selected for displacements. Represents a multiple of
    :any:`eon.schema.StructureComparisonConfig.neighbor_cutoff`.
    """
    mass_weight_factor: float = Field(
        default=1000.0,
        description="Atoms that are not within the active region have their atomic masses multiplied by this factor.",
    )
    use_sb_recycling: bool = Field(
        default=False, description="Turn superbasin recycling on and off."
    )


class CoarseGrainingConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    use_mcamc: bool = Field(
        default=False,
        description="This option determines whether the Monte Carlo with Absorbing Markov Chains (MCAMC) coarse graining method will be used.",
    )
    """
    Mutually exclusive with :any:`eon.schema.CoarseGrainingConfig.use_askmc`.
    """
    state_file: str = Field(
        default="superbasin",
        description="File name for the state-specific data stored within each of the state directories.",
    )
    superbasin_scheme: Literal["energy_level", "transition_counting", "rate"] = Field(
        default="transition_counting",
        description="MCAMC provides a method for calculating transition rates across superbasins. An additional method is needed in order to decide when to combine states into a superbasin.",
    )
    """
    Options:
     - ``transition_counting``: Counts the number of times that the simulation has transitioned between a given pair of states. After a critical number of transitions have occurred, the pair of states are merged to form a superbasin.
     - ``energy_level``: States are merged based on energy levels filling up existing basins.
     - ``rate``: States are merged based only on rate criteria.
    """
    max_size: int = Field(
        default=0,
        description="The maximal number of states that will be merged together. If 0, there is no limit.",
    )
    number_of_transitions: int = Field(
        default=5,
        description="This is the number of transitions that must occur between two states before they are merged into a superbasin.",
    )
    """
    Only used if :any:`eon.schema.CoarseGrainingConfig.superbasin_scheme` is
    ``transition_counting``.
    """
    energy_increment: float = Field(
        default=0.01,
        description="The first time each state is visited, it is assigned an energy level first equal to the energy of the minimum. Every time the state is visited again by the Monte Carlo simulation, the energy level is increased by this amount.",
    )
    """
    Only used if :any:`eon.schema.CoarseGrainingConfig.superbasin_scheme` is
    ``energy_level``.
    """
    rate_threshold: float = Field(
        default=1e9,
        description="Any state with a rate (in 1/time) greater than this is merged into a single basin.",
    )
    """
    Only used if :any:`eon.schema.CoarseGrainingConfig.superbasin_scheme` is
    ``rate``.
    """
    superbasin_confidence: bool = Field(
        default=True,
        description="Superbasin KMC confidence.",
    )
    """
    Superbasin KMC steps only consider exit processes from the superbasin. As
    fast processes get absorbed more and more into the superbasin, the relevant
    exit processes have higher and higher barriers. The confidence to have found
    all relevant processes leading away from a state is heavily influenced by
    the fast processes, that no longer exit the superbasin. To make sure that
    the confidence is adequately high for the barrier energies leaving the
    superbasin, additional searches need to be performed. This setting (which
    defaults to true) enables these additional searches. The searches are marked
    in states/<state_number>/search_results.txt by appending the number of the
    superbasin in which the search was performed in brackets. In general, this
    option should not be disabled! It exists only for debug purposes and cases
    where the user is sure that the additional searches are not needed.
    """
    use_askmc: bool = Field(
        default=False,
        description="This option determines whether the AS-KMC coarse graining method will be used.",
    )
    """
    Mutually exclusive with :any:`eon.schema.CoarseGrainingConfig.use_mcamc`.
    """
    askmc_confidence: float = Field(
        default=0.9,
        description="The confidence for AS-KMC. This value determines the accuracy of the direction of the dynamics trajectory.",
    )
    askmc_barrier_raise_param: float = Field(
        default=1.5,
        description="This parameter sets how much the barriers are raised during AS-KMC.",
    )
    askmc_high_barrier_def: float = Field(
        default=2.0,
        description="This parameter sets how high a barrier must be to be considered 'high' in AS-KMC.",
    )
    askmc_barrier_test_on: bool = Field(
        default=True,
        description="This test verifies that no low-barrier process in a superbasin has been missed, considering even processes which have not been visited.",
    )
    """
    Because the implemented Superbasin Criterion actually only considers
    processes which have been passed over at least once, there is some chance
    that a low-barrier process in a superbasin might have not been visited at
    all while all other low-barrier processes have been visited at least
    :math:`N_f` times. This is unlikely, but this test verifies that such has
    not happened, considering even processes which have not been visited when
    determining if the Superbasin Criterion has first, because the implemented
    Superbasin Criterion actually only considers processes which have been
    passed over at least once, there is some chance that a low-barrier process
    in a superbasin might have not been visited at all while all other
    low-barrier processes have been visited at least :math:`N_f` times. This is
    unlikely, but this test verifies that such has not happened, considering
    even processes which have not been visited when determining if the
    Superbasin Criterion has passed. This check should not add significant
    overhead.
    """
    askmc_connections_test_on: Optional[bool] = Field(
        default=False,
        description="Ensures that there are no processes which connect states in the defined superbasin which have not been visited yet and which have a low-barrier.",
    )
    """
    This parameter determines whether to ensure that there are no processes
    which connect states in the defined superbasin which have not been visited
    yet and which have a low-barrier. This check is somewhat more
    computationally expensive than the previous because structure comparisons
    have to be made when finding product states of unvisited processes.
    """


class OptimizerConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    opt_method: Literal["box", "cg", "qm", "lbfgs", "fire"] = Field(
        default="cg",
        description="The optimization method to use.",
    )
    """
    Options:
      - ``box``: Optimizes the atom positions and box using quickmin
      - ``cg``: Conjugate gradient
      - ``qm``: Quickmin
      - ``lbfgs``: Limited Memory Broyden-Fletcher-Goldfarb-Shanno QuasiNewton optimizer
      - ``fire``: Fast inertial relaxation engine
    """
    convergence_metric: Literal["norm", "max_atom", "max_component"] = Field(
        default="norm",
        description="The metric to use to determine when an optimization is complete.",
    )
    """
    Options:
      - ``norm``: The norm of the entire force vector
      - ``max_atom``: The maximum force on any non-frozen atom
      - ``max_component``: The maximum force on any non-frozen degree of freedom
    """
    converged_force: float = Field(
        default=0.01,
        description="When the convergence_metric is smaller than this value (eV/A), the structure is considered minimized.",
    )
    max_move: float = Field(
        default=0.2,
        description="Maximum distance that an atom may be moved in a single optimization step (Angstroms).",
    )
    time_step: float = Field(
        default=1.0,
        description="The dynamical timestep for the quickmin algorithm (fs).",
    )
    max_iterations: int = Field(
        default=1000,
        description="The maximum number of optimization iterations that will be performed.",
    )


class QuickMinConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    time_step: float = Field(
        default=1.0, description="Time step for QuickMin, in femtoseconds."
    )
    qm_steepest_descent: bool = Field(
        default=False,
        description="If true, use the steepest descent method for QuickMin.",
    )


class FIREConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    time_step: float = Field(
        default=1.0, description="Time step for FIRE, in femtoseconds."
    )
    time_step_max: float = Field(
        default=1.0, description="Maximum time step for FIRE, in femtoseconds."
    )


class LBFGSConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    lbfgs_memory: int = Field(
        default=5,
        description="Number of previous gradients and positions to store for the LBFGS method.",
    )
    lbfgs_inverse_curvature: float = Field(
        default=1e-3, description="Initial inverse curvature value for LBFGS."
    )
    lbfgs_max_inverse_curvature: float = Field(
        default=1e-2, description="Maximum inverse curvature value for LBFGS."
    )
    lbfgs_auto_scale: bool = Field(
        default=True, description="If true, auto-scale the inverse curvature in LBFGS."
    )
    lbfgs_angle_reset: bool = Field(
        default=False, description="If true, reset the LBFGS angle."
    )
    lbfgs_distance_reset: bool = Field(
        default=False, description="If true, reset the LBFGS distance."
    )


class CGConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    cg_no_overshooting: bool = Field(
        default=False, description="If true, prevent overshooting in CG."
    )
    cg_knock_out_max_move: bool = Field(
        default=False, description="If true, limit the maximum move in CG."
    )
    cg_line_search: bool = Field(
        default=False, description="If true, perform a line search in CG."
    )
    cg_line_converged: float = Field(
        default=1e-4, description="Convergence criterion for the line search in CG."
    )
    cg_max_iter_before_reset: int = Field(
        default=20, description="Maximum number of iterations before reset in CG."
    )
    cg_max_iter_line_search: int = Field(
        default=10,
        description="Maximum number of iterations for the line search in CG.",
    )


class SDConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    sd_alpha: float = Field(default=0.1, description="Alpha value for SD.")
    sd_twopoint: bool = Field(
        default=False, description="If true, use the two-point method in SD."
    )


class RefineConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    refine_opt_method: Literal["none", "cg", "lbfgs", "fire", "box", "qm"] = Field(
        default="none",
        description="The optimization method to use for refinement.",
    )
    """
    Options:
      - ``none``: No refinement optimization
      - ``cg``: Conjugate gradient
      - ``lbfgs``: Limited Memory Broyden-Fletcher-Goldfarb-Shanno QuasiNewton optimizer
      - ``fire``: Fast inertial relaxation engine
      - ``box``: Optimizes the atom positions and box using quickmin
      - ``qm``: Quickmin
    """
    refine_threshold: float = Field(
        default=0.5, description="Threshold for refinement optimization."
    )


class DebugConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    save_stdout: bool = Field(
        default=False,
        description="Save the standard output from the client to a file named stdout_0.dat.",
    )
    interactive_shell: bool = Field(
        default=True, description="Enable or disable the interactive shell."
    )
    result_files_path: str = Field(
        default="./debug_results/", description="Where to store all result files."
    )
    use_mean_time: bool = Field(
        default=False,
        description="Select transition times from the mean of the exponential distribution of escape times.",
    )
    register_extra_results: bool = Field(
        default=False,
        description="Register processes found for a state after leaving that state.",
    )
    target_trajectory: bool = Field(
        default=False,
        description="Follow the state-to-state trajectory of another AKMC simulation.",
    )
    keep_bad_saddles: bool = Field(
        default=False,
        description="Keep data about bad saddles. If true, the result files for failed saddle searches are kept in the badprocdata directory within the state directory for that search.",
    )
    keep_all_result_files: bool = Field(
        default=False, description="Stores all result files in main_directory/results."
    )
    estimate_neb_eigenvalues: bool = Field(
        default=False,
        description="Write out the estimated lowest eigenvalue for each image of the NEB. This can significantly slow down the calculation!",
    )
    neb_mmf_estimator: Literal["dimer", "lanczos"] = Field(
        default="dimer",
        description="Estimator used to get the lowest eigenvalue, only used if estimate_neb_eigenvalues is true.",
    )
    write_movies: bool = Field(
        default=False, description="Output a movie of the calculation."
    )
    write_movies_interval: int = Field(
        default=1, description="Write a movie frame every write_movies_interval steps."
    )
    # TODO(rg): Document, from config.yaml
    stop_criterion: float = Field(
        default=1e8, description="Criterion to stop the calculation."
    )


class DimerConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    opt_method: Literal["sd", "cg", "lbfgs"] = Field(
        default="cg",
        description="Optimization algorithm to choose the dimer rotation direction.",
    )
    """
    Options:
     - 'sd': Steepest descent, rotate along the rotational force.
     - 'cg': Conjugate gradient, rotate along conjugate directions.
     - 'lbfgs': Limited memory Broyden-Fletcher-Goldfarb-Shanno Quasi-Newton optimizer.
    """

    converged_angle: float = Field(
        default=5.0,
        description="The dimer is considered converged if it will be rotated fewer degrees than this angle.",
    )
    rotations_max: int = Field(
        default=10,
        description="This is the maximum number of rotations allowed for the dimer for each minimum mode estimation.",
    )

    dimer_rotation_angle: float = Field(
        default=0.005, description="Finite angle for dimer rotation."
    )
    dimer_improved: bool = Field(
        default=True, description="Indicates if the improved dimer method is used."
    )
    dimer_max_iterations: int = Field(
        default=1000, description="Maximum number of iterations allowed for the dimer."
    )
    dimer_rotations_min: int = Field(
        default=1,
        description="Minimum number of rotations for the dimer. [not improved]",
    )
    dimer_torque_min: float = Field(
        default=0.1, description="Minimum torque for the dimer. [not improved]"
    )
    dimer_torque_max: float = Field(
        default=1.0, description="Maximum torque for the dimer. [not improved]"
    )
    dimer_remove_rotation: bool = Field(
        default=False, description="Indicates if the rotation should be removed."
    )
    """
    Implements the method of :cite:t:`dm-melanderRemovingExternalDegrees2015`
    """


class NudgedElasticBandConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    images: int = Field(
        default=5, description="Number of NEB images between the fixed endpoints."
    )
    spring: float = Field(
        default=5.0, description="The spring constant, in eV/Ang^2 between the images."
    )
    climbing_image_method: bool = Field(
        default=False, description="Indicates if the climbing image method is used."
    )
    """
    As discussed in :cite:t:`neb-henkelmanClimbingImageNudged2000`.
    """
    old_tangent: bool = Field(
        default=False, description="Indicates if the old tangent method is used."
    )
    """
    From :cite:t:`neb-millsQuantumThermalEffects1994`, before :cite:t:`neb-henkelmanImprovedTangentEstimate2000`.
    """
    neb_max_iterations: int = Field(
        default=1000, description="Maximum number of iterations allowed for the NEB."
    )
    neb_climbing_image_converged_only: bool = Field(
        default=True,
        description="Indicates if only the climbing image converged is used.",
    )
    neb_doubly_nudged: bool = Field(
        default=False, description="Indicates if the doubly nudged method is used."
    )
    neb_doubly_nudged_switching: bool = Field(
        default=False,
        description="Indicates if the doubly nudged switching method is used.",
    )
    """
    Method as demonstrated in :cite:t:`neb-trygubenkoDoublyNudgedElastic2004`.
    """
    neb_elastic_band: bool = Field(
        default=False, description="Indicates if the elastic band method is used."
    )
    neb_converged_force: float = Field(
        default=0.01,  # Assuming `optConvergedForce` is 0.0 as it's not provided
        description="Converged force threshold for the NEB.",
    )
    """
    This defaults to being the same as :any:`eon.schema.OptimizerConfig.converged_force`
    """
    energy_weighted: bool = Field(
        default=False, description="Indicates if the energy-weighted method is used."
    )
    """
    Method as demonstrated in :cite:t:`neb-asgeirssonNudgedElasticBand2021`.
    """
    ew_ksp_min: float = Field(
        default=0.5, description="Minimum value for KSP in the energy-weighted method."
    )
    ew_ksp_max: float = Field(
        default=5.0, description="Maximum value for KSP in the energy-weighted method."
    )
    initial_path_in: str = Field(
        default="",
        description="File from which the initial path is read.",
    )
    """
    This file must contain a list of .con files, one per image on the path.
    """
    minimize_endpoints: bool = Field(
        default=True,
        description="Minimize the reactant and product before the NEB.",
    )
    ci_after: float = Field(
        default=math.inf,
        description="Convergence before the CI is turned on.",
    )
    ci_mmf: bool = Field(
        default=False,
        description="Use an MMF method for a few steps at the CI.",
    )
    ci_mmf_after: float = Field(
        default=0.5,
        description="Convergence before the CI is turned on.",
    )
    ci_mmf_nsteps: int = Field(
        default=10,
        description="Number of steps for which the MMF is run at the CI image.",
    )


class LanczosConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    tolerance: float = Field(
        default=0.01,
        description="This is the convergence criteria for relative change in the estimated lowest eigenvalue.",
    )
    max_iterations: int = Field(
        default=20,
        description="The maximum number of refinement iterations when calculating the minimum eigenvalue.",
    )
    quit_early: bool = Field(
        default=True,
        description="If the relative change between the previous lowest eigenvalue and the curvature along the initial direction is less than the tolerance, terminate.",
    )


class HessianConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)
    atom_list: Union[str, list[int]] = Field(
        default="All",
        description="The atoms that will be displaced in the calculation of the Hessian: a comma delimited list of atom indices, e.g. 0,1,2. Default is 'All'.",
    )
    zero_freq_value: float = Field(
        default=1e-6, description="The value assigned to zero frequencies."
    )


class DynamicsConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    time_step: float = Field(
        default=1.0, description="The duration of each MD step, in femtoseconds."
    )
    time: float = Field(
        default=1000.0,
        description="Total MD time, in femtoseconds.",
    )
    thermostat: Literal["none", "andersen", "langevin", "nose_hoover"] = Field(
        default="none", description="Thermostat to use for the dynamics simulation."
    )
    """
    Options:
    - ``none``: NVE dynamics with the verlet algorithm. Initial velocities set by temperature.
    - ``andersen``: Andersen thermostat with the Verlet algorithm.
    - ``langevin``: Langevin thermostat with the Verlet algorithm.
    - ``nose_hoover``: Nosé-Hoover thermostat with the Verlet algorithm.
    """
    andersen_collision_period: float = Field(
        default=100.0,
        description="The collision period (in fs) for the Andersen thermostat.",
    )
    andersen_alpha: float = Field(
        default=1.0, description="The collision strength in the Andersen thermostat."
    )
    nose_mass: float = Field(
        default=1.0,
        description="The effective mass of the additional degree of freedom in the Nosé-Hoover thermostat, which determines the coupling frequency of the thermostat.",
    )
    langevin_friction: float = Field(
        default=0.01,
        description="The damping coefficient for Langevin dynamics (1/fs).",
    )


class ParallelReplicaConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    dephase_time: float = Field(
        default=1000.0,
        description="Time (in fs) to decorrelate the replica trajectories.",
    )
    """
    The momenta will be inverted when reaching the dividing surface to prevent
    transitions occurring during this time.
    """

    state_check_interval: float = Field(
        default=1000.0, description="Frequency (in fs) to check the system state."
    )
    """
    Every ``state_check_interval``, the current configuration is minimized and
    compared to the initial state minimum to decide if a transition to a new
    state has occurred. When :any:`refine_transition` is set to true, the code
    will keep a buffer of configurations so that the transition time can be
    determined to more precision than ``state_check_interval``.
    """

    refine_transition: bool = Field(
        default=True, description="Whether the transition time is refined."
    )
    """
    When this option is true, an array of (:any:`state_check_interval`/
    :any:`state_save_interval`) configurations along the history of the
    trajectory is saved. A binary search is then used to determine the
    transition time. When this option is false, the transition time is taken to
    be when the new state was found. This function reduces the need for small
    values of :any:`state_check_interval`; the precision of the transition time
    is :any:`state_save_interval` * :any:`dynamics time_step
    <eon.schema.DynamicsConfig.time_step>`.
    """

    state_save_interval: float = Field(
        default=0.1 * 1000.0,  # 0.1 * state_check_interval
        description="Frequency to record system state.",
    )
    """
    How often the system is recorded to the buffer array when the
    :any:`refine_transition` option is
    activated. Increasing the value of
    :any:`state_save_interval` lowers the
    precision of the transition time estimate but also reduces memory usage and
    speeds up refinement of the transition step.
    """

    post_transition_time: float = Field(
        default=1000.0,
        description="Additional time (in fs) after a new state is found.",
    )
    """
    A state check will be employed after this post_transition_time to confirm
    that the state is stable. This additional check helps avoid meta-stable
    states. A value similar to
    :any:`dephase_time` is recommended.
    """

    stop_after_transition: bool = Field(
        default=False, description="Whether to stop the job when a new state is found."
    )

    dephase_loop_stop: bool = Field(
        default=False, description="Whether to stop the dephase loop."
    )

    dephase_loop_max: int = Field(
        default=5, description="Maximum number of dephase loops."
    )


class HyperdynamicsConfig(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    bias_potential: Literal["none", "bond_boost"] = Field(
        default="none", description="Type of bias potential to use."
    )
    """
    Options:
    - ``none``: With no bias potential, run regular MD.
    - ``bond_boost``: Bond boost method.
    """

    bb_dvmax: float = Field(
        default=0.0, description="The magnitude of the bond-boost bias potential."
    )

    bb_rmd_time: float = Field(
        default=100.0,
        description="Regular MD duration (in fs) to determine equilibrium bond length before adding bias potential.",
    )

    bb_boost_atom_list: Union[str, list[int]] = Field(
        default="ALL",
        description="The atoms that will be displaced in the calculation: a comma delimited list of atom indices, e.g. 0,1,2. Default is 'All'.",
    )

    bb_rcut: float = Field(
        default=3.0,
        description="Cutoff distance (in Angstroms) for bonds included in the bond-boost potential.",
    )

    bb_stretch_threshold: float = Field(
        default=0.2, description="Defines the bond-boost dividing surface."
    )
    """
    It should be smaller than the maximum fractional nearest-neighbor bond stretch or compression at any transition state.
    """

    bb_ds_curvature: float = Field(
        default=0.95, description="Curvature near the bond-boost dividing surface."
    )
    """
    It should have a value <= 1. We recommend the value to be 0.9-0.98.
    """


class GprHyperparameterOptConfig(BaseModel):
    """
    Configuration settings for Gaussian Process Regression (GPR)
    hyperparameter optimization.
    """

    model_config = ConfigDict(use_attribute_docstrings=True)

    hyperparameter_opt_method: Literal["scg_opt", "adam_opt"] = Field(
        default="scg_opt", description="Hyperparameter optimization algorithm to use."
    )

    check_derivatives: bool = Field(
        default=False,
        description="Indicates if derivatives should be checked during optimization.",
    )

    opt_tol_func: float = Field(
        default=1e-4, description="The tolerance for the function value."
    )

    opt_tol_sol: float = Field(
        default=1e-4, description="The tolerance for the solution."
    )

    opt_max_iterations: int = Field(
        default=400, description="Maximum number of optimization iterations."
    )

    scg_lambda_init: float = Field(
        default=10.0, description="The initial lambda value for the SCG algorithm."
    )

    scg_lambda_limit: float = Field(
        default=1e20,
        description="The upper limit for the lambda value in the SCG algorithm.",
    )

    adam_learning_rate: float = Field(
        default=1e-3, description="The learning rate for the ADAM algorithm."
    )

    adam_learning_rate_decay: float = Field(
        default=0.999, description="The learning rate decay for the ADAM algorithm."
    )

    adam_beta1: float = Field(
        default=0.9,
        description="The beta1 value (exponential decay rate for the first moment estimates) for the ADAM algorithm.",
    )

    adam_beta2: float = Field(
        default=0.99,
        description="The beta2 value (exponential decay rate for the second moment estimates) for the ADAM algorithm.",
    )

    adam_epsilon: float = Field(
        default=1e-8,
        description="A small constant for numerical stability in the ADAM algorithm.",
    )

    adam_weight_decay: float = Field(
        default=0.0, description="The weight decay for the ADAM algorithm."
    )

    adam_amsgrad: bool = Field(
        default=True,
        description="Indicates if the AMSGrad variant of ADAM should be used.",
    )


class GPRDimerConfig(BaseModel):
    """
    Configuration settings for the Gaussian Process Regression (GPR) Dimer method.
    """

    model_config = ConfigDict(use_attribute_docstrings=True)

    # Dimer Core Parameters
    finite_angle: float = Field(
        default=0.005, description="Finite angle for dimer rotation."
    )
    converged_angle: float = Field(
        default=0.08,
        description="The dimer is considered converged if it rotates fewer degrees than this angle.",
    )
    relaxation_converged_angle: float = Field(
        default=0.001, description="Convergence angle for relaxation (T_anglerot_gp)."
    )
    max_initial_rotation_iterations: int = Field(
        default=6,
        description="Maximum number of initial rotation iterations (should be DoF).",
    )
    max_relaxation_rotation_iterations: int = Field(
        default=10,
        description="Maximum number of relaxation rotation iterations (num_iter_rot_gp).",
    )
    divisor_t_dimer: int = Field(default=10, description="Divisor for T_dimer_gp.")
    max_outer_iterations: int = Field(
        default=300, description="Maximum number of outer iterations (num_bigiter)."
    )
    max_inner_iterations: int = Field(
        default=1000, description="Maximum number of inner iterations (num_iter)."
    )
    max_midpoint_displacement: float = Field(
        default=0.5, description="Maximum midpoint displacement (disp_max)."
    )
    rotation_opt_method: str = Field(
        default="lbfgs", description="Optimization method for rotation (method_rot)."
    )
    translation_opt_method: str = Field(
        default="lbfgs",
        description="Optimization method for translation (method_trans).",
    )
    active_radius: float = Field(
        default=5.0, description="Active radius for atom selection (actidst_fro)."
    )
    dimer_separation: float = Field(default=0.01, description="Dimer separation.")
    convex_region_step_size: float = Field(
        default=0.1, description="Convex region step size."
    )
    max_step_size: float = Field(default=0.1, description="Maximum step size.")
    rotation_removal_projection_threshold: float = Field(
        default=0.1, description="Threshold below which torques are removed."
    )
    ratio_at_limit: float = Field(default=0.66667, description="Ratio at limit.")
    nogp_initial_rotations: bool = Field(
        default=False,
        description="Flag to indicate whether GP is used for initial rotations.",
    )
    nogp_init_translations: bool = Field(
        default=False,
        description="Flag to indicate whether GP is used for initial translations.",
    )
    has_many_iterations: bool = Field(
        default=True,
        description="Flag to indicate if many iterations are used (islarge_num_iter).",
    )

    # GPR Parameters
    gpr_variance: float = Field(default=1e-8, description="GPR variance (gp_sigma2).")
    gpr_jitter_variance: float = Field(
        default=0.0, description="GPR jitter variance (jitter_sigma2)."
    )
    gpr_noise_variance: float = Field(
        default=1e-8, description="GPR noise variance (sigma2)."
    )
    prior_mean: float = Field(default=0.0, description="Prior mean.")
    prior_variance: float = Field(default=1.0, description="Prior variance (prior_s2).")
    prior_degrees_of_freedom: float = Field(
        default=20.0, description="Prior degrees of freedom (prior_nu)."
    )

    # GPR Optimization Parameters (nested model)
    gpr_hypopt: GprHyperparameterOptConfig = Field(
        default_factory=GprHyperparameterOptConfig,
        description="Hyperparameter optimization settings for GPR.",
    )

    # Prune Parameters
    use_prune: bool = Field(
        default=False, description="Flag to indicate if pruning should be used."
    )
    start_prune_at: int = Field(
        default=8, description="The iteration at which to start pruning."
    )
    nprune_vals: int = Field(default=3, description="The number of values to prune.")
    prune_threshold: float = Field(default=0.5, description="The pruning threshold.")
    es_threshold: float = Field(default=1.2, description="Early stopping threshold.")
    es_dist_metric: Literal["emd", "rmsd", "max1DLog"] = Field(
        default="emd", description="The distance for early stopping."
    )

    # Debugging Parameters
    report_level: int = Field(default=1, description="The reporting level.")
    debug_level: int = Field(default=2, description="The debugging level.")
    debug_output_directory: str = Field(
        default="output", description="The debug output directory."
    )
    debug_position_basename: str = Field(
        default="position", description="The debug position file basename."
    )
    debug_energy_basename: str = Field(
        default="energy", description="The debug energy file basename."
    )
    debug_gradient_basename: str = Field(
        default="gradient", description="The debug gradient file basename."
    )
    debug_out_ext: str = Field(
        default="dat", description="The debug output file extension."
    )
    debug_midpoint_offset: float = Field(
        default=3.0, description="The debug offset from the midpoint."
    )
    debug_y_step: float = Field(default=0.1, description="The debug y step.")
    debug_z_step: float = Field(default=0.1, description="The debug z step.")


class Config(BaseModel):
    model_config = ConfigDict(use_attribute_docstrings=True)

    main: MainConfig
    structure_comparison: StructureComparisonConfig
    akmc: AKMCConfig
    basin_hopping: BasinHoppingConfig
    paths: PathsConfig
    dimer: DimerConfig
    gprd: GPRDimerConfig
    neb: NudgedElasticBandConfig
    lanczos: LanczosConfig
    hessian: HessianConfig
    communicator: CommunicatorConfig
    process_search: ProcessSearchConfig
    prefactor: PrefactorConfig
    saddle_search: SaddleSearchConfig
    potential: PotentialConfig
    refine: RefineConfig
    kdb: KDBConfig
    dynamics: DynamicsConfig
    parrep: ParallelReplicaConfig
    hyperdyn: HyperdynamicsConfig
    recycling: RecyclingConfig
    coarse_graining: CoarseGrainingConfig
    optimizer: OptimizerConfig
    debug: DebugConfig

    @validator("communicator")
    def check_communicator(cls, v, values):
        if v.type == "local":
            assert v.client_path is not None and v.number_of_CPUs is not None, (
                "local communicator requires client_path and number_of_CPUs"
            )
        if v.type == "cluster":
            assert (
                v.script_path is not None
                and v.name_prefix is not None
                and v.queued_jobs is not None
                and v.cancel_job is not None
                and v.submit_job is not None
            ), (
                "cluster communicator requires script_path, name_prefix, queued_jobs, cancel_job, and submit_job"
            )
        return v

    @validator("coarse_graining")
    def check_coarse_graining(cls, v, values):
        if v.superbasin_scheme == "transition_counting":
            assert v.number_of_transitions is not None, (
                "transition_counting scheme requires number_of_transitions"
            )
        if v.superbasin_scheme == "energy_level":
            assert v.energy_increment is not None, (
                "energy_level scheme requires energy_increment"
            )
        if v.use_askmc:
            assert (
                v.askmc_confidence is not None
                and v.askmc_barrier_raise_param is not None
                and v.askmc_high_barrier_def is not None
                and v.askmc_barrier_test_on is not None
                and v.askmc_connections_test_on is not None
            ), (
                "askmc requires askmc_confidence, askmc_barrier_raise_param, askmc_high_barrier_def, askmc_barrier_test_on, and askmc_connections_test_on"
            )
        return v

    @validator("process_search", always=True)
    def update_process_search(cls, v, values):
        optimizer_config = values.get("optimizer")
        if optimizer_config and v.minimization_offset is None:
            v.minimization_offset = optimizer_config.max_move
        return v
