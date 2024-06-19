from pydantic import BaseModel, Field, ValidationError, validator
from typing import Optional


class MainConfig(BaseModel):
    job: str
    temperature: float
    checkpoint: bool
    random_seed: int


class StructureComparisonConfig(BaseModel):
    energy_difference: float
    distance_difference: float
    indistinguishable_atoms: bool
    check_rotation: bool
    brute_neighbors: bool
    neighbor_cutoff: float
    use_covalent: bool
    covalent_scale: float
    remove_translation: bool


class AKMCConfig(BaseModel):
    confidence: float
    server_side_process_search: bool
    thermally_accessible_window: float
    thermally_accessible_buffer: float
    max_kmc_steps: int
    confidence_scheme: str
    confidence_correction: bool
    max_rate: float
    eq_rate: float


class BasinHoppingConfig(BaseModel):
    initial_random_structure_probability: float
    initial_state_pool_size: int


class PathsConfig(BaseModel):
    main_directory: str
    jobs_out: str
    jobs_in: str
    incomplete: str
    states: str
    results: str
    potential_files: str
    bh_minima: str
    kdb_scratch: str
    kdb: str
    superbasins: str
    superbasin_recycling: str
    scratch: str


class CommunicatorConfig(BaseModel):
    type: str
    jobs_per_bundle: int
    num_jobs: int
    max_jobs: int
    client_path: Optional[str]
    number_of_CPUs: Optional[int]
    script_path: Optional[str]
    name_prefix: Optional[str]
    queued_jobs: Optional[str]
    cancel_job: Optional[str]
    submit_job: Optional[str]


class ProcessSearchConfig(BaseModel):
    minimization_offset: float


class PrefactorConfig(BaseModel):
    default_value: float


class PotentialConfig(BaseModel):
    mpi_poll_period: float


class SaddleSearchConfig(BaseModel):
    method: str
    max_iterations: int
    dynamics_temperature: float
    displace_random_weight: float
    displace_listed_atom_weight: float
    displace_listed_type_weight: float
    displace_all_listed: bool
    displace_under_coordinated_weight: float
    displace_least_coordinated_weight: float
    displace_not_FCC_HCP_weight: float
    displace_not_TCP_BCC_weight: float
    displace_not_TCP_weight: float
    displace_water_weight: float
    stdev_translation: float
    stdev_rotation: float
    molecule_list: str
    disp_at_random: int
    displace_magnitude: float
    displace_radius: float
    displace_min_norm: float
    void_bias_fraction: float
    displace_max_coordination: int
    random_mode: bool
    displace_1d: bool
    dynamics_max_init_curvature: float
    zero_mode_abort_curvature: float
    displace_atom_list: str = None
    displace_type_list: str = None


class KDBConfig(BaseModel):
    use_kdb: bool
    kdb_only: bool
    kdb_scratch_path: str
    kdb_path: str
    remove_duplicates: bool
    kdb_name: str
    kdb_nf: str
    kdb_dc: str
    kdb_mac: str


class RecyclingConfig(BaseModel):
    use_recycling: bool
    save_suggestions: bool
    displace_moved_only: Optional[bool]
    move_distance: float
    active_region: float
    mass_weight_factor: float
    use_sb_recycling: bool
    superbasin_recycling: Optional[str]


class CoarseGrainingConfig(BaseModel):
    use_mcamc: bool
    state_file: str
    superbasin_scheme: str
    max_size: int
    number_of_transitions: Optional[int]
    energy_increment: Optional[float]
    superbasin_confidence: bool
    use_askmc: bool
    askmc_confidence: Optional[float]
    askmc_barrier_raise_param: Optional[float]
    askmc_high_barrier_def: Optional[float]
    askmc_barrier_test_on: Optional[bool]
    askmc_connections_test_on: Optional[bool]


class OptimizerConfig(BaseModel):
    max_iterations: int


class DebugConfig(BaseModel):
    interactive_shell: bool
    keep_bad_saddles: bool
    keep_all_result_files: bool
    result_files_path: str
    register_extra_results: bool
    use_mean_time: bool
    target_trajectory: str
    stop_criterion: float


class Config(BaseModel):
    main: MainConfig
    structure_comparison: StructureComparisonConfig
    akmc: AKMCConfig
    basin_hopping: BasinHoppingConfig
    paths: PathsConfig
    communicator: CommunicatorConfig
    process_search: ProcessSearchConfig
    prefactor: PrefactorConfig
    saddle_search: SaddleSearchConfig
    potential: PotentialConfig
    kdb: KDBConfig
    recycling: RecyclingConfig
    coarse_graining: CoarseGrainingConfig
    optimizer: OptimizerConfig
    debug: DebugConfig

    @validator("communicator")
    def check_communicator(cls, v, values):
        if v.type == "local":
            assert (
                v.client_path is not None and v.number_of_CPUs is not None
            ), "local communicator requires client_path and number_of_CPUs"
        if v.type == "cluster":
            assert (
                v.script_path is not None
                and v.name_prefix is not None
                and v.queued_jobs is not None
                and v.cancel_job is not None
                and v.submit_job is not None
            ), "cluster communicator requires script_path, name_prefix, queued_jobs, cancel_job, and submit_job"
        return v

    @validator("coarse_graining")
    def check_coarse_graining(cls, v, values):
        if v.superbasin_scheme == "transition_counting":
            assert (
                v.number_of_transitions is not None
            ), "transition_counting scheme requires number_of_transitions"
        if v.superbasin_scheme == "energy_level":
            assert (
                v.energy_increment is not None
            ), "energy_level scheme requires energy_increment"
        if v.use_askmc:
            assert (
                v.askmc_confidence is not None
                and v.askmc_barrier_raise_param is not None
                and v.askmc_high_barrier_def is not None
                and v.askmc_barrier_test_on is not None
                and v.askmc_connections_test_on is not None
            ), "askmc requires askmc_confidence, askmc_barrier_raise_param, askmc_high_barrier_def, askmc_barrier_test_on, and askmc_connections_test_on"
        return v
