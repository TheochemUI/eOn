"""pyeonclient — in-process eOn client (nanobind).

Class-first surface (import as ``pyec``)::

    import pyeonclient as pyec
    from pyeonclient.models import DimerSpec, Accelerant

    # Min-mode: default method=improved; GP only with improved (validated)
    d = pyec.Dimer(matter, params, pot)
    d = pyec.Dimer(matter, params, pot, accelerant="gp")
    d = pyec.Dimer(matter, params, pot, spec=DimerSpec(accelerant=Accelerant.gp))
    d.compute(matter, direction)

    # NEB: path init + optional GP accelerant (validated NebSpec)
    path = pyec.neb_idpp_path(initial, final, n, params)
    neb = pyec.NEB(path, params, pot)
    neb = pyec.NEB(initial, final, params, pot, accelerant="gp")

    # Jobs
    out = pyec.TAD(matter, params, pot).run(inplace=False)
    out = pyec.MolecularDynamics(matter, params).run()

Docs: https://eondocs.org/user_guide/pyeonclient.html
"""

from __future__ import annotations

try:
    from pyeonclient._core import (  # type: ignore F401
        IoStatus,
        Job,
        JobType,
        Matter,
        NEBInit,
        OptType,
        Parameters,
        PbcConvention,
        Potential,
        PotType,
        RunStatus,
        __version__,
        append_results_timing,
        get_process_times,
        io_ok,
        io_status_name,
        job_type_from_name,
        job_type_name,
        make_job,
        make_potential,
        make_potential_from_ase,
        ase_potential as _ase_potential_core,
        attach_ase_calculator as _attach_ase_calculator_core,
        bind_ase_matter,
        matter_from_ase,
        pot_type_from_name,
        pot_type_name,
        steady_clock_now,
        write_potcall_summary,
        # NEB
        NudgedElasticBand,
        NEB as _NEBCore,
        NEBStatus,
        neb_read_file_paths,
        neb_load_path_from_files,
        neb_linear_path,
        neb_idpp_path,
        neb_idpp_collective_path,
        neb_sidpp_path,
        neb_initial_path,
        neb_write_results,
        pot_registry_total_force_calls,
        # Min-mode engines (chemist Dimer/NEB come from api after validation)
        Dimer as _DimerCore,
        ClassicDimer,
        ImprovedDimer,
        Lanczos,
        Davidson,
        resolve_mobile_atoms,
        free_atom_indices,
        built_with_gprd,
        # Saddle
        MinModeSaddleSearch,
        SaddleStatus,
        saddle_status_message,
        min_mode_saddle_search,
        # Analysis
        Hessian,
        get_prefactors,
        moved_atoms,
        all_free_atoms,
        # Class-first jobs
        MolecularDynamics,
        MonteCarlo,
        BasinHopping,
        ProcessSearch,
        TAD,
        ParallelReplica,
        SafeHyperdynamics,
        ReplicaExchange,
        structures_equal,
        structure_distance,
        # Build probes
        built_with_metatomic,
        built_with_rgpot,
        built_with_gp_surrogate,
    )
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "pyeonclient._core is not built. Configure with "
        "-Dwith_pyeonclient=true (nanobind; stable ABI / free-threaded)."
    ) from e

from pyeonclient.api import Dimer, NEB

# Typed specs optional — install pydantic via pyeonclient[models]
try:
    from pyeonclient.models import (  # type: ignore F401
        Accelerant,
        DimerSpec,
        MinModeMethod,
        NebSpec,
        PathInit,
    )
except ImportError:  # pragma: no cover
    Accelerant = None  # type: ignore[misc, assignment]
    DimerSpec = None  # type: ignore[misc, assignment]
    MinModeMethod = None  # type: ignore[misc, assignment]
    NebSpec = None  # type: ignore[misc, assignment]
    PathInit = None  # type: ignore[misc, assignment]

from pyeonclient.bridge import (
    from_structure,
    matter_to_structure,
    structure_to_matter,
    to_structure,
)
from pyeonclient.ase_bridge import (  # noqa: E402
    ase_potential,
    ase_to_matter,
    attach_ase_calculator,
    potential_from_ase,
    ase_to_structure,
    conframe_to_matter,
    matter_to_ase,
    matter_to_conframe,
    structure_to_ase,
)

from pyeonclient.backends import (  # noqa: E402
    ensure_metatomic_load_compat,
    list_backends,
    load_atomistic_model_compat,
    make_backend,
    make_metatomic_ase_calculator,
    register as register_backend,
)

from pyeonclient.steps import (  # noqa: E402
    append_timing,
    load_parameters,
    minimize_workdir,
    run_eon_cwd,
    run_job,
    run_job_in_directory,
    write_minimization_results,
    write_neb_results,
    neb_workdir,
    rgpot_metatomic_workdir,
    rgpot_metatomic_neb_workdir,
)

# Optional gpr_optim overlay (GP-OIE / OCINEB via gpr_optim.session)
try:
    from pyeonclient.gpr_overlay import (  # type: ignore F401
        compute_gpneb_from_matters,
        gpr_optim_available,
        gpneb_session_from_matters,
    )
except ImportError:  # pragma: no cover
    gpr_optim_available = lambda: False  # type: ignore[misc, assignment]
    gpneb_session_from_matters = None  # type: ignore[misc, assignment]
    compute_gpneb_from_matters = None  # type: ignore[misc, assignment]

to_ase = matter_to_ase
from_ase = ase_to_matter

__all__ = [
    "IoStatus",
    "Job",
    "JobType",
    "Matter",
    "Parameters",
    "PbcConvention",
    "Potential",
    "PotType",
    "OptType",
    "RunStatus",
    "NEBInit",
    "NEBStatus",
    "NudgedElasticBand",
    "NEB",
    "neb_linear_path",
    "neb_idpp_path",
    "neb_idpp_collective_path",
    "neb_sidpp_path",
    "neb_initial_path",
    "neb_load_path_from_files",
    "neb_read_file_paths",
    "neb_write_results",
    "Dimer",
    "PathInit",
    "Accelerant",
    "MinModeMethod",
    "NebSpec",
    "DimerSpec",
    "ClassicDimer",
    "ImprovedDimer",
    "Lanczos",
    "Davidson",
    "resolve_mobile_atoms",
    "free_atom_indices",
    "built_with_gprd",
    "MinModeSaddleSearch",
    "SaddleStatus",
    "saddle_status_message",
    "min_mode_saddle_search",
    "Hessian",
    "get_prefactors",
    "moved_atoms",
    "all_free_atoms",
    "MolecularDynamics",
    "MonteCarlo",
    "BasinHopping",
    "ProcessSearch",
    "TAD",
    "ParallelReplica",
    "SafeHyperdynamics",
    "ReplicaExchange",
    "structures_equal",
    "structure_distance",
    "make_job",
    "make_potential",
    "make_potential_from_ase",
    "ase_potential",
    "attach_ase_calculator",
    "bind_ase_matter",
    "matter_from_ase",
    "potential_from_ase",
    "make_backend",
    "list_backends",
    "register_backend",
    "ensure_metatomic_load_compat",
    "load_atomistic_model_compat",
    "make_metatomic_ase_calculator",
    "load_parameters",
    "run_job",
    "run_job_in_directory",
    "run_eon_cwd",
    "minimize_workdir",
    "neb_workdir",
    "rgpot_metatomic_workdir",
    "rgpot_metatomic_neb_workdir",
    "write_minimization_results",
    "write_neb_results",
    "write_potcall_summary",
    "append_results_timing",
    "append_timing",
    "get_process_times",
    "steady_clock_now",
    "pot_registry_total_force_calls",
    "io_ok",
    "io_status_name",
    "job_type_from_name",
    "job_type_name",
    "pot_type_from_name",
    "pot_type_name",
    "built_with_metatomic",
    "built_with_rgpot",
    "built_with_gp_surrogate",
    "from_ase",
    "to_ase",
    "ase_to_matter",
    "matter_to_ase",
    "ase_to_structure",
    "structure_to_ase",
    "matter_to_conframe",
    "conframe_to_matter",
    "from_structure",
    "to_structure",
    "matter_to_structure",
    "structure_to_matter",
    "gpr_optim_available",
    "gpneb_session_from_matters",
    "compute_gpneb_from_matters",
    "__version__",
]
