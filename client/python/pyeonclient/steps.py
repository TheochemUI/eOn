"""Optional workdir composition helpers for ClientEON stages.

Prefer Matter-first composition for new code (no package alias required)::

    from pyeonclient import Matter, NEB, NebSpec, Parameters, PathInit
    from pyeonclient import (
        append_timing,
        pot_registry_total_force_calls,
        steady_clock_now,
        write_neb_results,
        write_potcall_summary,
    )
    from pyeonclient.backends import make_backend

    params = Parameters()
    NebSpec(n_images=10, path_init=PathInit.idpp, energy_weighted=True).apply_to_parameters(
        params
    )
    pot = make_backend("metatomic", model_path="pet-mad.pt", params=params)
    initial, final = Matter(pot, params), Matter(pot, params)
    # … con2matter endpoints …
    t0, f0 = steady_clock_now(), pot_registry_total_force_calls()
    neb = NEB(initial, final, params, pot)
    status = neb.compute()
    write_neb_results(neb, params, pot_registry_total_force_calls() - f0)
    write_potcall_summary("_potcalls.json")
    append_timing("results.dat", t0)

Workdir helpers below still load ``config.ini`` for eonclient-style batch
dirs (``make_potential`` from PotType).

Typical job path (any JobType including NEB)::

    params = load_parameters("config.ini")
    t0 = steady_clock_now()
    files = run_job(params)          # make_job + Job.run + drop Job
    write_potcall_summary()
    append_timing("results.dat", t0)

Matter path (minimization without Job wrapper)::

    params = load_parameters("config.ini")  # or Parameters() + make_backend
    pot = make_potential(params.potential, params)  # INI / PotType path
    # Class-first: pot = make_backend("metatomic", model_path=..., params=params)
    m = Matter(pot, params)
    m.con2matter("pos.con")
    _m, ok = m.relax(inplace=True, write_movie=True,
                     prefix_movie="minimization", prefix_checkpoint="pos")
    m.matter2con("min.con")
    write_minimization_results(params, pot, m, converged=ok)
    del m, pot
    write_potcall_summary()
    append_timing("results.dat", t0)

See the user guide: https://eondocs.org/user_guide/pyeonclient.html
"""

from __future__ import annotations

import os
from contextlib import chdir
from pathlib import Path
from typing import Any, Sequence

from pyeonclient._core import (
    JobType,
    Matter,
    Parameters,
    Potential,
    PotType,
    RunStatus,
    append_results_timing,
    get_process_times,
    make_job,
    make_potential,
    steady_clock_now,
    write_potcall_summary,
)


def load_parameters(path: str | Path = "config.ini") -> Parameters:
    """Step: load ``config.ini`` (or JSON via Parameters.load_json)."""
    p = Parameters()
    p.load(str(path))
    return p


def write_neb_results(
    neb: Any,
    params: Parameters,
    force_calls_neb: int = 0,
) -> list[str]:
    """Step: write NEB artifacts (``results.dat``, ``neb.con``, ``sp.con``, …).

    Accepts either a chemist :class:`~pyeonclient.NEB` (uses ``.band`` after
    ``compute()``) or a raw :class:`~pyeonclient.NudgedElasticBand`.
    """
    from pyeonclient._core import neb_write_results

    band = getattr(neb, "band", neb)
    return list(neb_write_results(band, params, int(force_calls_neb)))


def run_job(params: Parameters) -> list[str]:
    """Step: ``make_job`` → ``Job.run`` → destroy Job (Potential released).

    Does not write potcall summary or timing; call those next.
    """
    job = make_job(params)
    try:
        return list(job.run())
    finally:
        # Explicit drop so PotRegistry sees destruction before summary.
        del job


def append_timing(
    results_path: str | Path = "results.dat",
    t0: float | None = None,
    *,
    elapsed: float | None = None,
) -> None:
    """Step: append ``time_seconds`` / ``user_time`` / ``system_time``.

    Lines use the results.dat ``<value> <key>`` contract (same as job writers
    and ``parse_results``). Pass *t0* from :func:`steady_clock_now` before the
    work, or *elapsed* wall seconds directly.
    """
    if elapsed is None:
        if t0 is None:
            raise ValueError("append_timing requires t0= or elapsed=")
        elapsed = steady_clock_now() - t0
    _real, user, system = get_process_times()
    append_results_timing(str(results_path), float(elapsed), float(user), float(system))


def write_minimization_results(
    params: Parameters,
    pot: Potential,
    matter: Matter,
    *,
    converged: bool,
    path: str | Path = "results.dat",
) -> None:
    """Step: write MinimizationJob-style results.dat body (no timing footer)."""
    status = RunStatus.GOOD if converged else RunStatus.FAIL_MAX_ITERATIONS
    # magic_enum-style names (same strings as MinimizationJob)
    pot_name = str(params.potential).rsplit(".", 1)[-1]
    status_name = str(status).rsplit(".", 1)[-1]
    status_i = 0 if status == RunStatus.GOOD else 1
    lines = [
        f"{status_i} termination_reason",
        f"{status_name} termination_reason_text",
        "minimization job_type",
        f"{pot_name} potential_type",
        f"{int(pot.force_call_counter)} total_force_calls",
    ]
    if status != RunStatus.FAIL_POTENTIAL_FAILED:
        lines.append(f"{matter.potential_energy} potential_energy")
    Path(path).write_text("\n".join(lines) + "\n")


def minimize_workdir(
    workdir: str | Path,
    *,
    pos: str = "pos.con",
    min_con: str = "min.con",
    movie_prefix: str = "minimization",
    write_movie: bool = True,
) -> list[str]:
    """Compose Matter steps for a minimization workdir (cookbook min_reactant).

    Steps: load_parameters → make_potential → Matter → con2matter → relax →
    matter2con → write_minimization_results → del handles → potcalls → timing.
    """
    work = Path(workdir).resolve()
    files: list[str] = []
    with chdir(work):
        t0 = steady_clock_now()
        params = load_parameters("config.ini")
        pot = make_potential(params.potential, params)
        matter = Matter(pot, params)
        from pyeonclient._core import io_ok

        st = matter.con2matter(pos)
        if not io_ok(st):
            raise RuntimeError(f"con2matter failed for {pos}: {st}")
        quiet = bool(getattr(params, "quiet", False))
        ckpt = bool(getattr(params, "checkpoint", False))
        _m, converged = matter.relax(
            inplace=True,
            quiet=quiet,
            write_movie=write_movie,
            checkpoint=ckpt,
            prefix_movie=movie_prefix,
            prefix_checkpoint="pos",
        )
        converged = bool(converged)
        st = matter.matter2con(min_con)
        if not io_ok(st):
            raise RuntimeError(f"matter2con failed for {min_con}: {st}")
        files.append(min_con)
        write_minimization_results(params, pot, matter, converged=converged)
        files.append("results.dat")
        del matter
        del pot
        write_potcall_summary("_potcalls.json")
        files.append("_potcalls.json")
        append_timing("results.dat", t0)
    return files


def run_job_in_directory(
    workdir: str | Path,
    params: Parameters | None = None,
) -> list[str]:
    """Compose ClientEON job steps in *workdir* (any job type, including NEB).

    Steps: (optional load config.ini) → run_job → write_potcall_summary →
    append_timing. Prefer :func:`minimize_workdir` when you want Matter.relax
    control rather than MinimizationJob.
    """
    work = Path(workdir).resolve()
    with chdir(work):
        t0 = steady_clock_now()
        if params is None:
            params = Parameters()
        ini = work / "config.ini"
        if ini.is_file():
            params.load(str(ini))
        files = run_job(params)
        write_potcall_summary("_potcalls.json")
        if "_potcalls.json" not in files:
            files.append("_potcalls.json")
        append_timing("results.dat", t0)
        return files


def run_eon_cwd() -> list[str]:
    """Cookbook helper: run whatever job config.ini requests in cwd."""
    return run_job_in_directory(Path("."), Parameters())


def neb_workdir(
    workdir: str | Path = ".",
    *,
    reactant: str = "reactant.con",
    product: str = "product.con",
) -> list[str]:
    """Compose NEB steps in *workdir* using :class:`NudgedElasticBand`.

    Steps (explicit)::

        load_parameters
        make_potential
        load endpoint Matter (or path-list endpoints when init=FILE)
        NudgedElasticBand(initial, final, params, pot)
        pot_registry_total_force_calls  # baseline
        neb.compute()
        find_extrema if GOOD
        neb_write_results
        del neb / pot
        write_potcall_summary
        append_timing
    """
    from pyeonclient._core import (
        NudgedElasticBand,
        NEBStatus,
        neb_write_results,
        pot_registry_total_force_calls,
    )

    work = Path(workdir).resolve()
    files: list[str] = []
    with chdir(work):
        t0 = steady_clock_now()
        params = load_parameters("config.ini")
        pot = make_potential(params.potential, params)

        initial = Matter(pot, params)
        final = Matter(pot, params)
        # Endpoints: path list first/last when FILE init, else reactant/product
        init_method = getattr(params, "neb_init_method", None)
        path_in = getattr(params, "neb_initial_path", "") or ""
        if path_in and Path(path_in).is_file():
            from pyeonclient._core import neb_read_file_paths

            plist = neb_read_file_paths(path_in)
            if len(plist) < 2:
                raise RuntimeError("NEB path list needs >= 2 frames")
            st0 = initial.con2matter(plist[0])
            st1 = final.con2matter(plist[-1])
        else:
            st0 = initial.con2matter(reactant)
            st1 = final.con2matter(product)
        from pyeonclient._core import io_ok

        if not io_ok(st0):
            raise RuntimeError("failed to load NEB reactant frame")
        if not io_ok(st1):
            raise RuntimeError("failed to load NEB product frame")

        if getattr(params, "neb_minimize_endpoints", False):
            initial.relax(
                inplace=True,
                quiet=bool(getattr(params, "quiet", False)),
                write_movie=bool(getattr(params, "write_movies", False)),
                checkpoint=bool(getattr(params, "checkpoint", False)),
                prefix_movie="react_neb",
                prefix_checkpoint="react_neb",
            )
            final.relax(
                inplace=True,
                quiet=bool(getattr(params, "quiet", False)),
                write_movie=bool(getattr(params, "write_movies", False)),
                checkpoint=bool(getattr(params, "checkpoint", False)),
                prefix_movie="prod_neb",
                prefix_checkpoint="prod_neb",
            )

        neb = NudgedElasticBand(initial, final, params, pot)
        f0 = pot_registry_total_force_calls()
        status = neb.compute()
        f_neb = pot_registry_total_force_calls() - f0
        if status == NEBStatus.GOOD:
            neb.find_extrema()
        files.extend(neb_write_results(neb, params, f_neb))
        del neb
        del initial
        del final
        del pot
        write_potcall_summary("_potcalls.json")
        files.append("_potcalls.json")
        append_timing("results.dat", t0)
    return files


def rgpot_metatomic_workdir(
    workdir: str | Path,
    *,
    engine_path: str | Path | None = None,
    model_path: str | Path | None = None,
    pos: str = "pos.con",
) -> list[str]:
    """Minimization via ``potential=RGPOT`` + ``backend=metatomic`` (dlopen engine).

    Thin-host path: host links Rgpot + loads ``libmetatomic_engine.so`` at
    runtime. Overwrites/extends config.ini fields if *engine_path* / *model_path*
    are given, then runs :func:`run_job_in_directory`.
    """
    work = Path(workdir).resolve()
    with chdir(work):
        params = load_parameters("config.ini") if Path("config.ini").is_file() else Parameters()
        params.potential = PotType.RGPOT
        params.job = JobType.Minimization
        params.rgpot_backend = "metatomic"
        if engine_path is not None:
            params.rgpot_engine_path = str(Path(engine_path).resolve())
        if model_path is not None:
            params.rgpot_model_path = str(Path(model_path).resolve())
            # also set Metatomic section for dual-read in RgpotPot
            params.metatomic_model_path = params.rgpot_model_path
        # Prefer Job path (RgpotPot) over Matter.relax (native Metatomic only)
        t0 = steady_clock_now()
        files = run_job(params)
        write_potcall_summary("_potcalls.json")
        if "_potcalls.json" not in files:
            files.append("_potcalls.json")
        append_timing("results.dat", t0)
        return files


def rgpot_metatomic_neb_workdir(
    workdir: str | Path,
    *,
    engine_path: str | Path | None = None,
    model_path: str | Path | None = None,
) -> list[str]:
    """NEB via RGPOT metatomic engine (same Job composition as eonclient)."""
    work = Path(workdir).resolve()
    with chdir(work):
        params = load_parameters("config.ini")
        params.potential = PotType.RGPOT
        params.rgpot_backend = "metatomic"
        if engine_path is not None:
            params.rgpot_engine_path = str(Path(engine_path).resolve())
        if model_path is not None:
            p = str(Path(model_path).resolve())
            params.rgpot_model_path = p
            params.metatomic_model_path = p
        t0 = steady_clock_now()
        files = run_job(params)
        write_potcall_summary("_potcalls.json")
        if "_potcalls.json" not in files:
            files.append("_potcalls.json")
        append_timing("results.dat", t0)
        return files
