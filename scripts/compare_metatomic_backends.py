#!/usr/bin/env python3
"""Compare metatomic force engines on one Matter geometry.

Backends (key → factory via pyeonclient.backends)::

  metatomic / metatomic_fat   — fat-linked PotType.METATOMIC
  rgpot_metatomic / metatomic_dlopen — RGPOT dlopen libmetatomic_engine.so
  ase_metatomic               — metatomic_ase.MetatomicCalculator via ASE wrap

Each backend runs in a **subprocess**. Fat-linked C++ and the Python
``metatomic-torch`` extension both register ``TORCH_LIBRARY(metatomic)``; loading
both in one process aborts.

Env::

  EON_PET_MAD_MODEL   path to .pt (required)
  RGPOT_METATOMIC_ENGINE  path to libmetatomic_engine.so (for dlopen)
  EON_PET_MAD_POS     optional pos.con (else cookbook default / d016)

Usage::

  python scripts/compare_metatomic_backends.py
  python scripts/compare_metatomic_backends.py --json /tmp/mta-compare.json
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


def _prepare_worker_env(backend: str, base_env: dict | None = None) -> dict:
    """Per-backend env so ASE Python metatomic and fat C++ do not double-load.

    Fat-linked ``libmetatomic_pot`` and pip ``metatomic-torch`` each register
    ``TORCH_LIBRARY(metatomic)``.  ASE workers therefore use a non-metatomic
    ``_core`` (``build-pyeon-methods``) when available.
    """
    import shutil
    import tempfile

    env = dict(base_env or os.environ)
    repo = Path(__file__).resolve().parents[1]
    root_candidates = [
        Path.cwd(),
        repo,
    ]
    if backend != "ase_metatomic":
        return env

    core = os.environ.get("EON_ASE_PYEON_CORE", "")
    core_path = Path(core) if core else None
    if core_path is None or not core_path.is_file():
        for root in root_candidates:
            for rel in (
                "bbdir-mta-bench-ase/client/python/_core.abi3.so",
                "build-pyeon-methods/client/python/_core.abi3.so",
                "bbdir-mta-bench/client/python/_core.abi3.so",
            ):
                cand = root / rel
                if cand.is_file():
                    # Prefer non-metatomic cores when discoverable via name.
                    core_path = cand
                    if "ase" in rel or "methods" in rel:
                        break
            if core_path is not None and core_path.is_file() and (
                "ase" in str(core_path) or "methods" in str(core_path)
            ):
                break
    if core_path is None or not core_path.is_file():
        # fall back to default import; may TORCH_LIBRARY-clash with fat-linked core
        return env

    pkg_src = None
    for root in root_candidates:
        src = root / "client" / "python" / "pyeonclient"
        if (src / "backends.py").is_file():
            pkg_src = src
            break
    if pkg_src is None:
        return env

    tmp = Path(tempfile.mkdtemp(prefix="pyec-ase-"))
    dest = tmp / "pyeonclient"
    dest.mkdir()
    for f in pkg_src.glob("*.py"):
        shutil.copy2(f, dest / f.name)
    shutil.copy2(core_path, dest / "_core.abi3.so")
    # client build dir for sibling .so via rpath $ORIGIN/..
    build_client = core_path.parent.parent  # .../client/python -> .../client
    ld = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = f"{dest}:{build_client}:{ld}"
    pp = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{tmp}:{pp}"
    env["EON_ASE_WORKER_TMP"] = str(tmp)
    return env



def _load_geometry(pos_con: Path):
    """Return (positions, Z int64, cell, masses) from a .con via ASE."""
    try:
        from ase.io import read

        atoms = read(str(pos_con))
        pos = np.ascontiguousarray(atoms.get_positions(), dtype=np.float64)
        z = np.ascontiguousarray(atoms.get_atomic_numbers(), dtype=np.int64)
        cell = np.ascontiguousarray(np.asarray(atoms.get_cell()), dtype=np.float64)
        masses = np.ascontiguousarray(atoms.get_masses(), dtype=np.float64)
        return pos, z, cell, masses, bool(np.any(atoms.pbc))
    except Exception as e:
        raise SystemExit(f"failed to read {pos_con}: {e}") from e


def _eval(name: str, pot, pos, z, cell, n_repeat: int = 3) -> dict:
    out: dict = {"backend": name, "ok": False}
    try:
        # warm-up
        E0, F0 = pot.get_ef(pos, z, cell)
        t0 = time.perf_counter()
        E, F = E0, F0
        for _ in range(n_repeat):
            E, F = pot.get_ef(pos, z, cell)
        elapsed = time.perf_counter() - t0
        F = np.asarray(F, dtype=np.float64)
        out.update(
            {
                "ok": True,
                "energy": float(E),
                "max_force": float(np.max(np.linalg.norm(F, axis=1))),
                "force_rms": float(np.sqrt(np.mean(F * F))),
                "n_atoms": int(pos.shape[0]),
                "seconds_per_call": elapsed / n_repeat,
                "n_repeat": n_repeat,
            }
        )
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def _worker_main(argv: list[str] | None = None) -> int:
    """Evaluate one backend; print JSON row to stdout."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--pos", required=True)
    ap.add_argument("--engine", default="")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--repeat", type=int, default=3)
    args = ap.parse_args(argv)

    pos, z, cell, _masses, _pbc = _load_geometry(Path(args.pos))
    from pyeonclient.backends import make_backend

    try:
        kw: dict = {"model_path": args.model, "device": args.device}
        if args.backend in ("rgpot_metatomic", "metatomic_dlopen"):
            kw["engine_path"] = args.engine
        pot = make_backend(args.backend, **kw)
        row = _eval(args.backend, pot, pos, z, cell, n_repeat=args.repeat)
    except Exception as e:
        row = {
            "backend": args.backend,
            "ok": False,
            "error": f"{type(e).__name__}: {e}",
        }
    sys.stdout.write(json.dumps(row) + "\n")
    return 0 if row.get("ok") else 1


def _resolve_paths(args: argparse.Namespace):
    repo = Path(__file__).resolve().parents[1]
    model = Path(args.model).expanduser() if args.model else None
    if not model or not model.is_file():
        for cand in (
            repo / "subprojects/gpr_optim/bench_data/petmad/pet-mad-s-v1.5.0.pt",
            repo.parent / "gpr_optim/bench_data/petmad/pet-mad-s-v1.5.0.pt",
            Path.cwd() / "subprojects/gpr_optim/bench_data/petmad/pet-mad-s-v1.5.0.pt",
        ):
            if cand.is_file():
                model = cand
                break
    if not model or not model.is_file():
        print("error: set --model / EON_PET_MAD_MODEL", file=sys.stderr)
        return None, None, None

    pos_con = Path(args.pos) if args.pos else None
    if not pos_con or not pos_con.is_file():
        for cand in (
            model.parent / "d016_pos.con",
            model.parent / "minimal_pt.con",
            repo / "subprojects/gpr_optim/bench_data/petmad/d016_pos.con",
            repo / "docs/lj13.con",
        ):
            if cand.is_file():
                pos_con = cand
                break
    if not pos_con or not pos_con.is_file():
        print("error: set --pos / EON_PET_MAD_POS", file=sys.stderr)
        return None, None, None

    engine = Path(args.engine) if args.engine else None
    if not engine or not engine.is_file():
        for cand in (
            repo / "bbdir-mta-bench/client/librgpot_metatomic_engine.so",
            repo / "bbdir-mta-bench/client/libmetatomic_engine.so",
            repo / "bbdir/client/librgpot_metatomic_engine.so",
            repo / "bbdir/client/libmetatomic_engine.so",
            Path("bbdir-mta-bench/client/librgpot_metatomic_engine.so"),
            Path("bbdir-mta-bench/client/libmetatomic_engine.so"),
            Path("build-pyeon-mta-rgpot/client/librgpot_metatomic_engine.so"),
            Path("build-pyeon-mta-rgpot/client/libmetatomic_engine.so"),
            Path("build-mta-compare/client/librgpot_metatomic_engine.so"),
            Path("build-mta-compare/client/libmetatomic_engine.so"),
            Path("build-pyeon-mta/client/librgpot_metatomic_engine.so"),
            Path("build-pyeon-mta/client/libmetatomic_engine.so"),
        ):
            if cand.is_file():
                engine = cand
                break
    return model, pos_con, engine


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if argv and argv[0] == "--worker":
        return _worker_main(argv[1:])

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--model",
        default=os.environ.get("EON_PET_MAD_MODEL", ""),
        help="metatomic .pt model path",
    )
    ap.add_argument(
        "--engine",
        default=os.environ.get("RGPOT_METATOMIC_ENGINE", ""),
        help="libmetatomic_engine.so for dlopen path",
    )
    ap.add_argument(
        "--pos",
        default=os.environ.get("EON_PET_MAD_POS", ""),
        help="pos.con geometry",
    )
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--repeat", type=int, default=3)
    ap.add_argument("--json", dest="json_out", default="")
    args = ap.parse_args(argv)

    model, pos_con, engine = _resolve_paths(args)
    if model is None:
        return 2

    import pyeonclient as pyec
    from pyeonclient.backends import list_backends

    print("pyeonclient", pyec.__version__)
    print("built_with_metatomic", pyec.built_with_metatomic())
    print("built_with_rgpot", pyec.built_with_rgpot())
    print("backends", list_backends())
    print("model", model)
    print("pos", pos_con)
    print("engine", engine)

    pos, z, cell, masses, pbc = _load_geometry(pos_con)
    print(f"geometry n={len(z)} pbc={pbc}")

    jobs: list[tuple[str, dict]] = [
        ("metatomic", {"model_path": str(model), "device": args.device}),
        (
            "rgpot_metatomic",
            {
                "model_path": str(model),
                "engine_path": str(engine) if engine else "",
                "device": args.device,
            },
        ),
        ("ase_metatomic", {"model_path": str(model), "device": args.device}),
    ]

    results = []
    script = str(Path(__file__).resolve())
    for key, kw in jobs:
        print(f"\n=== {key} ===")
        cmd = [
            sys.executable,
            script,
            "--worker",
            "--backend",
            key,
            "--model",
            str(model),
            "--pos",
            str(pos_con),
            "--device",
            args.device,
            "--repeat",
            str(args.repeat),
        ]
        if engine:
            cmd.extend(["--engine", str(engine)])
        proc = None
        try:
            proc = subprocess.run(
                cmd,
                check=False,
                capture_output=True,
                text=True,
                env=_prepare_worker_env(key),
            )
            # last non-empty stdout line is JSON row
            lines = [ln for ln in (proc.stdout or "").splitlines() if ln.strip()]
            if not lines:
                err = ((proc.stderr or "") + "\n" + (proc.stdout or "")).strip()[:800]
                row = {
                    "backend": key,
                    "ok": False,
                    "error": f"worker exit {proc.returncode}: {err}",
                }
            else:
                row = json.loads(lines[-1])
        except Exception as e:
            row = {"backend": key, "ok": False, "error": f"{type(e).__name__}: {e}"}
        results.append(row)
        if row.get("ok"):
            print(
                f"  E={row['energy']:.8f}  "
                f"max|F|={row['max_force']:.6e}  "
                f"rmsF={row['force_rms']:.6e}  "
                f"t={row['seconds_per_call']*1e3:.2f} ms/call"
            )
        else:
            print(f"  FAIL {row.get('error')}")
            if proc is not None and proc.stderr:
                tail = proc.stderr.strip().splitlines()[-8:]
                for line in tail:
                    print(f"    | {line}")

    # pairwise energy deltas among successes
    ok = [r for r in results if r.get("ok")]
    if len(ok) >= 2:
        print("\n=== energy deltas (row - col) ===")
        for i, a in enumerate(ok):
            for j, b in enumerate(ok):
                if j <= i:
                    continue
                dE = a["energy"] - b["energy"]
                dF = abs(a["max_force"] - b["max_force"])
                print(
                    f"  {a['backend']} - {b['backend']}: "
                    f"dE={dE:.6e}  d(max|F|)={dF:.6e}"
                )

    payload = {
        "model": str(model),
        "pos": str(pos_con),
        "engine": str(engine) if engine else None,
        "device": args.device,
        "built_with_metatomic": bool(pyec.built_with_metatomic()),
        "built_with_rgpot": bool(pyec.built_with_rgpot()),
        "results": results,
    }
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(payload, indent=2) + "\n")
        print("wrote", args.json_out)
    return 0 if all(r.get("ok") for r in results) else (0 if any(r.get("ok") for r in results) else 1)


if __name__ == "__main__":
    raise SystemExit(main())
