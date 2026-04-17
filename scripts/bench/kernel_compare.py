"""Benchmark harness: compare amsel MCAMC kernels head-to-head.

Runs the three `sb_kernel` variants (fpta_bundle, mrm_direct, ngt)
through eon.mcamc on a deterministic superbasin fixture and prints a
timing + result table. The point is to let reviewers and potential
users see that amsel delivers the same answer across kernels and
quantify the per-call cost.

Invoke:

    python scripts/bench/kernel_compare.py [--states N] [--seed S]

Tracks bead amsel-yo8 Step 4.
"""
from __future__ import annotations

import argparse
import time

import numpy as np

import eon.mcamc as kernels


def build_deterministic_superbasin(
    n_transient: int, n_absorbing: int, seed: int
) -> tuple[list[int], list[int], list[tuple[int, int, float]]]:
    """Build (transient, absorbing, rate-edge-list) for a dense
    superbasin with log-uniform rates. Seed-deterministic."""
    rng = np.random.default_rng(seed)
    transient = list(range(n_transient))
    absorbing = list(range(n_transient, n_transient + n_absorbing))
    edges: list[tuple[int, int, float]] = []
    for i in transient:
        for j in transient:
            if i == j:
                continue
            if rng.random() < 0.5:
                rate = float(10.0 ** rng.uniform(-3.0, 3.0))
                edges.append((i, j, rate))
        for a in absorbing:
            if rng.random() < 0.5:
                rate = float(10.0 ** rng.uniform(-3.0, 3.0))
                edges.append((i, a, rate))
    # Guarantee entry (state 0) can reach SOMETHING.
    reaches_any_abs = any(f == 0 and t in absorbing for f, t, _ in edges)
    if not reaches_any_abs:
        edges.append((0, absorbing[0], 1.0))
    return transient, absorbing, edges


def build_q_r_c(
    transient: list[int],
    absorbing: list[int],
    rates: list[tuple[int, int, float]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_t, n_a = len(transient), len(absorbing)
    t_index = {sid: i for i, sid in enumerate(transient)}
    a_index = {sid: i for i, sid in enumerate(absorbing)}
    Q = np.zeros((n_t, n_t))
    R = np.zeros((n_t, n_a))
    c = np.zeros(n_t)
    for f, t, r in rates:
        if f not in t_index:
            continue
        c[t_index[f]] += r
        if t in t_index:
            Q[t_index[f], t_index[t]] += r
        else:
            R[t_index[f], a_index[t]] += r
    return Q, R, c


def bench_one(fn, repeats: int = 5):
    """Time a zero-arg callable. Returns (mean_ms, result)."""
    timings = []
    result = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = fn()
        t1 = time.perf_counter()
        timings.append((t1 - t0) * 1e3)
    return float(np.mean(timings)), result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transient", type=int, default=8)
    parser.add_argument("--absorbing", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260417)
    parser.add_argument("--repeats", type=int, default=50)
    args = parser.parse_args()

    transient, absorbing, edges = build_deterministic_superbasin(
        args.transient, args.absorbing, args.seed
    )
    print(
        f"benchmark: {args.transient} transient / {args.absorbing} absorbing, "
        f"{len(edges)} edges, seed {args.seed}"
    )

    # fpta_bundle (Q, R, c)
    Q, R, c = build_q_r_c(transient, absorbing, edges)
    ms_fb, (t_fb, B_fb, res_fb) = bench_one(
        lambda: kernels.mcamc(Q, R, c), repeats=args.repeats
    )
    mfpt_fb = float(t_fb[0])
    b_row_fb = B_fb[0, :]

    # mrm_direct
    ms_mrm, (tau_mrm, rates_mrm, _x) = bench_one(
        lambda: kernels.mrm_direct(
            transient=transient, absorbing=absorbing, rates=edges, entry=0
        ),
        repeats=args.repeats,
    )
    b_row_mrm = np.asarray(rates_mrm) * tau_mrm

    # ngt (one call per absorbing slot)
    def _ngt_row():
        committors = np.zeros(len(absorbing))
        mfpt = None
        for slot_idx, a in enumerate(absorbing):
            try:
                _rate, mfpt_i, committor = kernels.ngt(
                    transient=transient,
                    absorbing=absorbing,
                    rates=edges,
                    source=[0],
                    target=[a],
                )
                committors[slot_idx] = committor
                if mfpt is None:
                    mfpt = mfpt_i
            except Exception:
                continue
        return mfpt, committors

    ms_ngt, (mfpt_ngt, b_row_ngt) = bench_one(_ngt_row, repeats=args.repeats)

    print()
    print("kernel       |  mean ms |         MFPT | first slot P")
    print("-" * 58)
    print(f"fpta_bundle  | {ms_fb:8.3f} | {mfpt_fb:12.6e} | {b_row_fb[0]:12.6e}")
    print(f"mrm_direct   | {ms_mrm:8.3f} | {tau_mrm:12.6e} | {b_row_mrm[0]:12.6e}")
    print(f"ngt (per-slot)| {ms_ngt:8.3f} | {mfpt_ngt:12.6e} | {b_row_ngt[0]:12.6e}")
    print()
    # Sanity: the three MFPTs and absorption rows should agree.
    diffs = [
        ("mrm  vs fpta_bundle MFPT", abs(mfpt_fb - tau_mrm) / max(abs(mfpt_fb), 1e-12)),
        ("ngt  vs fpta_bundle MFPT", abs(mfpt_fb - mfpt_ngt) / max(abs(mfpt_fb), 1e-12)),
        ("mrm  vs fpta_bundle b[0]", abs(b_row_fb[0] - b_row_mrm[0]) / max(abs(b_row_fb[0]), 1e-12)),
        ("ngt  vs fpta_bundle b[0]", abs(b_row_fb[0] - b_row_ngt[0]) / max(abs(b_row_fb[0]), 1e-12)),
    ]
    print("cross-kernel relative differences:")
    for label, rel in diffs:
        print(f"  {label} = {rel:.3e}")


if __name__ == "__main__":
    main()
