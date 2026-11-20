#!/usr/bin/env python3
"""Render the CuH2 line-integral tables and their quoted numbers.

The SSoT is the CSV written by ``test_cuh2_line_integral`` into
``$TMPDIR/eon_cuh2_line_integral/`` and committed under
``docs/source/fig/data/``. Every number the user guide quotes about this
regression is computed here from those rows, so a table and its prose cannot
disagree. The Markdown fragments are generated at docs build time via
``docs/source/conf.py``; do not commit them.

Usage::

  python scripts/line_integral_cuh2_tables.py
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_DATA = _REPO_ROOT / "docs/source/fig/data"
_OUT = _REPO_ROOT / "docs/source/fig/generated"


def _rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as handle:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    if not rows:
        raise ValueError(f"{path} has no rows")
    return rows


def _pair(values: list[float], digits: int) -> str:
    return " and ".join(f"{value:.{digits}f}" for value in values)


def render_rates(rates: list[dict[str, float]]) -> str:
    lines = [
        "| {math}`T` (K) | {math}`\\ln k_{\\rm LI}-\\ln k_{\\rm eOn}` "
        "| estimated determinant error | {math}`k_{\\rm LI}/k_{\\rm eOn}` |",
        "|---:|---:|---:|---:|",
    ]
    for row in rates:
        difference = row["log_rate_li"] - row["log_rate_instanton"]
        lines.append(
            f"| {row['temperature_K']:.1f} | {difference:.3f} "
            f"| {row['log_rate_li_error']:.3f} "
            f"| {math.exp(difference):.3f} |"
        )
    return "\n".join(lines) + "\n"


def render_convergence(convergence: list[dict[str, float]]) -> str:
    lines = [
        "| n=10 {math}`T` (K) | n=10, P=150 | n=20, P=150 | n=40, P=150 "
        "| n=80, P=150 | n=80 {math}`T` (K) | n=80, P=600 "
        "| matched eOn instanton, P=600 |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in convergence:
        lines.append(
            f"| {row['temperature_n10_K']:.1f} "
            f"| {row['log_rate_n10_P150']:.4f} "
            f"| {row['log_rate_n20_P150']:.4f} "
            f"| {row['log_rate_n40_P150']:.4f} "
            f"| {row['log_rate_n80_P150']:.4f} "
            f"| {row['temperature_n80_K']:.1f} "
            f"| {row['log_rate_n80_P600']:.4f} "
            f"| {row['log_rate_instanton_P600']:.4f} |"
        )

    def change(column: str, base: str) -> list[float]:
        return [abs(row[column] - row[base]) for row in convergence]

    images = change("log_rate_n80_P150", "log_rate_n10_P150")
    shifts = change("temperature_n80_K", "temperature_n10_K")
    replicas10 = change("log_rate_n10_P600", "log_rate_n10_P150")
    replicas80 = change("log_rate_n80_P600", "log_rate_n80_P150")
    bars = [row["log_rate_n80_P600_error"] for row in convergence]
    ratios = [
        math.exp(row["log_rate_n80_P600"] - row["log_rate_instanton_P600"])
        for row in convergence
    ]
    prose = (
        f"Raising n from 10 to 80 changes {{math}}`\\ln k` by "
        f"{_pair(images, 3)} and shifts the temperatures by "
        f"{_pair(shifts, 2)} K at the two fixed energies. Raising P from 150 "
        f"to 600 changes {{math}}`\\ln k` by {_pair(replicas10, 3)} at n=10 "
        f"and by {_pair(replicas80, 3)} at n=80. With n=80 and P=600, "
        f"{{math}}`k_{{\\rm LI}}/k_{{\\rm eOn}}` is {_pair(ratios, 4)}, "
        f"and the stochastic log-determinants give diagnostic geometric-tail "
        f"error estimates of {_pair(bars, 4)} for those {{math}}`\\ln k`; "
        f"these are not certified bounds."
    )
    return "\n".join(lines) + "\n\n" + prose + "\n"


def render_calls(calls: list[dict[str, float]]) -> str:
    counts = [row["force_calls"] for row in calls]
    products = [row["hessian_vector_products"] for row in calls]
    continuation = counts[1:]
    return (
        f"The regression reaches {calls[-1]['temperature_K']:.1f} K after "
        f"{len(counts)} energy points, with 0.01 eV as the maximum energy "
        f"step and the figure's ±0.5 log-unit reading uncertainty converted "
        f"to its {{math}}`\\ln k` continuation tolerance. It uses the same "
        f"{{math}}`0.01` force tolerance and velocity-projection optimizer "
        f"stated in the 2018 paper. It requires {sum(counts) / len(counts):.0f}"
        f" energy/force calls per path on average ({sum(counts):.0f} in "
        f"total), compared with the paper's approximate 900. The initial "
        f"MEP-seeded path costs {counts[0]:.0f} calls; the remaining "
        f"{len(continuation)} continuation paths average "
        f"{sum(continuation) / len(continuation):.0f} calls. The P=150 HQTST "
        f"rates take {sum(products) / len(products):.0f} Hessian-vector "
        f"products per energy point on average ({sum(products):.0f} in "
        f"total, {min(products):.0f} to {max(products):.0f}).\n"
    )


def render_measurements(
    benchmark: list[dict[str, float]], calls: list[dict[str, float]]
) -> str:
    """Render cost figures whose values come only from committed CSV rows."""
    levels = len(benchmark)
    total_force_calls = benchmark[0]["total_force_calls"]
    total_products = benchmark[0]["total_hessian_vector_products"]
    time_seconds = benchmark[0]["time_seconds"]
    for row in benchmark:
        if (
            row["total_force_calls"] != total_force_calls
            or row["total_hessian_vector_products"] != total_products
            or row["time_seconds"] != time_seconds
        ):
            raise ValueError("seven-level benchmark summary columns disagree")

    lines = [
        f"The fresh one-core, one-thread benchmark completed {levels} levels "
        f"with {total_force_calls:.0f} force calls and {total_products:.0f} "
        f"Hessian-vector products in {time_seconds:.2f} s.",
        "",
        "| level | T (K) | ln(k/s^-1) | estimated determinant error "
        "| path force calls | Hessian-vector products |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in benchmark:
        path_calls = (
            row["force_calls_li_neb"] + row["force_calls_continuation"]
        )
        lines.append(
            f"| {row['level']:.0f} | {row['temperature_K']:.2f} "
            f"| {row['log_rate_per_s']:.4f} | {row['log_rate_error']:.4f} "
            f"| {path_calls:.0f} | {row['hessian_vector_products']:.0f} |"
        )

    products = [row["hessian_vector_products"] for row in calls]
    force_calls = [row["force_calls"] for row in calls]
    lines.extend(
        [
            "",
            f"The committed {len(calls)}-level CuH2 regression CSV records "
            f"{sum(force_calls):.0f} LI-NEB force calls and "
            f"{sum(products):.0f} Hessian-vector products. The per-level "
            f"product count ranges from {min(products):.0f} to "
            f"{max(products):.0f}, averaging "
            f"{sum(products) / len(products):.0f}.",
        ]
    )
    return "\n".join(lines) + "\n"


def generate_docs_assets(
    data: Path | None = None, out: Path | None = None
) -> list[Path]:
    """Write the three Markdown fragments and return their paths."""
    data = Path(data or _DATA)
    out = Path(out or _OUT)
    out.mkdir(parents=True, exist_ok=True)
    fragments = {
        "line_integral_cuh2_rates.md": render_rates(
            _rows(data / "cuh2_line_integral_rates.csv")
        ),
        "line_integral_cuh2_convergence.md": render_convergence(
            _rows(data / "cuh2_line_integral_convergence.csv")
        ),
        "line_integral_cuh2_calls.md": render_calls(
            _rows(data / "cuh2_line_integral_calls.csv")
        ),
        "line_integral_cuh2_measurements.md": render_measurements(
            _rows(data / "cuh2_line_integral_7level.csv"),
            _rows(data / "cuh2_line_integral_calls.csv"),
        ),
    }
    written = []
    for name, text in fragments.items():
        path = out / name
        path.write_text(text)
        written.append(path)
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=_DATA)
    parser.add_argument("--out", type=Path, default=_OUT)
    args = parser.parse_args(argv)
    for path in generate_docs_assets(args.data, args.out):
        print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
