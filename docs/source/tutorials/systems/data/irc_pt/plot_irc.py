"""Plot the energy along an IRC written by a saddle_search job with irc_after.

Usage: python plot_irc.py JOB_DIR OUTPUT.png

Reads irc.dat from JOB_DIR and, when present, the relaxed end energies from
results.dat.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_results(path):
    values = {}
    if path.is_file():
        for line in path.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2:
                values[parts[1]] = parts[0]
    return values


def main(job_dir, output):
    job_dir = Path(job_dir)
    data = np.loadtxt(job_dir / "irc.dat", skiprows=1)
    s, energy = data[:, 1], data[:, 2]
    minus, plus = s < 0, s > 0

    fig, ax = plt.subplots(figsize=(6, 3.6), constrained_layout=True)
    ax.plot(s[minus], energy[minus], ".-", ms=3, lw=1, color="tab:blue",
            label="against the mode (irc_min1)")
    ax.plot(s[plus], energy[plus], ".-", ms=3, lw=1, color="tab:orange",
            label="along the mode (irc_min2)")
    ax.plot([0.0], [0.0], marker="*", ms=14, color="gold", mec="black",
            ls="none", label="saddle")

    # Open circles: the ends after eOn's minimizer relaxed them.
    results = read_results(job_dir / "results.dat")
    if "potential_energy_saddle" in results:
        e0 = float(results["potential_energy_saddle"])
        ends = (("potential_energy_irc_min1", s[0], "tab:blue"),
                ("potential_energy_irc_min2", s[-1], "tab:orange"))
        for key, s_end, color in ends:
            if key in results:
                ax.plot([s_end], [float(results[key]) - e0], marker="o",
                        ms=8, mfc="none", color=color, ls="none")

    ax.axvline(0.0, color="0.7", lw=0.8, zorder=0)
    ax.set_xlabel(r"distance along the IRC from the saddle ($\AA$)")
    ax.set_ylabel("energy relative to the saddle (eV)")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(output, dpi=150)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
