from pypotlib.systems import cu_slab_h2 as cuh2slab
import ase
import ase.io
from ase.mep import NEB, NEBTools
from ase.optimize import BFGS, MDMin, FIRE, LBFGS
import matplotlib.pyplot as plt
import torch

import itertools as it
import copy
import os
import sys
import subprocess
from pathlib import Path
import configparser

from pyprotochemgp.transforms._ase import ASE_MoveCoordTrain
from pyprotochemgp.systems.cuh2slab import prepare_scatter_points


# Kanged from rgoswami.me
def getstrform(pathobj):
    return str(pathobj.absolute())


gitroot = Path(
    subprocess.run(
        ["git", "rev-parse", "--show-toplevel"], check=True, capture_output=True
    )
    .stdout.decode("utf-8")
    .strip()
)

cuh2_min1 = ase.io.read("min1_0Cu.con")
cuh2_min1.calc = cuh2slab.CuH2PotSlab()
true_e_dat = cuh2slab.plt_data(
    cuh2_min1,
    hh_range=cuh2slab.PltRange(low=-0.05, high=5),
    h2slab_range=cuh2slab.PltRange(low=-0.05, high=5),
    n_points=cuh2slab.PlotPoints(x_npt=40, y_npt=40)
)

# print([x.energy for x in true_e_dat.pltpts])

def plot_band(_index, _band, _k, _method="EON", _opt="QM", _ci="False"):
    plot_last = [ASE_MoveCoordTrain().transform(x) for x in _band]
    cuh2slab.contour_plot(
        true_e_dat.pltpts,
        scatter_points=prepare_scatter_points(torch.tensor(plot_last), cuh2_min1),
        title=f"({_opt}, CI: {_ci}, {_k})\n Band {_index} ({_method})",
    )
    oname = f"{_method}_{_opt}_{_ci}"
    plt.savefig(f"neb_path_{oname}_{_index:04d}.png")
    plt.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Script to convert .con files to .png files for CuH2"
    )
    parser.add_argument(
        "--ffmpeg",
        help="Generate video via subprocess",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    eon_conf = configparser.ConfigParser()
    eon_conf.read("config.ini")
    ec_k = eon_conf["Nudged Elastic Band"]["spring"]
    ec_opt = eon_conf["Optimizer"]["opt_method"]
    ec_ci = eon_conf["Nudged Elastic Band"]["climbing_image_method"]

    cwd = Path(".")
    for idx, con_file in enumerate(cwd.glob("neb_path_*.con")):
        plot_band(
            idx,
            ase.io.eon.read_eon(con_file),
            _k=ec_k,
            _opt=ec_opt,
            _ci=ec_ci,
        )

    if args.ffmpeg:
        oname = f"EON_{ec_opt}_{ec_ci}"
        subprocess.run(
            [
                "ffmpeg",
                "-framerate",
                "1",
                "-pattern_type",
                "glob",
                "-i",
                f"neb_path_{oname}_*.png",
                "-c:v",
                "libx264",
                "-r",
                "30",
                "-pix_fmt",
                "yuv420p",
                f"neb_path_{oname}.mp4",
            ]
        )
