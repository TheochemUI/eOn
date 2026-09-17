#!/usr/bin/env python

import os
import sys
import shutil
import tempfile
import configparser
import optparse
from pathlib import Path

import pathfix

if __name__ == "__main__":

    op = optparse.OptionParser(usage = "%prog [options] <input con file> <potential> <output con file>")
    op.add_option("--box", action="store_true", dest="box", default=False,
                  help="relax the box along with the atomic coordinates")
    (options, args) = op.parse_args()

# get the input file
    posfile = "pos.con"
    if len(args) > 0:
        posfile = args[0]

    config = configparser.ConfigParser()
    if Path("config.ini").is_file():
        isconfig = True
        config.read("config.ini")
    else:
        config.read(Path(pathfix.path) / "default_config.ini")

# set the potential
    if isconfig:
        potential = config.get("Potential", "potential", "none")
    if len(args) > 1:
        potential = args[1]
    if potential == "none":
        op.print_help()
        sys.exit()

    cwd = os.getcwd()
    td = tempfile.mkdtemp()

    shutil.copyfile(posfile, Path(td) / "pos.con")

    if Path("potfiles").exists():
        potfiles = os.listdir("potfiles")
        for potfile in potfiles:
            src = Path("potfiles") / potfile
            if src.is_file():
                shutil.copyfile(src, Path(td) / potfile)

    config.set("Main", "job", "minimization")
    config.set("Potential", "potential", potential)
    if options.box:
        config.add_section('Optimizers')
        config.set("Optimizers", "opt_method", "box")

    cf = open(Path(td) / "config.ini", 'w')
    config.write(cf)
    cf.close()

    os.chdir(td)
    os.system("eonclient")
    shutil.copyfile(Path(td) / "min.con", Path(cwd) / "min.con")
