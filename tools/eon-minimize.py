#!/usr/bin/env python

import os
import sys
import shutil
import tempfile
import configparser
import optparse

import pathfix

if __name__ == "__main__":

    op = optparse.OptionParser(usage = "%prog [options] <input con file> <potential> <output con file>")
    op.add_option("--box", action="store_true", dest="box", default=False,
                  help="relax the box along with the atomic coordinates")
    (options, args) = op.parse_args()
#    if len(args) < 3:
#        op.print_help()
#        sys.exit()

# get the input file
    posfile = "pos.con"
    if len(args) > 0:
        posfile = args[0]

    config = configparser.SafeConfigParser()
    if os.path.isfile("config.ini"):
        isconfig = True
        config.read("config.ini")
    else:
        config.read(os.path.join(pathfix.path, "default_config.ini"))

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

    shutil.copyfile(posfile, os.path.join(td, "pos.con"))

    if os.path.exists("potfiles"):
        potfiles = os.listdir("potfiles")
        for potfile in potfiles:
            if os.path.isfile(os.path.join("potfiles", potfile)):
                shutil.copyfile(os.path.join("potfiles", potfile), os.path.join(td, potfile))

    config.set("Main", "job", "minimization")
    config.set("Potential", "potential", potential)
    if options.box:
        config.add_section('Optimizers')
        config.set("Optimizers", "opt_method", "box")

    cf = open(os.path.join(td, "config.ini"), 'w')
    config.write(cf)
    cf.close()

    os.chdir(td)
#    os.system(os.path.join(pathfix.path, "../client", "eonclient"))
    os.system("eonclient")
    shutil.copyfile(os.path.join(td, "min.con"), os.path.join(cwd, "min.con"))
