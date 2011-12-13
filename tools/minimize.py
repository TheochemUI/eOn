#!/usr/bin/env python

import os
import sys
import shutil
import tempfile
import ConfigParser
import optparse

import pathfix

if __name__ == "__main__":

    op = optparse.OptionParser(usage = "%prog [options] <input con file> <potential> <output con file>")
    op.add_option("--box", action="store_true", dest="box", default=False,
                  help="relax the box along with the atomic coordinates")
    (options, args) = op.parse_args()
    if len(args) < 3:
        op.print_help()
        sys.exit()

    config = ConfigParser.SafeConfigParser()
    config.read(os.path.join(pathfix.path, "default_config.ini"))

    cwd = os.getcwd()
    td = tempfile.mkdtemp()

    shutil.copyfile(args[0], os.path.join(td, "reactant_passed.con"))

    config.set("Main", "potential", args[1])
    config.set("Main", "job", "minimization")
    if options.box:
        config.add_section('Optimizers')
        config.set("Optimizers", "opt_method", "box")

    cf = open(os.path.join(td, "config_passed.ini"), 'w')
    config.write(cf)
    cf.close()

    os.chdir(td)
    os.system(os.path.join(pathfix.path, "client", "client"))
    shutil.copyfile(os.path.join(td, "reactant.con"), os.path.join(cwd, args[2]))

