#!/usr/bin/env python

# Author: Ian Johnson

import pathfix
import config
import sys
import os

try:
    os.makedirs('../docs/_autogen')
except:
    pass


class documentation():
    for i in range(len(config.format)):
        section = config.format[i].name.lower().replace(" ", "_")
        filename = '../docs/_autogen/config_%s.txt' % section
        file(filename, 'w')
        f =open(filename, 'w')
        f.write("%s\n" %config.format[i].name)
        f.write("-"*len(config.format[i].name))
        f.write("\n")
        f.write("\n")
        f.write("%s\n" %config.format[i].description)
        f.write("\n")

        for j in range(len(config.format[i].keys)):
            f.write("**%s**\n" %config.format[i].keys[j].name)
            f.write("\n")
            if config.format[i].keys[j].kind == "float":
                try:
                    f.write("   default: %g\n" %config.format[i].keys[j].default)
                except:
                    pass
            else:
                f.write("   default: %s\n" %config.format[i].keys[j].default)
            f.write("\n")
            f.write("   %s\n" %config.format[i].keys[j].description)
            f.write("\n")
            if config.format[i].keys[j].default != None:
                if len(config.format[i].keys[j].values) != 0:
                    f.write("   options:\n")
                    f.write("\n")
                    for k in range(len(config.format[i].keys[j].values)):
                        f.write("       **%s**: %s\n" %(config.format[i].keys[j].values[k].name, config.format[i].keys[j].values[k].description))
                        f.write("\n")
