#!/usr/bin/env python

import os
import configparser

import pathfix
import yaml

yaml_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'eon', 'config.yaml')
yaml_file = open(yaml_path)
y = yaml.load(yaml_file)
yaml_file.close()

for section in y:
    print()
    print("[%s]" % section)
    for option in y[section]['options']:
	    print("%s = %s" % (option, str(y[section]['options'][option]['default'])))
