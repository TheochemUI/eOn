#!/usr/bin/env python

import configparser
from pathlib import Path

import pathfix
import yaml

yaml_path = Path(__file__).resolve().parent.parent / "eon" / "config.yaml"
y = yaml.load(yaml_path.read_text())

for section in y:
    print()
    print("[%s]" % section)
    for option in y[section]['options']:
	    print("%s = %s" % (option, str(y[section]['options'][option]['default'])))
