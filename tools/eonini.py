#!/usr/bin/env python

import os
import ConfigParser

import pathfix
import yaml

if os.path.isfile('config.ini'):
	print 'Error: config.ini file already exists. Aborting.'
	import sys
	sys.exit()

yaml_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'eon', 'config.yaml')
yaml_file = open(yaml_path)
y = yaml.load(yaml_file)
yaml_file.close()

f = open('config.ini', 'w')

for section in y:
	f.write('\n') 
	f.write("[%s]\n" % section)
	for option in y[section]['options']:
		f.write("%s = %s\n" % (option, str(y[section]['options'][option]['default'])))

f.close()

