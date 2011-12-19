#!/usr/bin/env python
#try:
#    from setuptools import setup
#except ImportError:
from distutils.core import setup

packages = ['eonpy']
package_dir = {'eonpy':'eonpy'}
package_data = {'eonpy': ['config.yaml']}

setup(name='eonpy',
      packages=packages,
      package_dir=package_dir,
      package_data=package_data,
      scripts=['bin/eon'])
