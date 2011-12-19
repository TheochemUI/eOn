#!/usr/bin/env python
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

packages = ['eon']
package_dir = {'eon':'eon'}
package_data = {'eon': ['config.yaml']}

setup(name='eon',
      packages=packages,
      package_dir=package_dir,
      package_data=package_data,
      scripts=['bin/eon'])
