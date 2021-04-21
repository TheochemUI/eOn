#!/usr/bin/env python
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

packages = ['eon']
package_dir = {'eon':'eon'}
package_data = {'eon': ['config.yaml']}

setup(name='eon',
    version='1.0.0',
    author='Henkelman Group',
    url='http://henkelmanlab/eon/',
    packages=packages,
    package_dir=package_dir,
    package_data=package_data,
    entry_points={'console_scripts': ['eon=eon.main:main']})
