# -*- coding: utf-8 -*-
"""
Install the aiida_vasp python package.

usage: pip install -e .[graphs]
"""

import os
import json
from setuptools import setup, find_packages

SETUP_JSON_PATH = os.path.join(os.path.dirname(__file__), 'setup.json')
README_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'README.rst')
with open(README_PATH, 'r') as readme:
    LONG_DESCRIPTION = readme.read()

if __name__ == '__main__':
    with open(SETUP_JSON_PATH, 'r') as info:
        SETUP_KWARGS = json.load(info)
    setup(packages=find_packages(exclude=['aiida']),
          keywords='vasp, aiida, wannier90, workflow, materials',
          long_description=LONG_DESCRIPTION,
          **SETUP_KWARGS)
