# -*- coding: utf-8 -*-
"""
setup: usage: pip install -e .[graphs]
"""

import os
import json
from setuptools import setup, find_packages

SETUP_JSON_PATH = os.path.join(os.path.dirname(__file__), 'setup.json')

if __name__ == '__main__':
    with open(SETUP_JSON_PATH, 'r') as info:
        SETUP_KWARGS = json.load(info)
    setup(packages=find_packages(exclude=['aiida']), **SETUP_KWARGS)
