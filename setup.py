# -*- coding: utf-8 -*-
"""
setup: usage: pip install -e .[graphs]
"""

import json
from setuptools import setup, find_packages

if __name__ == '__main__':
    with open('setup.json', 'r') as info:
        SETUP_KWARGS = json.load(info)
    setup(packages=find_packages(exclude=['aiida']), **SETUP_KWARGS)
