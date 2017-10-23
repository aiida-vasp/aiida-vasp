# -*- coding: utf-8 -*-
"""
setup: usage: pip install -e .[graphs]
"""

import re
from setuptools import setup, find_packages

# Get the version number
with open('./aiida_vasp/__init__.py') as f:
    MATCH_EXPR = "__version__[^'\"]+(['\"])([^'\"]+)"
    VERSION = re.search(MATCH_EXPR, f.read()).group(2).strip()

if __name__ == '__main__':
    setup(
        name='aiida-vasp',
        version=VERSION,
        description='AiiDA Plugin for running VASP calculations.',
        url='https://github.com/DropD/aiida-vasp',
        author='Rico Häuselmann',
        author_email='haeuselm@epfl.ch',
        license='MIT License, see LICENSE.txt file.',
        classifiers=[
            'Development Status :: 3 - Alpha', 'Environment :: Plugins',
            'Framework :: AiiDA', 'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 2.7',
            'Topic :: Scientific/Engineering :: Physics'
        ],
        keywords='vasp aiida wannier90 workflows',
        packages=find_packages(exclude=['aiida']),
        include_package_data=True,
        setup_requires=['reentry'],
        reentry_register=True,
        install_requires=[
            'aiida-core[atomic_tools]', 'ase', 'pymatgen', 'subprocess32',
            'click', 'chainmap'
        ],
        extras_require={
            'graphs': ['matplotlib'],
            'regen_default_paws': ['lxml'],
            'test': ['aiida-pytest'],
            'wannier': ['aiida-wannier90'],
            'dev': [
                'pre-commit', 'prospector', 'pylint', 'flake8', 'yapf',
                'coverage', 'pytest', 'pytest-cov', 'pgtest >= 1.1.0'
            ]
        },
        entry_points={
            'aiida.calculations': [
                'vasp.vasp = aiida_vasp.calcs.vasp:VaspCalculation',
                'vasp.vasp2w90 = aiida_vasp.calcs.vasp2w90:Vasp2w90Calculation'
            ],
            'aiida.data': [
                'vasp.archive = aiida_vasp.data.archive:ArchiveData',
                'vasp.chargedensity = aiida_vasp.data.chargedensity:ChargedensityData',
                'vasp.paw = aiida_vasp.data.paw:PawData',
                'vasp.wavefun = aiida_vasp.data.wavefun:WavefunData'
            ],
            'aiida.parsers': [
                'vasp.vasp = aiida_vasp.parsers.vasp:VaspParser',
                'vasp.vasp2w90 = aiida_vasp.parsers.vasp2w90:Vasp2w90Parser'
            ],
            'aiida.workflows': [
                'vasp.scf = aiida_vasp.workflows.scf:ScfWorkflow',
                'vasp.nscf = aiida_vasp.workflows.nscf:NscfWorkflow',
                'vasp.projections = aiida_vasp.workflows.projections:ProjectionsWorkflow',
                'vasp.autowindows = aiida_vasp.workflows.autowindows:AutowindowsWorkflow',
                'vasp.wannier = aiida_vasp.workflows.wannier:WannierWorkflow',
                'vasp.windows= aiida_vasp.workflows.windows:WindowsWorkflow',
            ],
            'aiida.cmdline.data': ['vasp.paw = aiida_vasp.commands.paw:_Paw']
        },
        scripts=['utils/runwf.py'])
