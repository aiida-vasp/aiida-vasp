# -*- coding: utf-8 -*-
"""
setup: usage: pip install -e .[graphs]
"""

from setuptools import setup, find_packages

if __name__ == '__main__':
    setup(
        name='aiida-vasp',
        version='0.9.0a1',
        description='AiiDA Plugin for running VASP -> Wannier workflows',
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
            'aiida-core', 'ase', 'pymatgen', 'subprocess32', 'click',
            'chainmap'
        ],
        extras_require={
            'graphs': ['matplotlib'],
            'regen_default_paws': ['lxml'],
            'test': ['aiida-pytest'],
            'wannier': ['aiida-wannier90']
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
