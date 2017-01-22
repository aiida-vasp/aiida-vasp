# -*- coding: utf-8 -*-
"""
setup: usage: pip install -e .[graphs]
"""

from setuptools import setup, find_packages

if __name__ == '__main__':
    setup(
        name='aiida_vasp',
        version='0.9.0a1',
        description='AiiDA Plugin for running VASP -> Wannier workflows',
        url='https://github.com/DropD/aiida-vasp',
        author='Rico Häuselmann',
        author_email='haeuselm@epfl.ch',
        license='MIT License, see LICENSE.txt file.',
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Plugins',
            'Framework :: AiiDA',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 2.7',
            'Topic :: Scientific/Engineering :: Physics'
        ],
        keywords='vasp aiida wannier90 workflows',
        packages=find_packages(),
        include_package_data=True,
        install_requires=[
            'aiida[ssh]',
            'ase',
            'pymatgen',
            'subprocess32'
        ],
        extras_require={
            'graphs': ['matplotlib']
            'regen_default_paws': ['lxml']
        },
        entry_points={
            'aiida.calculations': [
                'vasp.scf = aiida_vasp.calcs.scf:ScfCalculation',
                'vasp.nscf = aiida_vasp.calcs.nscf:NscfCalculation',
                'vasp.amn = aiida_vasp.calcs.amn:AmnCalculation',
                'vasp.wannier = aiida_vasp.calcs.wannier:WannierCalculation',
                'vasp.vasp5 = aiida_vasp.calcs.vasp5:Vasp5Calculation'
            ],
            'aiida.data': [
                'vasp.archive = aiida_vasp.data.archive:ArchiveData',
                'vasp.chargedensity = aiida_vasp.data.chargedensity:ChargedensityData',
                'vasp.paw = aiida_vasp.data.paw:PawData',
                'vasp.wavefun = aiida_vasp.data.wavefun:WavefunData'
            ],
            'aiida.parsers': [
                'vasp.scf = aiida_vasp.parsers.scf:ScfParser',
                'vasp.nscf = aiida_vasp.parsers.nscf:NscfParser',
                'vasp.amn = aiida_vasp.parsers.amn:AmnParsaer',
                'vasp.wannier = aiida_vasp.parsers.wannier:WannierParser',
                'vasp.vasp5 = aiida_vasp.parsers.vasp5:Vasp5Parser'
            ],
            'aiida.workflows': [
                'vasp.scf = aiida_vasp.workflows.scf:ScfWorkflow',
                'vasp.nscf = aiida_vasp.workflows.nscf:NscfWorkflow',
                'vasp.projections = aiida_vasp.workflows.projections:ProjectionsWorkflow',
                'vasp.autowindows = aiida_vasp.workflows.autowindows:AutowindowsWorkflow',
                'vasp.wannier = aiida_vasp.workflows.wannier:WannierWorkflow',
                'vasp.windows= aiida_vasp.workflows.windows:WindowsWorkflow',
            ],
            'aiida.cmdline': [
                'vasp.paw = aiida_vasp.commands.paw:_Paw'
            ]
        },
        scripts=['utils/runwf.py']
    )
