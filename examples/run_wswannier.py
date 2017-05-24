#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>

from __future__ import division, print_function, unicode_literals

import itertools

from aiida.orm import Code, CalculationFactory, Computer, QueryBuilder
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.vasp.archive import ArchiveData

def get_input_archive():
    archive_description = u'InSb unstrained Wannier90 input'
    qb = QueryBuilder()
    qb.append(
        ArchiveData,
        filters={'description': {'==': archive_description}}
    )
    res = qb.all()
    if len(res) == 0:
        # create archive
        res = ArchiveData()
        res.add_file('reference_input/wannier90.mmn')
        res.add_file('reference_input/wannier90.amn')
        res.add_file('reference_input/wannier90.eig')
        res.description = archive_description
        res.store()
    elif len(res) > 1:
        raise ValueError('Query returned more than one matching ArchiveData instance.')
    else:
        res = res[0][0]
    return res
    
def run_wswannier():
    input_archive = get_input_archive()
    code = Code.get_from_string('Wannier90_2.1.0')
    calc = CalculationFactory('vasp.wswannier')()
    calc.use_code(code)
    # Monch 
    calc.set_resources(dict(
        num_machines=1,
        tot_num_mpiprocs=1
    ))
    calc.set_computer(Computer.get('Monch'))
    calc.set_queue_name('express_compute')
    calc.use_data(input_archive)
    # hard-code parameters for InSb test
    calc.use_settings(ParameterData(dict=dict(
        num_wann=14,
        num_bands=36,
        dis_num_iter=1000,
        num_iter=0,
        dis_win_min=-4.5,
        dis_win_max=16.,
        dis_froz_min=-4.,
        dis_froz_max=6.5,
        write_hr=True,
        use_ws_distance=True,
        write_xyz=True,
        write_tb=True,
        projections=[
            ['In : s; px; py; pz'],
            ['Sb : px; py; pz']
        ],
        spinors=True,
        unit_cell_cart=[
            [0., 3.2395, 3.2395],
            [3.2395, 0., 3.2395],
            [3.2395, 3.2395, 0.]
        ],
        bands_plot=True,
        atoms_cart=[
            ['In', 0., 0., 0.],
            ['Sb', 1.61975, 1.61975, 1.61975]
        ],
        mp_grid='6 6 6',
        kpoints=[
            list(reversed(x)) for x in 
            itertools.product(
                (0., 1 / 6, 2 / 6, 0.5, -2 / 6, -1 / 6),
                repeat=3
            )
        ],
        kpoint_path=[
            ['G', 0, 0, 0, 'L', 0.5, 0.5, 0.5],
        ]
    )))
    calc.store_all()
    calc.submit()
    print('Submitted calculation', calc.pk)
    

if __name__ == '__main__':
    run_wswannier()
