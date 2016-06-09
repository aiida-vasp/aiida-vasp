# ~ from aiida.orm.workflow import Workflow
import sys


# ~ class WorkflowHelper(Workflow):
class WorkflowHelper(object):
    '''Base Class for AiiDA-VASP workflows'''
    def __init__(self, **kwargs):
        # ~ super(WorkflowHelper, self).__init__(**kwargs)
        self.parent = kwargs['parent']

    def get_parameters(self):
        return self.parent.get_parameters()

    def _get_calc_maker(self, calc_type, **kwargs):
        from aiida.orm import Code
        from aiida.orm.calculation.job.vasp.maker import VaspMaker
        params = self.get_parameters()
        good_init = bool(kwargs.get('continue_from'))
        good_init |= bool(kwargs.get('copy_from'))
        good_init |= bool(kwargs.get('structure'))
        if not good_init:
            kwargs['structure'] = params['structure']
            kwargs['paw_family'] = kwargs.get('paw_family',
                                              params['paw_family'])
            kwargs['paw_map'] = kwargs.get('paw_map') or params['paw_map']
        maker = VaspMaker(calc_cls=calc_type, **kwargs)
        if params.get('settings'):
            maker.rewrite_settings(**params['settings'])

        kpoints = params.get('kpoints')
        kp_valid, log = self._verify_kpoints(params)
        if kpoints and kp_valid:
            if kpoints.get('mesh'):
                maker.set_kpoints_mesh(kpoints['mesh'])
            elif kpoints.get('list'):
                maker.set_kpoints_list(kpoints['list'])
            elif kpoints.get('path'):
                maker.set_kpoints_path(kpoints['path'])
        elif kwargs.get('continue_from') and kp_valid:
            pass
        else:
            raise ValueError(log)

        maker.code = Code.get_from_string(params['vasp_code'])
        maker.label = params.get('label', self._default_label(maker.structure))
        maker.queue = params['queue']
        maker.resources = params['resources']
        return maker

    def _default_label(self, structure):
        label = '[{}]: VASP scf run for {}'
        label = label.format(self.__class__.__name__,
                             structure.get_formula())
        return label

    def _calc_start_msg(self, name, calc):
        msg = 'Calculation started: {name}, PK={calc.pk}, uuid={calc.uuid}'
        return msg.format(name=name, calc=calc)

    def _calc_invalid_outs_msg(self, calc, links):
        msg = ('Calculation {} does not have all required output nodes. '
               'The required outputs are: {}')
        return msg.format(calc.pk, links)

    def _get_first_step_calc(self, step):
        step_calcs = self.parent.get_step_calculations(step)
        if not bool(step_calcs):
            self.parent.append_to_report(
                'no calculation found for step {}'.format(step.__name__))
            return None
        first_calc = step_calcs[0]
        return first_calc

    def _verify_calc_output(self, calc, output_keys):
        if not bool(calc):
            return False
        out = calc.get_outputs_dict()
        valid_out = True
        for k in output_keys:
            valid_out &= bool(out.get(k))
        return valid_out

    # ~ def set_params(self, params, **kwargs):
        # ~ self._verify_params(params)
        # ~ self.parent.set_params(params, **kwargs)

    def _verify_params(self, params, silent=False):
        valid = True
        log = []

        par_dict = self.parent.__class__.__dict__
        verify_funcs = {k: v for k, v in par_dict.iteritems()
                        if '_verify_param_' in k}

        for name, func in verify_funcs.iteritems():
            valid_i, log_i = func(self.parent, params)
            valid &= valid_i
            log.append(log_i)

        log = '\n'.join(log)
        print >> sys.stderr, log
        return valid, log

    def _verify_kpoints(self, params):
        log = ''
        kpoints = params.get('kpoints')
        if kpoints:
            if len(kpoints) != 1:
                valid = False
                log += ('{}: parameters: kpoints dict must '
                        'contain exactly one item').format(
                            self.parent.__class__.__name__)
            else:
                valid = bool(kpoints.get('mesh'))
                valid |= bool(kpoints.get('list'))
                valid |= bool(kpoints.get('path'))
                if not valid:
                    log += ('{}: parameters: kpoints dict must '
                            'contain one of the keys '
                            '"mesh", "list" or "path"').format(
                                self.parent.__class__.__name__)

        elif params.get('continue_from'):
            valid = True
        else:
            log += ('{}: parameters: kpoints dict must be set, if '
                    'not continuing from a finished calculation').format(
                        self.parent.__class__.__name__)
            valid = False
        return valid, log

    def _verify_paws(self, params):
        from aiida.orm import DataFactory
        PawData = DataFactory('vasp.paw')
        log = ''
        valid = True
        if not params.get('continue_from'):
            paw_fam = params.get('paw_family')
            paw_map = params.get('paw_map')
            fam_valid = bool(paw_fam) and PawData.check_family(paw_fam)
            if not fam_valid:
                log += ('{}: parameters: paw_family must be set '
                        'to the name of a PAW family containing '
                        'PAWs necessary for this calculation, if '
                        'not continuing from a finished calculation').format(
                    self.parent.__class__.__name__)
                valid = False
            map_given = bool(paw_map)
            if not map_given:
                log += ('{}: parameters: paw_map must be a dict '
                        'containing a 1 to 1 mapping of elements '
                        'to symbols, if not continuing from a finished '
                        'calculation.').format(self.parent.__class__.__name__)
                valid = False
        return valid, log

    def get_params_template(cls, continuation=False):
        tmpl = {}
        tmpl['vasp_code'] = 'code@computer'
        tmpl['queue'] = 'queue name on the remote computer'
        tmpl['resources'] = {'num_machines': 'int',
                             'num_mpiprocs_per_machine': 'int'}
        tmpl['kpoints'] = {'mesh | list | path': ['...']}
        tmpl['label'] = 'optional label for calculations'
        tmpl['description'] = 'optional description for calculations'
        tmpl['extras'] = {'explanation': ('dict with extra attributes you '
                          'want to give to all calcs run by this workflow.')}
        tmpl['settings'] = {'explanation': ('incar keys for the calculation')}
        if continuation:
            tmpl['continue_from'] = 'uuid of a finished calculation'
            tmpl['kpoints'] = ['default: same as previous calc)']
            tmpl['settings'] = {'explanation': ('additional / override incar keys')}
        else:
            tmpl['paw_family'] = ('name of a PAW family as imported from the '
                                  'commandline')
            tmpl['paw_map'] = {'Element': 'Symbol'}
            tmpl['structure'] = 'POSCAR or CIF file'
        return tmpl

    # ~ @classmethod
    def get_template(self, path=None):
        import json
        from os.path import abspath, expanduser, exists
        result = None
        tpl_dict = self.parent.get_params_template()
        if path:
            if not exists(path):
                path = abspath(expanduser(path))
                with open(path, 'w') as input_tpl:
                    json.dump(tpl_dict, input_tpl, indent=4, sort_keys=True)
                result = path
        else:
            result = json.dumps(tpl_dict, indent=2, sort_keys=True)
        return result
