from aiida.orm.workflow import Workflow
import sys


class WorkflowBase(Workflow):
    '''Base Class for AiiDA-VASP workflows'''
    def __init__(self, **kwargs):
        super(WorkflowBase, self).__init__(**kwargs)

    def _get_calc_maker(self, calc_type, **kwargs):
        from aiida.orm import Code
        from aiida.orm.calculation.job.vasp.maker import VaspMaker
        params = self.get_parameters()
        good_init = bool(kwargs.get('continue_from'))
        good_init |= bool(kwargs.get('copy_from'))
        good_init |= bool(kwargs.get('structure'))
        if not good_init:
            kwargs['structure'] = params['structure']
        maker = VaspMaker(calc_cls=calc_type, **kwargs)
        if params.get('settings'):
            maker.rewrite_settings(**params['settings'])

        kpoints = params.get('kpoints')
        kp_valid, log = self._verify_kpoints(params)
        if kp_valid:
            if kpoints.get('mesh'):
                maker.set_kpoints_mesh(kpoints['mesh'])
            elif kpoints.get('list'):
                maker.set_kpoints_list(kpoints['list'])
            elif kpoints.get('path'):
                maker.set_kpoints_path(kpoints['path'])
        else:
            raise ValueError(log)

        maker.code = Code.get_from_string(params['vasp_code'])
        maker.label = params.get('label', self.default_label(maker.structure))
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
        msg.format(calc.pk, links)

    def _get_first_step_calc(self, step):
        step_calcs = self.get_step_calculations(step)
        if not bool(step_calcs):
            self.append_to_report('no calculation found for step {}'.format(
                step.__name__))
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

    def set_params(self, params, **kwargs):
        self._verify_params(params)
        super(WorkflowBase, self).set_params(params, **kwargs)

    def _verify_params(self, params):
        valid = True
        log = []

        kp = self._verify_kpoints(self, params)
        valid &= kp[0]
        log.append(kp[1])

        paw = self._verify_paws(self, params)
        valid &= paw[0]
        log.append(paw[1])

        print >> sys.stderr, '\n'.join(log)
        return valid

    def _verify_kpoints(self, params):
        log = ''
        kpoints = params.get('kpoints')
        if kpoints:
            if len(kpoints) != 1:
                valid = False
                log += ('{}: parameters: kpoints dict must '
                        'contain exactly one item').format(
                            self.__class__.__name__)
            else:
                valid = bool(kpoints.get('mesh'))
                valid |= bool(kpoints.get('list'))
                valid |= bool(kpoints.get('path'))
                if not valid:
                    log += ('{}: parameters: kpoints dict must '
                            'contain one of the keys '
                            '"mesh", "list" or "path"').format(
                                self.__class__.__name__)

        elif params.get('continue_from'):
            valid = True
        else:
            log += ('{}: parameters: kpoints dict must be set, if '
                    'not continuing from a finished calculation').format(
                        self.__class__.__name__)
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
                        self.__class__.__name__)
                valid = False
            map_given = bool(paw_map)
            if not map_given:
                log += ('{}: parameters: paw_map must be a dict '
                        'containing a 1 to 1 mapping of elements '
                        'to symbols, if not continuing from a finished '
                        'calculation.').format(self.__class__.__name__)
                valid = False
        return valid, log

    def get_params_template(self, continuation=False):
        # ~ item = lambda desc, typ, req: '[{req}] ({type}) {desc}'.format(
            # ~ desc=desc, type=typ, req=req)
        # ~ required = 'required'
        # ~ optional = 'optional'
        # ~ req_if_cont = continuation and required, optional
        # ~ req_if_orig = continuation and optional, required
        tmpl = {}
        tmpl['vasp_code'] = 'code in the db that runs vasp'
        tmpl['queue'] = 'queue name on the remote computer'
        tmpl['resources'] = 'aiida style resource dict'
        tmpl['kpoints'] = '{("mesh" | "list" | "path"): [...]}'
        tmpl['label'] = 'optional: override default label'
        tmpl['description'] = 'optional: override default description'
        tmpl['extras'] = ('dict with extra attributes you want to give '
                          'to all calcs run by this workflow.')
        if continuation:
            tmpl['continue_from'] = 'uuid of a finished calculation'
            tmpl['paw_family'] = 'name of a PAW family as imported from the commandline'
            tmpl['paw_map'] = 'dict with element keys mapping to symbol values'
            tmpl['kpoints'] += ' (default: same as previous calc)'
        return tmpl
