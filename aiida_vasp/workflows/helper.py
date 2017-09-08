"""Common code for all workflows"""


class WorkflowHelper(object):
    """Helper Class for AiiDA-VASP workflows"""

    def __init__(self, **kwargs):
        self.parent = kwargs['parent']

    def get_parameters(self):
        return self.parent.get_parameters()

    def _get_calc_maker(self, calc_type, **kwargs):
        """
        Instatiate a :py:class:`VaspMaker <aiida_vasp.calcs.maker.VaspMaker>`
        and feed it the common input parameters for all vasp workflows.
        """
        from aiida.orm import Code
        from aiida_vasp.calcs.maker import VaspMaker
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
                maker.kpoints.labels = params.get('kpoint_labels', '')
            elif kpoints.get('path'):
                maker.set_kpoints_path(kpoints['path'])
        elif kp_valid:
            pass
        else:
            raise ValueError(log)

        maker.code = Code.get_from_string(params['vasp_code'])
        maker.label = params.get('label', self._default_label(maker.structure))
        maker.queue = params['queue']
        maker.resources = params['resources']
        return maker

    def _default_label(self, structure):
        label = '[{}]: Calculation for {}'
        label = label.format(self.parent.__class__.__name__,
                             structure.get_formula())
        return label

    @staticmethod
    def _calc_start_msg(name, calc):
        """returns a message that the calculation was started"""
        msg = 'Calculation started: {name}, PK={calc.pk}, uuid={calc.uuid}'
        return msg.format(name=name, calc=calc)

    @staticmethod
    def _subwf_start_msg(name, workflow):
        """returns a message that the subworkflow was started"""
        msg = 'Workflow started: {name}, PK={wf.pk}, uuid={wf.uuid}'
        return msg.format(name=name, wf=workflow)

    def _wf_start_msg(self):
        """returns a message that the workflow was started"""
        params = self.get_parameters()
        msg = 'Workflow started: {wf_class} **{wf_label}**'
        wfclass = self.parent.__class__.__name__
        wflabel = self.parent.label or params.get('label', 'unlabeled')
        return msg.format(wf_class=wfclass, wf_label=wflabel)

    @staticmethod
    def _calc_invalid_outs_msg(calc, links):
        msg = ('Calculation {} does not have all required output nodes. '
               'The required outputs are: {}')
        return msg.format(calc.pk, links)

    def _get_first_step_calc(self, step):
        """return either the first calculation for a step or
        None, if there isn't any"""
        step_calcs = self.parent.get_step_calculations(step)
        if not bool(step_calcs):
            self.parent.append_to_report(
                'no calculation found for step {}'.format(step.__name__))
            return None
        first_calc = step_calcs[0]
        return first_calc

    @staticmethod
    def _verify_calc_output(calc, output_keys):
        """returns True only if all output link names given in
        output_keys are present on calc"""
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
        """finds _verify_param_* methods on the parent instance
        and loops through them to check parameter consistency"""
        valid = True
        log = []

        par_dict = self.parent.__class__.__dict__
        verify_funcs = {
            k: v
            for k, v in par_dict.iteritems() if '_verify_param_' in k
        }

        for _, func in verify_funcs.iteritems():
            if isinstance(func, classmethod):
                func = func.__func__
            valid_i, log_i = func(self.parent, params)
            valid &= valid_i
            log.append(log_i)

        log = filter(None, log)
        log = '\n'.join(log)
        if not silent:
            if any(log):
                raise ValueError(log)
        return valid, log

    @classmethod
    def _verify_params_cls(cls, wf_cls, params, silent=False):
        """same as _verify_params, but looks for classmethods on
        the given workflow class"""
        valid = True
        log = []

        par_dict = wf_cls.__class__.__dict__
        verify_funcs = {
            k: v
            for k, v in par_dict.iteritems()
            if '_verify_param_' in k and isinstance(v, classmethod)
        }

        for _, func in verify_funcs.iteritems():
            valid_i, log_i = func.__func__(wf_cls, params)
            valid &= valid_i
            log.append(log_i)

        log = filter(None, log)
        log = '\n'.join(log)
        if not silent:
            if any(log):
                raise ValueError(log)
        return valid, log

    def _verify_kpoints(self, params):
        """common kpoints parameter consistency check"""
        log = ''
        kpoints = params.get('kpoints')
        if kpoints:
            if len(kpoints) != 1:
                valid = False
                log += (
                    '{}: parameters: kpoints dict must '
                    'contain exactly one item'
                    ).format(self.parent.__class__.__name__)  # yapf: disable
            else:
                valid = bool(kpoints.get('mesh'))
                valid |= bool(kpoints.get('list'))
                valid |= bool(kpoints.get('path'))
                if not valid:
                    log += (
                        '{}: parameters: kpoints dict must '
                        'contain one of the keys '
                        '"mesh", "list" or "path"'
                        ).format(self.parent.__class__.__name__)  # yapf: disable

        elif params.get('continue_from'):
            valid = True
        else:
            log += (
                '{}: parameters: kpoints dict must be set, if '
                'not continuing from a finished calculation'
                ).format(self.parent.__class__.__name__)  # yapf: disable
            valid = False
        return valid, log

    def _verify_paws(self, params):
        """check the paw input parameter for availability of the indicated
        paw nodes in the database."""
        from aiida.orm import DataFactory
        paw_cls = DataFactory('vasp.paw')
        log = ''
        valid = True
        if not params.get('continue_from'):
            paw_fam = params.get('paw_family')
            paw_map = params.get('paw_map')
            fam_valid = bool(paw_fam) and paw_cls.check_family(paw_fam)
            if not fam_valid:
                log += (
                    '{}: parameters: paw_family must be set '
                    'to the name of a PAW family containing '
                    'PAWs necessary for this calculation, if '
                    'not continuing from a finished calculation'
                    ).format(self.parent.__class__.__name__)  # yapf: disable
                valid = False
            map_given = bool(paw_map)
            if not map_given:
                log += ('{}: parameters: paw_map must be a dict '
                        'containing a 1 to 1 mapping of elements '
                        'to symbols, if not continuing from a finished '
                        'calculation.').format(self.parent.__class__.__name__)
                valid = False
        return valid, log

    @classmethod
    def get_params_template(cls, continuation=False):
        """returns a dictionary with some common keys and explanations"""
        tmpl = {}
        tmpl['vasp_code'] = 'code@computer'
        tmpl['queue'] = 'queue name on the remote computer'
        tmpl['resources'] = {
            'num_machines': 'int',
            'num_mpiprocs_per_machine': 'int'
        }
        tmpl['kpoints'] = {'mesh | list | path': ['...']}
        tmpl['label'] = 'optional label for calculations'
        tmpl['description'] = 'optional description for calculations'
        tmpl['extras'] = {
            'explanation': ('dict with extra attributes you '
                            'want to give to all calcs run by this workflow.')
        }
        tmpl['settings'] = {'explanation': ('incar keys for the calculation')}
        if continuation:
            tmpl['continue_from'] = 'uuid of a finished calculation'
            tmpl['kpoints'] = ['default: same as previous calc)']
            tmpl['settings'] = {
                'explanation': ('additional / '
                                'override incar keys')
            }
        else:
            tmpl['paw_family'] = ('name of a PAW family as imported from the '
                                  'commandline')
            tmpl['paw_map'] = {'Element': 'Symbol'}
            tmpl['structure'] = 'POSCAR or CIF file'
        return tmpl

    @classmethod
    def get_template(cls, wf_class=None, path=None):
        """
        Stores the output of the parent's get_params_template method
        in a JSON file
        """
        import json
        from os.path import abspath, expanduser, exists
        result = None
        wf_class = wf_class or cls
        tpl_dict = wf_class.get_params_template()
        if path:
            path = abspath(expanduser(path))
            write = True
            if exists(path):
                overwrite = raw_input(path + ' exists! Overwrite? [y/N]: ')
                return bool(overwrite.lower() in ['y', 'yes'])
            if write:
                with open(path, 'w') as input_tpl:
                    json.dump(tpl_dict, input_tpl, indent=4, sort_keys=True)
                result = path
        else:
            result = json.dumps(tpl_dict, indent=2, sort_keys=True)
        return result
