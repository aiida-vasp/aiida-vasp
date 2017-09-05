from aiida.orm import Workflow, WorkflowFactory
from helper import WorkflowHelper
from aiida_vasp.utils import compare_bands as bcp


class AutowindowsWorkflow(Workflow):
    '''Try different inner and outer windows with wannier,
    compare them and choose the best one according to simplistic criteria'''
    Helper = WorkflowHelper
    ScfWf = WorkflowFactory('vasp.scf')
    NscfWf = WorkflowFactory('vasp.nscf')
    ProjWf = WorkflowFactory('vasp.projections')
    WannierWf = WorkflowFactory('vasp.wannier')

    def __init__(self, **kwargs):
        self.helper = self.Helper(parent=self)
        super(AutowindowsWorkflow, self).__init__(**kwargs)

    @Workflow.step
    def start(self):
        self.append_to_report(self.helper._wf_start_msg())
        params = self.get_parameters()
        kp = params['kpoints']

        scfpar = self.get_vasp_params(params)
        scfpar['settings'] = params['settings']
        scfpar['structure'] = params['structure']
        scfpar['kpoints'] = {'mesh': kp['mesh']}
        scfpar['paw_family'] = params['paw_family']
        scfpar['paw_map'] = params['paw_map']

        wf = self.ScfWf(params=scfpar)
        wf.label = params.get('label')
        wf.start()
        self.attach_workflow(wf)
        self.append_to_report(self.helper._subwf_start_msg('Scf', wf))

        self.next(self.get_win)

    @Workflow.step
    def get_win(self):
        start_wf = self.get_step(self.start).get_sub_workflows()[0]
        scf_calc = start_wf.get_result('calc')

        params = self.get_parameters()

        winpar = self.get_vasp_params(params)
        winpar['continue_from'] = scf_calc.uuid
        winpar['use_wannier'] = True

        wf = self.NscfWf(params=winpar)
        wf.label = params.get('label')
        wf.start()
        self.attach_workflow(wf)
        self.append_to_report(self.helper._subwf_start_msg('Win', wf))
        self.next(self.get_projections)

    @Workflow.step
    def get_projections(self):
        win_wf = self.get_step(self.get_win).get_sub_workflows()[0]
        win_calc = win_wf.get_result('calc')

        params = self.get_parameters()
        kppath = self._kppath_vasp_to_wannier(params['kpoints']['path'])

        projpar = self.get_vasp_params(params)
        projpar['continue_from'] = win_calc.uuid
        projpar['projections'] = params['projections']
        projpar['wannier_settings'] = {
            'num_wann': params['wannier_settings']['num_wann'],
            'use_bloch_phases': False,
            'bands_plot': True,
            'hr_plot': True,
            'kpoint_path': kppath
        }

        wf = self.ProjWf(params=projpar)
        wf.label = params.get('label')
        wf.start()
        self.attach_workflow(wf)
        self.append_to_report(self.helper._subwf_start_msg('Proj', wf))
        self.next(self.get_bands_preview)

    @classmethod
    def _kppath_vasp_to_wannier(cls, kppath):
        import itertools
        wannier_kpp = []
        for segment in kppath:
            # flatten the segment list
            wannier_kpp.append(
                list(itertools.chain.from_iterable(segment[:4])))
        return wannier_kpp

    @Workflow.step
    def get_bands_preview(self):
        params = self.get_parameters()

        start_wf = self.get_step(self.start).get_sub_workflows()[0]
        scf_calc = start_wf.get_result('calc')

        bandpar = self.get_vasp_params(params)
        bandpar['continue_from'] = scf_calc.uuid
        bandpar['kpoints'] = {'path': params['kpoints']['path']}
        bandpar['use_wannier'] = False

        wf = self.NscfWf(params=bandpar)
        wf.label = params.get('label')
        wf.start()
        self.attach_workflow(wf)
        self.append_to_report(
            self.helper._subwf_start_msg('Preview-Bands', wf))

        self.next(self.make_windows)

    @Workflow.step
    def make_windows(self):
        params = self.get_parameters()
        num_wann = params['wannier_settings']['num_wann']

        bandpr_wf = self.get_step(
            self.get_bands_preview).get_sub_workflows()[0]
        bandpr_calc = bandpr_wf.get_result('calc')
        bandpr_bands = bandpr_calc.out.bands
        b, o = bandpr_bands.get_bands(also_occupations=True)
        b = bcp._firstspin(b)
        o = bcp._firstspin(o)
        ef = bandpr_calc.out.results.get_attr('efermi')

        bandgap = bcp.band_gap(b, o, efermi=ef)
        bmin = [b[:, i].min() for i in range(b.shape[1])]
        bmax = [b[:, i].max() for i in range(b.shape[1])]

        # find minimum innter window by looking at band gap energies
        import math
        iw_min = math.floor(bandgap['vector'][0][1])
        iw_max = math.ceil(bandgap['vector'][1][1])

        # find outer widow by counting nwann / 2 bands down and nwann / 2 up
        # from efermi
        n = int(math.ceil(num_wann / 2))

        lower_max = [i for i in bmax if i < ef]
        lower_max.sort(reverse=True)
        ow_min = math.floor(bmin[bmax.index(lower_max[n - 1])])

        upper_min = sorted([i for i in bmin if i > ef])
        ow_max = math.ceil(bmax[bmin.index(upper_min[n - 1])])

        windows = []
        for i in range(params['num_owindows']):
            ow_offset = i * params['owindows-increment']
            owindow = [ow_min - ow_offset, ow_max + ow_offset]
            for j in range(params['num_iwindows']):
                iw_offset = j * params['iwindows-increment']
                iwindow = [iw_min - iw_offset, iw_max + iw_offset]
                windows.append({'outer': owindow, 'inner': iwindow})

        self.add_attribute('windows', windows)

        self.next(self.get_tbmodel)

    @Workflow.step
    def get_tbmodel(self):
        proj_wf = self.get_step(self.get_projections).get_sub_workflows()[0]
        proj_calc = proj_wf.get_result('calc')

        params = self.get_parameters()

        count = 0
        wpar = self.get_wannier_params(params)
        wpar['continue_from'] = proj_calc.uuid
        wpar['settings']['bands_plot'] = True
        wpar['settings']['hr_plot'] = True

        wfpk = []
        for window in self.get_attribute('windows'):
            wpar['settings']['dis_win_min'] = window['outer'][0]
            wpar['settings']['dis_win_max'] = window['outer'][1]
            wpar['settings']['dis_froz_min'] = window['inner'][0]
            wpar['settings']['dis_froz_max'] = window['inner'][1]

            wf = self.WannierWf(params=wpar)
            wf.label = params.get('label')
            wf.start()
            self.attach_workflow(wf)
            count += 1
            wfpk.append(wf.pk)

        self.append_to_report('running tbmodels for {} windows'.format(count))
        self.append_to_report(
            'tbmodels pk-range: {} - {}'.format(wfpk[0], wfpk[-1]))

        self.next(self.get_reference_bands)

    @Workflow.step
    def get_reference_bands(self):
        wannier_wf = self.get_step(self.get_tbmodel).get_sub_workflows()[0]
        wannier_bands = wannier_wf.get_result('bands')
        start_wf = self.get_step(self.start).get_sub_workflows()[0]
        scf_calc = start_wf.get_result('calc')
        kplist = wannier_bands.get_kpoints().tolist()
        kplabels = wannier_bands.labels

        params = self.get_parameters()

        bandpar = self.get_vasp_params(params)
        bandpar['continue_from'] = scf_calc.uuid
        bandpar['kpoints'] = {'list': kplist}
        bandpar['kpoint_labels'] = kplabels
        bandpar['use_wannier'] = False

        wf = self.NscfWf(params=bandpar)
        wf.label = params.get('label')
        wf.start()
        self.attach_workflow(wf)
        self.append_to_report(self.helper._subwf_start_msg('Ref-Bands', wf))

        self.next(self.gather_results)

    @Workflow.step
    def gather_results(self):
        self.append_to_report('retrieving and compiling results')
        wannier_wf_list = self.get_step(self.get_tbmodel).get_sub_workflows()
        band_wf = self.get_step(
            self.get_reference_bands).get_sub_workflows()[0]
        self.add_result('reference_bands',
                        band_wf.get_result('calc').out.bands)
        self.add_result('reference_calc', band_wf.get_result('calc'))

        wbands_list = []
        for wf in wannier_wf_list:
            try:
                calc = wf.get_result('calc')
                bands = wf.get_result('bands')
                wbands_list.append(bands.uuid)
                self.add_result('bands_{}'.format(calc.pk), bands)
            except Exception as e:
                wset = wf.get_parameters()['settings']
                window = 'inner: {}-{}, outer: {}-{}'.format(
                    wset['dis_froz_min'], wset['dis_froz_max'],
                    wset['dis_win_min'], wset['dis_win_max'])
                self.append_to_report(('workflow {pk} with window {window} '
                                       'did not yield the expected results: \n'
                                       '{error}').format(
                                           pk=wf.pk,
                                           window=window,
                                           error=repr(e)))

        self.add_attribute('wbands_list', wbands_list)

        self.next(exit)

    @classmethod
    def get_general_params(cls, params):
        genpar = {}
        genpar['extras'] = params.get('extras', {}).copy()
        # ~ genpar['extras']['wf_uuid'] = unicode(cls.uuid)
        genpar['label'] = params.get('label')
        genpar['description'] = params.get('description')
        return genpar

    @classmethod
    def get_vasp_params(cls, params):
        vasppar = cls.get_general_params(params)
        vasppar['vasp_code'] = params['vasp_code']
        vasppar['resources'] = params['resources']
        vasppar['queue'] = params['queue']
        return vasppar

    @classmethod
    def get_wannier_params(cls, params):
        resources = params.get('wannier_resources', params['resources'].copy())
        kppath = cls._kppath_vasp_to_wannier(params['kpoints']['path'])
        queue = params.get('wannier_queue', params['queue'])

        wpar = cls.get_general_params(params)
        wpar['wannier_code'] = params['wannier_code']
        wpar['settings'] = params['wannier_settings'].copy()
        wpar['settings']['kpoint_path'] = kppath
        wpar['resources'] = resources
        wpar['queue'] = queue
        return wpar

    @classmethod
    def get_template(cls, *args, **kwargs):
        '''returns a JSON formatted string that could be stored
        in a file, edited, loaded and used as parameters to run
        this workflow.'''
        return cls.Helper.get_template(*args, wf_class=cls, **kwargs)

    @classmethod
    def get_params_template(cls):
        '''returns a dictionary with the necessary keys to
        run this workflow and explanations to each key as values'''
        tmpl = cls.Helper.get_params_template()
        wtpl = cls.WannierWf.get_params_template()
        ptpl = cls.ProjWf.get_params_template()
        tmpl['wannier_settings'] = {'num_wann': 'int', 'hr_plot': True}
        tmpl['wannier_code'] = wtpl['wannier_code']
        tmpl['projections'] = ptpl['projections']
        tmpl['iwindows-increment'] = ('eV increment between '
                                      'different inner windows')
        tmpl['num_iwindows'] = 'number of inner windows to try'
        tmpl['owindows-increment'] = ('eV increment between '
                                      'different outer windows')
        tmpl['num_owindows'] = 'number of outer windows to try'
        tmpl['kpoints'] = {'mesh': [], 'path': []}
        tmpl['#kpoints'] = (
            'mesh for everything up until wannier_setup'
            'path for bands, format: [["A", [...], "B", [...]], [...]]')
        return tmpl

    @classmethod
    def _verify_param_resources(cls, params):
        valid = True
        log = ''
        nbands = params['settings'].get('nbands')
        if nbands:
            res = params['resources']
            nproc = res['num_machines'] * res['num_mpiprocs_per_machine']
            if (nbands % nproc) != 0:
                valid = False
                log += ('nbands is not divisible by '
                        'num_machines * num_mpiprocs_per_machine')
        return valid, log

    def set_params(self, params):
        self.helper._verify_params(params)
        super(AutowindowsWorkflow, self).set_params(params)
