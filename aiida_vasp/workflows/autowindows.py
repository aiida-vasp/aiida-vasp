"""
AiiDA - Workflow for investigating optimal Wannier90 window parameters
"""
from aiida.orm import Workflow, WorkflowFactory

from aiida_vasp.utils import compare_bands as bcp
from .helper import WorkflowHelper


class AutowindowsWorkflow(Workflow):
    """Try different inner and outer windows with wannier,
    compare them and choose the best one according to simplistic criteria"""
    Helper = WorkflowHelper
    ScfWf = WorkflowFactory('vasp.scf')
    NscfWf = WorkflowFactory('vasp.nscf')
    ProjWf = WorkflowFactory('vasp.projections')
    WannierWf = WorkflowFactory('vasp.wannier')

    def __init__(self, **kwargs):
        self.helper = self.Helper(parent=self)
        super(AutowindowsWorkflow, self).__init__(**kwargs)

    @Workflow.step
    # pylint: disable=protected-access
    def start(self):
        """Launch the initial SCF VASP sub workflow"""
        self.append_to_report(self.helper._wf_start_msg())
        params = self.get_parameters()
        kpoints = params['kpoints']

        scfpar = self.get_vasp_params(params)
        scfpar['parameters'] = params['parameters']
        scfpar['structure'] = params['structure']
        scfpar['kpoints'] = {'mesh': kpoints['mesh']}
        scfpar['paw_family'] = params['paw_family']
        scfpar['paw_map'] = params['paw_map']

        workflow = self.ScfWf(params=scfpar)
        workflow.label = params.get('label')
        workflow.start()
        self.attach_workflow(workflow)
        self.append_to_report(self.helper._subwf_start_msg('Scf', workflow))

        self.next(self.get_win)

    @Workflow.step
    def get_win(self):
        """Launch the Vasp2Wannier workflow to get the initial Wannier90 input file"""
        start_wf = self.get_step(self.start).get_sub_workflows()[0]
        scf_calc = start_wf.get_result('calc')

        params = self.get_parameters()

        winpar = self.get_vasp_params(params)
        winpar['continue_from'] = scf_calc.uuid
        winpar['use_wannier'] = True

        workflow = self.NscfWf(params=winpar)
        workflow.label = params.get('label')
        workflow.start()
        self.attach_workflow(workflow)
        self.append_to_report(self.helper._subwf_start_msg('Win', workflow))  # pylint: disable=protected-access
        self.next(self.get_projections)

    @Workflow.step
    # pylint: disable=protected-access
    def get_projections(self):
        """Start the Wannier90 subworkflow to obtain projections"""
        win_wf = self.get_step(self.get_win).get_sub_workflows()[0]
        win_calc = win_wf.get_result('calc')

        params = self.get_parameters()
        kppath = self._kppath_vasp_to_wannier(params['kpoints']['path'])

        projpar = self.get_vasp_params(params)
        projpar['continue_from'] = win_calc.uuid
        projpar['projections'] = params['projections']
        projpar['wannier_parameters'] = {
            'num_wann': params['wannier_parameters']['num_wann'],
            'use_bloch_phases': False,
            'bands_plot': True,
            'hr_plot': True,
            'kpoint_path': kppath
        }

        workflow = self.ProjWf(params=projpar)
        workflow.label = params.get('label')
        workflow.start()
        self.attach_workflow(workflow)
        self.append_to_report(self.helper._subwf_start_msg('Proj', workflow))
        self.next(self.get_bands_preview)

    @classmethod
    def _kppath_vasp_to_wannier(cls, kppath):
        """Convert a kpoint path obtained from VASP to Wannier90 format"""
        import itertools
        wannier_kpp = []
        for segment in kppath:
            # flatten the segment list
            wannier_kpp.append(
                list(itertools.chain.from_iterable(segment[:4])))
        return wannier_kpp

    @Workflow.step
    def get_bands_preview(self):
        """Run a VASP DFT calculation to get a rough band structure to automatically set reasonable window parameters"""
        params = self.get_parameters()

        start_wf = self.get_step(self.start).get_sub_workflows()[0]
        scf_calc = start_wf.get_result('calc')

        bandpar = self.get_vasp_params(params)
        bandpar['continue_from'] = scf_calc.uuid
        bandpar['kpoints'] = {'path': params['kpoints']['path']}
        bandpar['use_wannier'] = False

        workflow = self.NscfWf(params=bandpar)
        workflow.label = params.get('label')
        workflow.start()
        self.attach_workflow(workflow)
        self.append_to_report(
            self.helper._subwf_start_msg('Preview-Bands', workflow))  # pylint: disable=protected-access

        self.next(self.make_windows)

    @Workflow.step
    # pylint: disable=too-many-locals,protected-access
    def make_windows(self):
        """Create the window parameter sets and store them into an attribute"""
        params = self.get_parameters()
        num_wann = params['wannier_parameters']['num_wann']

        bandpr_wf = self.get_step(
            self.get_bands_preview).get_sub_workflows()[0]
        bandpr_calc = bandpr_wf.get_result('calc')
        bandpr_bands = bandpr_calc.out.bands
        bands, occupations = bandpr_bands.get_bands(also_occupations=True)
        bands = bcp._firstspin(bands)
        occupations = bcp._firstspin(occupations)
        e_fermi = bandpr_calc.out.results.get_attr('efermi')

        bandgap = bcp.band_gap(bands, occupations, efermi=e_fermi)
        bmin = [bands[:, i].min() for i in range(bands.shape[1])]
        bmax = [bands[:, i].max() for i in range(bands.shape[1])]

        # find minimum inner window by looking at band gap energies
        import math
        iw_min = math.floor(bandgap['vector'][0][1])
        iw_max = math.ceil(bandgap['vector'][1][1])

        # find outer widow by counting nwann / 2 bands down and nwann / 2 up
        # from efermi
        nwann_half = int(math.ceil(num_wann / 2))

        lower_max = [i for i in bmax if i < e_fermi]
        lower_max.sort(reverse=True)
        ow_min = math.floor(bmin[bmax.index(lower_max[nwann_half - 1])])

        upper_min = sorted([i for i in bmin if i > e_fermi])
        ow_max = math.ceil(bmax[bmin.index(upper_min[nwann_half - 1])])

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
        """Get the tight binding models and band structures for them for all the requested window parameter sets"""
        proj_wf = self.get_step(self.get_projections).get_sub_workflows()[0]
        proj_calc = proj_wf.get_result('calc')

        params = self.get_parameters()

        count = 0
        wpar = self.get_wannier_params(params)
        wpar['continue_from'] = proj_calc.uuid
        wpar['parameters']['bands_plot'] = True
        wpar['parameters']['hr_plot'] = True

        wfpk = []
        for window in self.get_attribute('windows'):
            wpar['parameters']['dis_win_min'] = window['outer'][0]
            wpar['parameters']['dis_win_max'] = window['outer'][1]
            wpar['parameters']['dis_froz_min'] = window['inner'][0]
            wpar['parameters']['dis_froz_max'] = window['inner'][1]

            workflow = self.WannierWf(params=wpar)
            workflow.label = params.get('label')
            workflow.start()
            self.attach_workflow(workflow)
            count += 1
            wfpk.append(workflow.pk)

        self.append_to_report('running tbmodels for {} windows'.format(count))
        self.append_to_report(
            'tbmodels pk-range: {} - {}'.format(wfpk[0], wfpk[-1]))

        self.next(self.get_reference_bands)

    @Workflow.step
    def get_reference_bands(self):
        """Calculate bandstructure with VASP by DFT as reference for the tight-binding Wannier90 bands structures"""
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

        workflow = self.NscfWf(params=bandpar)
        workflow.label = params.get('label')
        workflow.start()
        self.attach_workflow(workflow)
        self.append_to_report(
            self.helper._subwf_start_msg('Ref-Bands', workflow))  # pylint: disable=protected-access

        self.next(self.gather_results)

    @Workflow.step
    def gather_results(self):
        """Set the results of the workflow and the 'wbands_list' attribute"""
        self.append_to_report('retrieving and compiling results')
        wannier_wf_list = self.get_step(self.get_tbmodel).get_sub_workflows()
        band_wf = self.get_step(
            self.get_reference_bands).get_sub_workflows()[0]
        self.add_result('reference_bands',
                        band_wf.get_result('calc').out.bands)
        self.add_result('reference_calc', band_wf.get_result('calc'))

        wbands_list = []
        for workflow in wannier_wf_list:
            try:
                calc = workflow.get_result('calc')
                bands = workflow.get_result('bands')
                wbands_list.append(bands.uuid)
                self.add_result('bands_{}'.format(calc.pk), bands)
            except Exception as err:  # pylint: disable=broad-except
                wset = workflow.get_parameters()['parameters']
                window = 'inner: {}-{}, outer: {}-{}'.format(
                    wset['dis_froz_min'], wset['dis_froz_max'],
                    wset['dis_win_min'], wset['dis_win_max'])
                self.append_to_report(('workflow {pk} with window {window} '
                                       'did not yield the expected results: \n'
                                       '{error}').format(
                                           pk=workflow.pk,
                                           window=window,
                                           error=repr(err)))

        self.add_attribute('wbands_list', wbands_list)

        self.next(exit)

    @classmethod
    def get_general_params(cls, params):
        """Get parameters that pertain to all runs"""
        genpar = {}
        genpar['extras'] = params.get('extras', {}).copy()
        # ~ genpar['extras']['wf_uuid'] = unicode(cls.uuid)
        genpar['label'] = params.get('label')
        genpar['description'] = params.get('description')
        return genpar

    @classmethod
    def get_vasp_params(cls, params):
        """Extract the parameters for the vasp runs"""
        vasppar = cls.get_general_params(params)
        vasppar['vasp_code'] = params['vasp_code']
        vasppar['resources'] = params['resources']
        vasppar['queue'] = params['queue']
        return vasppar

    @classmethod
    def get_wannier_params(cls, params):
        """Extract the parameters for the wannier runs"""
        resources = params.get('wannier_resources', params['resources'].copy())
        kppath = cls._kppath_vasp_to_wannier(params['kpoints']['path'])
        queue = params.get('wannier_queue', params['queue'])

        wpar = cls.get_general_params(params)
        wpar['wannier_code'] = params['wannier_code']
        wpar['parameters'] = params['wannier_parameters'].copy()
        wpar['parameters']['kpoint_path'] = kppath
        wpar['resources'] = resources
        wpar['queue'] = queue
        return wpar

    @classmethod
    def get_template(cls, *args, **kwargs):
        """returns a JSON formatted string that could be stored
        in a file, edited, loaded and used as parameters to run
        this workflow."""
        return cls.Helper.get_template(*args, wf_class=cls, **kwargs)

    @classmethod
    def get_params_template(cls):
        """returns a dictionary with the necessary keys to
        run this workflow and explanations to each key as values"""
        tmpl = cls.Helper.get_params_template()
        wtpl = cls.WannierWf.get_params_template()
        ptpl = cls.ProjWf.get_params_template()
        tmpl['wannier_parameters'] = {'num_wann': 'int', 'hr_plot': True}
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
        """Make sure the resources parameter matches with the 'nbands' setting"""
        valid = True
        log = ''
        nbands = params['parameters'].get('nbands')
        if nbands:
            res = params['resources']
            nproc = res['num_machines'] * res['num_mpiprocs_per_machine']
            if (nbands % nproc) != 0:
                valid = False
                log += ('nbands is not divisible by '
                        'num_machines * num_mpiprocs_per_machine')
        return valid, log

    def set_params(self, params, force=False):
        self.helper._verify_params(params)  # pylint: disable=protected-access
        super(AutowindowsWorkflow, self).set_params(params, force=force)
