"""
Bands workchain.

----------------
Intended to be used to extract the band structure using SeeKpath as a preprossesor
to extract the k-point path.
"""

# pylint: disable=attribute-defined-outside-init, import-outside-toplevel
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import WorkChain, append_, calcfunction
from aiida.plugins import WorkflowFactory
from aiida_vasp.utils.aiida_utils import get_data_class, get_data_node
from aiida_vasp.utils.workchains import prepare_process_inputs, compose_exit_code
from aiida_vasp.assistant.parameters import inherit_and_merge_parameters


class BandsWorkChain(WorkChain):
    """Extract the band structure using k-point paths fetched from SeeKpath."""

    _verbose = False
    _next_workchain_string = 'vasp.vasp'
    _next_workchain = WorkflowFactory(_next_workchain_string)

    @classmethod
    def define(cls, spec):
        super(BandsWorkChain, cls).define(spec)
        spec.expose_inputs(cls._next_workchain, exclude=('parameters', 'settings', 'kpoints'))
        spec.input('parameters', valid_type=get_data_class('dict'), required=False)
        spec.input('settings', valid_type=get_data_class('dict'), required=False)
        spec.input('smearing.gaussian', valid_type=get_data_class('bool'), required=False, default=lambda: get_data_node('bool', True))
        spec.input('smearing.sigma', valid_type=get_data_class('float'), required=False, default=lambda: get_data_node('float', 0.05))
        spec.input('restart_folder',
                   valid_type=get_data_class('remote'),
                   required=True,
                   help="""
            The folder to restart in, which contains the files from the prerun to extract the charge density.
            """)
        spec.input('bands.kpoints_distance',
                   valid_type=get_data_class('float'),
                   required=False,
                   default=lambda: get_data_node('float', 0.05),
                   help="""
            The distance between each k-point along each high-symmetry line.
            """)
        spec.input('bands.decompose_bands',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
            Decompose the band structure on each atom.
            """)
        spec.input('bands.decompose_wave',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
            Decompose the wave function.
            """)
        spec.input('bands.lm',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
            Further decompose the decomposition into l- and m-states.
            """)
        spec.input('bands.phase',
                   valid_type=get_data_class('bool'),
                   required=False,
                   default=lambda: get_data_node('bool', False),
                   help="""
            Further decompose the l- and m-state decomposition into phases.
            """)
        spec.input('bands.wigner_seitz_radius',
                   valid_type=get_data_class('list'),
                   required=False,
                   default=lambda: get_data_node('list', list=[False]),
                   help="""
            The Wigner-Seitz radius for each atom type in AA as a list. If set, the internal projectors are not utilzed.
            """)
        spec.outline(
            cls.initialize,
            cls.get_kpoints_path,
            cls.init_next_workchain,
            cls.run_next_workchain,
            cls.verify_next_workchain,
            cls.results,
            cls.finalize
        )  # yapf: disable

        spec.expose_outputs(cls._next_workchain)
        spec.output('bands', valid_type=get_data_class('array.bands'))
        spec.exit_code(0, 'NO_ERROR', message='the sun is shining')
        spec.exit_code(420, 'ERROR_NO_CALLED_WORKCHAIN', message='no called workchain detected')
        spec.exit_code(500, 'ERROR_UNKNOWN', message='unknown error detected in the bands workchain')
        spec.exit_code(2001, 'ERROR_BANDSDATA_NOT_FOUND', message='BandsData not found in exposed_ouputs')

    def initialize(self):
        """Initialize."""
        self._init_context()
        self._init_inputs()
        self._init_settings()

    def _init_context(self):
        """Initialize context variables."""
        self.ctx.exit_code = self.exit_codes.ERROR_UNKNOWN  # pylint: disable=no-member
        self.ctx.inputs = AttributeDict()

    def _init_settings(self):
        """Initialize the settings."""
        # Make sure we parse the bands
        if 'settings' in self.inputs:
            settings = AttributeDict(self.inputs.settings.get_dict())
        else:
            settings = AttributeDict({'parser_settings': {}})
        dict_entry = {'add_bands': True}
        try:
            settings.parser_settings.update(dict_entry)
        except AttributeError:
            settings.parser_settings = dict_entry
        self.ctx.inputs.settings = settings

    def _init_inputs(self):
        """Initialize inputs."""
        self.ctx.inputs.parameters = self._init_parameters()

        # Do not put the SeeKPath parameters in the inputs to avoid port checking
        # of the next workchain
        self.ctx.seekpath_parameters = get_data_node('dict', dict={'reference_distance': self.inputs.bands.kpoints_distance.value})

        try:
            self._verbose = self.inputs.verbose.value
            self.ctx.inputs.verbose = self.inputs.verbose
        except AttributeError:
            pass

    def _init_parameters(self):
        """Collect input to the workchain in the relax namespace and put that into the parameters."""

        # At some point we will replace this with possibly input checking using the PortNamespace on
        # a dict parameter type. As such we remove the workchain input parameters as node entities. Much of
        # the following is just a workaround until that is in place in AiiDA core.
        parameters = inherit_and_merge_parameters(self.inputs)

        # Now we need to make sure we keep the charge density fixed. When executing this
        # workchain we already receive a restart folder, where the charge density resides from the
        # previous run.
        parameters.charge = AttributeDict()
        parameters.charge.constant_charge = True

        return parameters

    def init_next_workchain(self):
        """Initialize the next workchain."""
        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError('No input dictionary was defined in self.ctx.inputs')

        # Add exposed inputs
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        # Make sure we do not have any floating dict (convert to Dict)
        self.ctx.inputs = prepare_process_inputs(self.ctx.inputs, namespaces=['dynamics'])

    def run_next_workchain(self):
        """Run the next workchain."""
        inputs = self.ctx.inputs
        running = self.submit(self._next_workchain, **inputs)

        self.report('launching {}<{}> '.format(self._next_workchain.__name__, running.pk))

        return self.to_context(workchains=append_(running))

    def get_kpoints_path(self):
        """
        Fetch the k-point path.

        Run SeeKpath to get the high symmetry lines of the given structure. This
        routine returns a new (potentially different to the input structure) primitive
        structure. It also returns the k-point path for this structure.
        """
        result = seekpath_structure_analysis(self.inputs.structure, self.ctx.seekpath_parameters)
        self.ctx.inputs.kpoints = result['explicit_kpoints']

    def verify_next_workchain(self):
        """Verify and inherit exit status from child workchains."""

        try:
            workchain = self.ctx.workchains[-1]
        except IndexError:
            self.report('There is no {} in the called workchain list.'.format(self._next_workchain.__name__))
            return self.exit_codes.ERROR_NO_CALLED_WORKCHAIN  # pylint: disable=no-member

        # Inherit exit status from last workchain (supposed to be
        # successfull)
        next_workchain_exit_status = workchain.exit_status
        next_workchain_exit_message = workchain.exit_message
        if not next_workchain_exit_status:
            self.ctx.exit_code = self.exit_codes.NO_ERROR  # pylint: disable=no-member
        else:
            self.ctx.exit_code = compose_exit_code(next_workchain_exit_status, next_workchain_exit_message)
            self.report('The called {}<{}> returned a non-zero exit status. '
                        'The exit status {} is inherited'.format(workchain.__class__.__name__, workchain.pk, self.ctx.exit_code))

        return self.ctx.exit_code

    def results(self):
        """Attach the remaining output results."""

        workchain = self.ctx.workchains[-1]
        out_dict = self.exposed_outputs(workchain, self._next_workchain)
        bands = out_dict.pop('bands', None)
        self.out_many(out_dict)
        if bands is None:
            self.ctx.exit_code = self.exit_codes.ERROR_BANDSDATA_NOT_FOUND  # pylint: disable=no-member
        else:
            self.out('bands', attach_labels(bands, self.ctx.inputs.kpoints))
            self.ctx.exit_code = self.exit_codes.NO_ERROR  # pylint: disable=no-member

        return self.ctx.exit_code

    def finalize(self):
        """Finalize the workchain."""
        return self.ctx.exit_code


@calcfunction
def seekpath_structure_analysis(structure, parameters):
    """
    Workfunction to extract k-points in the reciprocal cell.

    This workfunction will take a structure and pass it through SeeKpath to get the
    primitive cell and the path of high symmetry k-points through its Brillouin zone.
    Note that the returned primitive cell may differ from the original structure in
    which case the k-points are only congruent with the primitive cell.
    """
    from aiida.tools import get_explicit_kpoints_path
    return get_explicit_kpoints_path(structure, **parameters.get_dict())


@calcfunction
def attach_labels(bands, kpoints):
    bands_with_labels = bands.clone()
    bands_with_labels.labels = kpoints.labels
    return bands_with_labels
