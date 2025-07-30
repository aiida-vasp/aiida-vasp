"""
VASP NEB workchain.

Contains the VaspNEBWorkChain class definition which uses the BaseRestartWorkChain.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from aiida import __version__ as aiida_version
from aiida import orm
from aiida.common.exceptions import InputValidationError, NotExistent
from aiida.common.extendeddicts import AttributeDict
from aiida.common.lang import override
from aiida.engine import ExitCode, ProcessSpec, while_
from aiida.engine.processes.workchains.restart import BaseRestartWorkChain, ProcessHandlerReport, process_handler
from aiida.plugins import CalculationFactory

# pylint: disable=too-many-branches, too-many-statements
from packaging import version

from aiida_vasp.assistant.parameters import ParametersMassage
from aiida_vasp.calcs.neb import VaspNEBCalculation
from aiida_vasp.data.potcar import PotcarData
from aiida_vasp.parsers.content_parsers.potcar import MultiPotcarIo
from aiida_vasp.utils.workchains import compose_exit_code

# Additional tags for VTST calculations - these are not the tags used by standard VASP
VTST_ADDITIONAL_TAGS = {
    'iopt': 'TAG for VTST',
    'maxmove': 'Maximum ionic movement',
    'ilbfgsmem': 'Number of steps saved when building the inverse Hessian matrix',
    'lglobal': 'Optimize the NEB globally instead of image-by-image',
    'lautoscale': 'Automatically determines INVCURV',
    'invcurv': 'Initial inverse curvature, used to construct the inverse Hessian matrix',
    'llineopt': 'Use a force based line minimizer for translation',
    'fdstep': 'Finite difference step size for line optimizer',
    'timestep': 'Dynamical time step',
    'sdalpha': 'Ratio between force and step size',
    'ftimemax': 'Maximum dynamical time step allowed',
    'ftimedec': 'Factor to decrease dt',
    'ftimeinc': 'Factor to increase dt',
    'falpha': 'Parameter that controls velocity damping',
    'fnmin': 'Minium number of iterations before adjust alpha and dt',
    'lclimb': 'Use climbing image mode',
    'ichain': 'Indicates which method to run. NEB (ICHAIN=0) is the default',
    'ltangentold': 'Flag to turn on the old central difference tangent',
    'ldneb': 'Flag to turn on modified double nudging',
    'lnebcell': 'Flag to turn on SS-NEB. Used with ISIF=3 and IOPT=3.',
    'jacobian': 'Controls weight of lattice to atomic motion. Ω is volume and N is the number of atoms.',
}


class VaspNEBWorkChain(BaseRestartWorkChain):
    """
    A NEB workchain

    Error handling enriched wrapper around VaspNEBCalculation.

    Deliberately conserves most of the interface (required inputs) of the VaspNEBCalculation class, but
    makes it possible for a user to interact with a workchain and not a calculation.

    In addition, implement restarts of calculation when the calculation is net full converged for error handling.

    """

    _verbose = False
    _process_class = CalculationFactory('vasp.neb')
    _norm_disp_threshold = 4.0

    @classmethod
    def define(cls, spec: ProcessSpec) -> None:
        super(VaspNEBWorkChain, cls).define(spec)
        spec.expose_inputs(
            cls._process_class,
            exclude=('potential', 'kpoints', 'dynamics', 'metadata'),
        )
        spec.input(
            'kpoints',
            valid_type=orm.KpointsData,
            required=False,
        )
        spec.input(
            'kpoints_spacing',
            valid_type=orm.Float,
            required=False,
        )
        spec.input(
            'potential_family',
            valid_type=orm.Str,
            required=True,
        )
        spec.input(
            'potential_mapping',
            valid_type=orm.Dict,
            required=True,
        )
        spec.input(
            'options',
            valid_type=orm.Dict,
            required=True,
        )
        spec.input(
            'max_iterations',
            valid_type=orm.Int,
            required=False,
            default=lambda: orm.Int(5),
            help="""
            The maximum number of iterations to perform.
            """,
        )
        spec.input(
            'clean_workdir',
            valid_type=orm.Bool,
            required=False,
            default=lambda: orm.Bool(False),
            help="""
            If True, clean the work dir upon the completion of a successful calculation.
            """,
        )
        spec.input(
            'verbose',
            valid_type=orm.Bool,
            required=False,
            default=lambda: orm.Bool(True),
            help="""
            If True, enable more detailed output during workchain execution.
            """,
        )
        spec.input(
            'dynamics.positions_dof',
            valid_type=orm.List,
            required=False,
            help="""
            Site dependent flag for selective dynamics when performing relaxation
            """,
        )
        spec.input(
            'ldau_mapping',
            valid_type=orm.Dict,
            required=False,
            help="Mappings, see the doc string of 'get_ldau_keys'",
        )
        spec.input(
            'kpoints_spacing',
            valid_type=orm.Float,
            required=False,
            help='Spacing for the kpoints in units A^-1 * 2pi (CASTEP style `kpoints_mp_spacing`)',
        )
        spec.input(
            'kpoints_spacing_vasp',
            valid_type=orm.Float,
            required=False,
            help='Spacing for the kpoints in units A^-1 (VASP style)',
        )
        spec.outline(
            cls.setup,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )

        spec.expose_outputs(cls._process_class)
        spec.exit_code(
            0,
            'NO_ERROR',
            message='the sun is shining',
        )
        spec.exit_code(
            700,
            'ERROR_NO_POTENTIAL_FAMILY_NAME',
            message='the user did not supply a potential family name',
        )
        spec.exit_code(
            701,
            'ERROR_POTENTIAL_VALUE_ERROR',
            message='ValueError was returned from get_potcars_from_structure',
        )
        spec.exit_code(
            702,
            'ERROR_POTENTIAL_DO_NOT_EXIST',
            message='the potential does not exist',
        )
        spec.exit_code(
            703,
            'ERROR_IN_PARAMETER_MASSAGER',
            message='the exception: {exception} was thrown while massaging the parameters',
        )
        spec.exit_code(
            501,
            'SUB_NEB_CALCULATION_ERROR',
            message='Unrecoverable error in launched NEB calculations.',
        )

    def setup(self) -> None | ExitCode:
        super().setup()

        # Setup the initial inputs
        self.ctx.inputs = self.exposed_inputs(self._process_class)

        # Stage the neb images
        self.ctx.neb_images = self.inputs.neb_images

        # Handle and convert additional inputs and store them in self.ctx.inputs
        exit_code = self._setup_vasp_inputs()
        if exit_code is not None:
            return exit_code

        # Sanity checks
        self._check_neb_inputs()
        return None

    @process_handler(priority=500, exit_codes=[VaspNEBCalculation.exit_codes.ERROR_IONIC_NOT_CONVERGED])  # pylint: disable=no-member
    def handle_unconverged(self, node: orm.WorkChainNode) -> ProcessHandlerReport | None:
        """
        Handle the problem where the NEB optimization is not converged.

        Note that VASP could reach NSW before the actual convergence.
        Hence this check is necessary even for finished runs.
        """
        if 'misc' not in node.outputs:
            self.report('Cannot found the `misc` output containing the NEB run data')
            return None
        misc_dict = node.outputs.misc.get_dict()

        neb_data = misc_dict.get('neb_data')
        if neb_data is None:
            self.report('Cannot found the `neb_data` dictionary containing the NEB run data')
            return None

        converged = [tmp.get('neb_converged', False) for tmp in neb_data.values()]
        if not all(converged):
            self.report('At least one image is not converged in the run. Restart required.')

            # Attach images
            out = self._attach_output_structure(node)
            if out is not None:
                return out

            self.report(f'Successfully handled unconverged calculation {node}.')
            return ProcessHandlerReport()
        self.report(f'Cannot handle ionically unconverged calculation {node}.')
        return None

    @process_handler(priority=900, exit_codes=[VaspNEBCalculation.exit_codes.ERROR_DID_NOT_FINISH])  # pylint: disable=no-member
    def handle_unfinished(self, node: orm.WorkChainNode) -> ProcessHandlerReport | None:
        """
        Handle the case where the calculations is not fully finished.
        This checks the existing of the run_stats field in the parsed per-image misc output
        """

        finished = []
        # Since 1.6.3 the nested namespaces are handled properly.
        if version.parse(aiida_version) >= version.parse('1.6.3'):
            if 'misc' not in node.outputs:
                self.report('Cannot found the `misc` output containing the parsed per-image data')
                return None

            misc_dict = node.outputs.misc.get_dict()
            if 'run_status' in misc_dict:
                finished = {key: value.get('finished', False) for key, value in misc_dict['run_status'].items()}

        if not all(finished.values()):
            self.report('At least one image did not reach the end of VASP execution - calculation not finished!')

            out = self._attach_output_structure(node)
            if out is not None:
                return out

            # No further process handling is needed
            self.report(f'Successfully handled unfinished calculation {node}.')
            return ProcessHandlerReport(do_break=True)
        self.report(f'Cannot handle unfinished calculation {node}.')
        return None

    def _attach_output_structure(self, node: orm.WorkChainNode) -> ProcessHandlerReport | None:
        """
        Attached the output structure of a children node as the inputs for the
        next workchain launch.
        """
        output_images = AttributeDict()  # A dictionary holding the structures with keys like 'image_xx'
        output_images = node.outputs['structure']

        nout = len(output_images)
        nexists = len(self.inputs.neb_images)
        if nout != nexists:
            self.report(f'Number of parsed images: {nout} does not equal to the images need to restart: {nexists}.')
            return ProcessHandlerReport(do_break=True, exit_code=self.exit_codes.SUB_NEB_CALCULATION_ERROR)  # pylint: disable=no-member
        self.report(f'Attached output structures from the previous calculation {node} as new inputs.')
        self.ctx.inputs.neb_images = output_images
        return None

    def _check_neb_inputs(self) -> None:
        """
        Perform some simple checks for the NEB inputs

        This method is called once by ``self.setup``
        """

        incar = self.ctx.inputs.parameters

        images = incar.get('images')

        if not images:
            raise InputValidationError('IMAGES parameters is not set in the INCAR inputs')

        nimages = len(self.ctx.inputs.neb_images)

        if nimages != images:
            raise InputValidationError('Mismatch between IMAGES and actual number supplied input structures.')

        # Check for NEB tags
        iopt = incar.get('iopt', 0)
        ibrion = incar.get('ibrion')
        potim = incar.get('potim')

        # Check the sanity of parameters
        if ibrion != 3:
            self.report('WARNING: IBRION should be set to 3 for VTST runs, proceed with caution.')
        elif potim != 0:
            self.report(
                'WARNING: Using VTST optimisors with IBRION=3, but POTIM is not set to zero, proceed with caution.'
            )
        if iopt == 0:
            self.report('WARNING: IOPT not set.')

        if ibrion == 2:
            raise InputValidationError('IBRION=2 should not be used for NEB optimization!!')

        # Check the displacement of atoms between the frames
        # the hope is that this may detect simple errors such as atoms going across the PBC or
        # the order of atoms are changed between different frames

        tmp = list(self.ctx.inputs.neb_images.items())
        tmp.sort(key=lambda x: x[0])
        frames = [x[1].get_ase() for x in tmp]
        frames = [self.ctx.inputs.initial_structure.get_ase()] + frames + [self.ctx.inputs.final_structure.get_ase()]

        last_frame = frames[0]
        # Function for computing the distance using the scaled positions
        rel_dist = np.vectorize(lambda x: x if x < 0.5 else 1.0 - x)
        for iframe, frame in enumerate(frames[1:]):
            # Relative displacements
            disp = abs(frame.get_scaled_positions() - last_frame.get_scaled_positions()) % 1.0
            # Apply convention
            disp = rel_dist(disp)
            # Convert back to absolute displacement
            disp = disp @ frame.cell
            norm_disp = np.linalg.norm(disp, axis=1)
            sort_idx = np.argsort(norm_disp)
            if norm_disp[sort_idx[-1]] > self._norm_disp_threshold:
                raise InputValidationError(
                    'Large displacement detected for atom {} at frame {} - please check the inputs images'.format(
                        sort_idx[-1], iframe + 1
                    )
                )
            last_frame = frame

    def _setup_vasp_inputs(self) -> None:
        """
        Setup the inputs for VASP calculation

        This method is called once by ``self.setup``
        TODO: merge with vasp.v2.vasp
        """

        # Set the kpoints (kpoints)
        if 'kpoints' in self.inputs:
            self.ctx.inputs.kpoints = self.inputs.kpoints
        elif 'kpoints_spacing' in self.inputs:
            kpoints = orm.KpointsData()
            kpoints.set_cell_from_structure(self.ctx.inputs.initial_structure)
            kpoints.set_kpoints_mesh_from_density(self.inputs.kpoints_spacing.value * np.pi * 2)
            self.ctx.inputs.kpoints = kpoints
        elif 'kpoints_spacing_vasp' in self.inputs:
            kpoints = orm.KpointsData()
            kpoints.set_cell_from_structure(self.ctx.inputs.initial_structure)
            kpoints.set_kpoints_mesh_from_density(self.inputs.kpoints_spacing.value)
            self.ctx.inputs.kpoints = kpoints
        else:
            raise InputValidationError("Must supply either 'kpoints' or 'kpoints_spacing' or 'kpoints_spacing_vasp")

        # Set settings

        unsupported_parameters = dict(VTST_ADDITIONAL_TAGS)
        if 'settings' in self.inputs:
            self.ctx.inputs.settings = self.inputs.settings
            # Also check if the user supplied additional tags that is not in the supported file.
            try:
                unsupported_parameters = self.ctx.inputs.settings.unsupported_parameters
            except AttributeError:
                pass

        # Perform inputs massage to accommodate generalization in higher lying workchains
        # and set parameters.
        try:
            parameters_massager = ParametersMassage(self.inputs.parameters, unsupported_parameters)
        except Exception as exception:  # pylint: disable=broad-except
            return self.exit_codes.ERROR_IN_PARAMETER_MASSAGER.format(exception=exception)  # pylint: disable=no-member
        try:
            # Only set if they exists
            # Set any INCAR tags
            self.ctx.inputs.parameters = parameters_massager.parameters.incar
            # Set any dynamics input (currently only for selective dynamics, e.g. custom write to POSCAR)
            self.ctx.inputs.dynamics = parameters_massager.parameters.dynamics
            # Here we could set additional override flags, but those are not relevant for this VASP plugin
        except AttributeError:
            pass

        # Setup LDAU keys
        if 'ldau_mapping' in self.inputs:
            ldau_settings = self.inputs.ldau_mapping.get_dict()
            ldau_keys = get_ldau_keys(self.ctx.inputs.initial_structure, **ldau_settings)
            # Directly update the raw inputs passed to VaspCalculation
            self.ctx.inputs.parameters.update(ldau_keys)

        # Set settings
        if 'settings' in self.inputs:
            self.ctx.inputs.settings = self.inputs.settings

        # Set options
        # Options is very special, not storable and should be
        # wrapped in the metadata dictionary, which is also not storable
        # and should contain an entry for options
        if 'options' in self.inputs:
            options = {}
            options.update(self.inputs.options)
            self.ctx.inputs.metadata = {}
            self.ctx.inputs.metadata['options'] = options
            # Override the parser name if it is supplied by the user.
            parser_name = self.ctx.inputs.metadata['options'].get('parser_name')
            if parser_name:
                self.ctx.inputs.metadata['options']['parser_name'] = parser_name
            # Also make sure we specify the entry point for the
            # Set MPI to True, unless the user specifies otherwise
            withmpi = self.ctx.inputs.metadata['options'].get('withmpi', True)
            self.ctx.inputs.metadata['options']['withmpi'] = withmpi
        else:
            raise InputValidationError('`options` not supplied')

        # Utilise default input/output selections
        self.ctx.inputs.metadata['options']['input_filename'] = 'INCAR'

        # Set the CalcJobNode to have the same label as the WorkChain
        self.ctx.inputs.metadata['label'] = self.inputs.metadata.get('label', '')
        self.report(self.ctx.inputs.metadata)

        # Verify and set potentials (potcar)
        if not self.inputs.potential_family.value:
            self.report(  # pylint: disable=not-callable
                'An empty string for the potential family name was detected.'
            )
            return self.exit_codes.ERROR_NO_POTENTIAL_FAMILY_NAME  # pylint: disable=no-member
        try:
            self.ctx.inputs.potential = PotcarData.get_potcars_from_structure(
                structure=self.inputs.initial_structure,
                family_name=self.inputs.potential_family.value,
                mapping=self.inputs.potential_mapping.get_dict(),
            )
        except ValueError as err:
            return compose_exit_code(self.exit_codes.ERROR_POTENTIAL_VALUE_ERROR.status, str(err))  # pylint: disable=no-member
        except NotExistent as err:
            return compose_exit_code(self.exit_codes.ERROR_POTENTIAL_DO_NOT_EXIST.status, str(err))  # pylint: disable=no-member

        self.ctx.verbose = bool(self.inputs.get('verbose', self._verbose))

        return None

    @override
    def results(self) -> ExitCode | None:
        """Attach the outputs specified in the output specification from the last completed process."""
        node = self.ctx.children[self.ctx.iteration - 1]

        # We check the `is_finished` attribute of the work chain and not the successfulness of the last process
        # because the error handlers in the last iteration can have qualified a "failed" process as satisfactory
        # for the outcome of the work chain and so have marked it as `is_finished=True`.
        max_iterations = self.inputs.max_iterations.value  # type: ignore[union-attr]
        if not self.ctx.is_finished and self.ctx.iteration >= max_iterations:
            self.report(
                f'reached the maximum number of iterations {max_iterations}: '
                f'last ran {self.ctx.process_name}<{node.pk}>'
            )
            return self.exit_codes.ERROR_MAXIMUM_ITERATIONS_EXCEEDED  # pylint: disable=no-member

        self.report(f'work chain completed after {self.ctx.iteration} iterations')

        # Simply attach the output of the last children
        self.out_many({key: node.outputs[key] for key in node.outputs})
        return None


# The code below should be moved for utility module, but I keep them here for now

FELEMS = [
    'La',
    'Ce',
    'Pr',
    'Nd',
    'Pm',
    'Sm',
    'Eu',
    'Gd',
    'Tb',
    'Dy',
    'Ho',
    'Er',
    'Tm',
    'Yb',
    'Lu',
    'Ac',
    'Th',
    'Pa',
    'U',
    'Np',
    'Pu',
    'Am',
    'Cm',
    'Bk',
    'Cf',
    'Es',
    'Fm',
    'Md',
    'No',
    'Lr',
]


def get_ldau_keys(
    structure: orm.StructureData,
    mapping: dict[str, list[float]],
    utype: int = 2,
    jmapping: dict[str, float] | None = None,
    felec: bool = False,
) -> dict[str, Any]:
    """
    Setup LDAU mapping.

    In VASP, the U for each species has to be defined in the order that they appear in POSCAR.
    This is a helper to make sure the values of U are associated to each species.

    :param structure: The structure, either StructureData or ase.Atoms is fine.
    :param mapping: A dictionary in the format of  ``{"Mn": [d, 4]...}`` for U.
    :param utype: The type of LDA+U, default to 2, which is the one with only one parameter.
    :param jmapping: A dictionary in the format of  ``{"Mn": [d, 4]...}`` but for J.
    :param felec: Whether we are dealing with f electrons; will increase lmaxmix if we are.

    :returns: A dictionary to be used to update the raw input parameters for VASP.
    """
    if isinstance(structure, orm.StructureData):
        species = MultiPotcarIo.potentials_order(structure)
    else:
        # For ASE atoms, we keep the order of species occurrence no sorting is done
        species = []
        for symbol in structure.get_chemical_symbols():
            if symbol not in species:
                species.append(symbol)

    lsymbols = {'d': 2, 'f': 3, 'p': 1}
    if jmapping is None:
        jmapping = {}
    # Setup U array
    ldauu = []
    ldauj = []
    ldaul = []
    count = 0
    for specie in species:
        if specie in mapping:
            uvalue = mapping[specie][1]
            j = jmapping.get(specie, 0.0)
            ldaul.append(lsymbols[mapping[specie][0]])
            ldauu.append(mapping[specie][1])

            j = jmapping.get(specie, 0.0)
            ldauj.append(j)

            if specie in FELEMS:
                felec = True
            # Count the number of valid mappings
            if uvalue != 0.0 or j != 0.0:
                count += 1

        else:
            ldauu.append(0.0)
            ldauj.append(0.0)
            ldaul.append(-1)

    if count > 0:
        # Only enable U is there is any non-zero value
        output = {
            'ldauu': ldauu,
            'ldauj': ldauj,
            'ldautype': utype,
            'lmaxmix': 6 if felec else 4,
            'ldaul': ldaul,
            'ldau': True,
        }
    else:
        output = {}
    return output
