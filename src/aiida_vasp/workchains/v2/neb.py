"""
VASP NEB workchain.

Contains the VaspNEBWorkChain class definition which uses the BaseRestartWorkChain.
"""

from __future__ import annotations

import numpy as np
from aiida import orm
from aiida.common.exceptions import InputValidationError, NotExistent
from aiida.common.extendeddicts import AttributeDict
from aiida.common.lang import override
from aiida.engine import ExitCode
from aiida.engine.processes.workchains.restart import ProcessHandlerReport, process_handler
from aiida.plugins import CalculationFactory

from aiida_vasp.calcs.neb import VaspNEBCalculation
from aiida_vasp.data.potcar import PotcarData
from aiida_vasp.utils.workchains import compose_exit_code

from .vasp import VaspWorkChain

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


class VaspNEBWorkChain(VaspWorkChain):
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
    _default_unsupported_parameters = VTST_ADDITIONAL_TAGS

    def init_inputs(self) -> None | ExitCode:
        exit_code = super().init_inputs()
        if exit_code is not None:
            return exit_code
        # Run some additional checks
        return self.check_neb_inputs()

    def setup_potcar(self) -> None:
        # Verify and set potentials (potcar)
        if not self.inputs.potential_family.value:
            self.report('An empty string for the potential family name was detected.')  # pylint: disable=not-callable
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

    def check_neb_inputs(self) -> None | ExitCode:
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

    @process_handler(priority=500, exit_codes=[VaspNEBCalculation.exit_codes.ERROR_IONIC_NOT_CONVERGED])  # pylint: disable=no-member
    def handle_ionic_conv(self, node: orm.WorkChainNode) -> ProcessHandlerReport | None:
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

    @process_handler(priority=910, exit_codes=[VaspNEBCalculation.exit_codes.ERROR_DID_NOT_FINISH])  # pylint: disable=no-member
    def handle_unfinished_calc_ionic(self, node: orm.WorkChainNode) -> ProcessHandlerReport | None:
        """
        Handle the case where the calculations is not fully finished.
        This checks the existing of the run_stats field in the parsed per-image misc output
        """

        finished = []
        # Since 1.6.3 the nested namespaces are handled properly.
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

    @process_handler(priority=4)
    def check_calc_is_finished(self, node: orm.CalcJobNode) -> ProcessHandlerReport | None:
        """
        Check if the calculation has reached the end of execution.
        """
        misc = node.outputs.misc.get_dict()
        run_status = misc['run_status']
        for key, value in run_status.items():
            if not value.get('finished'):
                self.report(f'The child calculation {node} - image {key} did not reach the end of execution.')
                return ProcessHandlerReport(exit_code=self.exit_codes.ERROR_CALCULATION_NOT_FINISHED, do_break=True)
            return None

    def _get_run_status(self, node):
        """Return the run status of the calculation."""
        return node.outputs.misc['run_status']['01']
