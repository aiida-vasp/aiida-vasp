"""
Module containing the OptionHolder class
"""

from typing import Optional

from aiida.orm import Dict
from pydantic import BaseModel, Field, ValidationError

# pylint:disable=raise-missing-from


class OptionContainer(BaseModel):
    """
    Base class for a container of options
    """

    def aiida_dict(self) -> Dict:
        """Return an ``aiida.orm.Dict`` presentation"""

        python_dict = self.model_dump()
        return Dict(dict=python_dict)

    @classmethod
    def aiida_validate(cls, input_dict: dict | Dict, namespace: None | str = None) -> None:  # pylint:disable=unused-argument
        """
        Validate a dictionary/Dict node, this can be used as the validator for
        the Port accepting the inputs

        This is used as validator for the `spec.input` call.
        """
        if isinstance(input_dict, Dict):
            input_dict = input_dict.get_dict()
        try:
            cls(**input_dict)
        except ValidationError as error:
            return str(error)
        return None

    @classmethod
    def aiida_serialize(cls, python_dict: dict) -> Dict:
        """
        serialize a dictionary into Dict

        This method can be passed as a `serializer` key word parameter of for the `spec.input` call.
        """
        obj = cls(**python_dict)
        return obj.aiida_dict()

    @classmethod
    def aiida_description(cls) -> str:
        """
        Return a string for the options of a OptionContains in a human-readable format.
        """

        obj = cls()
        template = '{:>{width_name}s}:  {:10s} \n{default:>{width_name2}}: {}'
        entries = []
        for name, field in obj.model_fields.items():
            # Each entry is name, type, doc, default value
            entries.append([name, str(field.annotation.__name__), field.description, field.default])
        max_width_name = max(len(entry[0]) for entry in entries) + 2

        lines = []
        for entry in entries:
            lines.append(
                template.format(
                    *entry,
                    width_name=max_width_name,
                    width_name2=max_width_name + 10,
                    default='Default',
                )
            )
        return '\n'.join(lines)


class CalcSettingsConfig(OptionContainer):
    """Schema for the .settings port used by both VaspCalculation and VaspWorkChain"""

    parser_setting: Optional[dict] = Field(description='Settings for the parser')
    ADDITIONAL_RETRIEVE_LIST: Optional[list] = Field(description='Additional list of files to be retrieved')
    ADDITIONAL_RETRIEVE_TEMPORARY_LIST: Optional[list] = Field(
        description=('Additional list of files to be retrieved, but not store in the storage')
    )
    PROVENANCE_EXCLUDE_LIST: Optional[list] = Field(
        description=('Additional list of files to be retrieved, but not store in the storage')
    )
    ALWAYS_STORE: Optional[list] = Field(
        description=('Additional list of files to be retrieved, but not store in the storage')
    )
    skip_param_validation: Optional[bool] = Field(
        description='Skip the validation of the input parameters', default=False
    )
    unsupported_parameters: Optional[dict] = Field(description='None-standard VASP parameters that are valid')


class RelaxOptions(OptionContainer):
    """Options for VaspRelaxWorkChain"""

    algo: str = Field(description='The algorithm to use for relaxation', examples=['cg', 'rd'], default='cg')
    energy_cutoff: Optional[float] = Field(
        description='The cut off energy difference when the relaxation is stopped (e.g. EDIFF)',
        default=None,
    )
    force_cutoff: float = Field(
        description='The maximum force when the relaxation is stopped (e.g. EDIFFG)',
        default=0.03,
    )
    steps: int = Field(description='Number of relaxation steps to perform (eg. NSW)', default=60)
    positions: bool = Field(description='If True, perform relaxation of the atomic positions', default=True)
    shape: bool = Field(description='If True, perform relaxation of the cell shape', default=True)
    volume: bool = Field(description='If True, perform relaxation of the cell volume', default=True)
    convergence_on: bool = Field(description='If True, perform convergence checks within the workchain', default=True)
    convergence_absolute: bool = Field(
        description='If True, use absolute values where possible when performing convergence checks',
        default=False,
    )
    convergence_max_iterations: int = Field(description='Maximum iterations for convergence checking', default=5)
    convergence_positions: float = Field(
        description=(
            'The cutoff value for the convergence check on positions in Angstram. A negative value by pass the check.'
        ),
        default=0.1,
    )
    convergence_volume: float = Field(
        description='The cutoff value for the convergence check on volume between the two structures.'
        ' A negative value by pass the check.',
        default=0.01,
    )
    convergence_shape_lengths: float = Field(
        description='The cutoff value for the convergence check on the lengths of the unit cell'
        ' vectors, between input and the outputs structure. A negative value by pass'
        ' the check.',
        default=0.1,
    )
    convergence_shape_angles: float = Field(
        description='The cutoff value for the convergence check on the angles of the unit cell vectors, '
        'between input and the outputs structure. A'
        ' negative value by pass the check.',
        default=0.1,
    )
    convergence_mode: str = Field(
        description="Mode of the convergence check for positions. 'inout' for checking input/output structure, "
        "or 'last' to check only the change of"
        ' the last step.',
        examples=['inout', 'last'],
        default='last',
    )
    reuse: bool = Field(
        description='Whether reuse the previous calculation by copying over the remote folder',
        default=False,
    )
    clean_reuse: bool = Field(
        description='Whether to perform a final cleaning of the reused calculations', default=True
    )
    keep_sp_workdir: bool = Field(
        description='Whether to keep the workdir of the final singlepoint calculation', default=False
    )
    perform: bool = Field(description="Do not perform any relaxation if set to 'False'", default=True)
    hybrid_calc_bootstrap: bool = Field(
        description='Whether to bootstrap hybrid calculation by performing standard DFT first', default=False
    )
    hybrid_calc_bootstrap_wallclock: int = Field(
        description='Wall time limit in second for the bootstrap calculation', default=3600
    )
    keep_magnetization: bool = Field(
        description='Whether to keep magnetization from the previous calculation if possible', default=False
    )
    double_relax_mode: bool = Field(
        description='Experimental: Run in double relax mode - launch of the sub workflow is only performed up to two '
        'times without checking convergence in the end. This is useful for cases where the convergence is difficult '
        'due to change of basis set with variable cell and high-throughput studies.',
        default=False,
    )

    # TODO: implement pyandtic checks

    # @classmethod
    # def validate_dict(cls, input_dict, port=None):
    #     """Check mutually exclusive fields"""
    #     super().validate_dict(input_dict, port)
    #     if isinstance(input_dict, orm.Dict):
    #         input_dict = input_dict.get_dict()
    #     force_cut = input_dict.get('force_cutoff')
    #     energy_cut = input_dict.get('energy_cutoff')
    #     if force_cut is None and energy_cut is None:
    #         raise InputValidationError("Either 'force_cutoff' or 'energy_cutoff' should be supplied")
    #     if (force_cut is not None) and (energy_cut is not None):
    #         raise InputValidationError("Cannot set both 'force_cutoff' and 'energy_cutoff'")


class ConvOptions(OptionContainer):
    """Template for the Dict node controlling the workchain behaviour"""

    cutoff_start: float = Field(description='The starting cut-off energy', default=300.0)
    cutoff_stop: float = Field(description='The Final cut-off energy', default=700.0)
    cutoff_step: float = Field(description='Step size of the cut-off energy', default=50.0)
    kspacing_start: float = Field(description='The starting kspacing', default=0.07)
    kspacing_stop: float = Field(description='The final kspacing', default=0.02)
    kspacing_step: float = Field(description='Step size of the cut-off energy', default=0.01)
    cutoff_kconv: float = Field(description='The cut-off energy used for kpoints convergence tests', default=450.0)
    kspacing_cutconv: float = Field(
        description='The kpoints spacing used for cut-off energy convergence tests', default=0.07
    )


class BandOptions(OptionContainer):
    """Options for VaspRelaxWorkChain"""

    symprec: float = Field(description='Precision of the symmetry determination', default=0.01)
    band_mode: str = Field(
        description=(
            'Mode for generating the band path. Choose from: bradcrack, pymatgen,seekpath-aiida and latimer-munro.'
        ),
        examples=['bradcrack', 'pymatgen', 'seekpath', 'seekpath-aiida', 'latimer-munro'],
        default='seekpath-aiida',
    )
    # TODO: enable explicit seekpath passing
    band_kpoints_distance: float = Field(
        description='Spacing for band distances for automatic kpoints generation, used by seekpath-aiida mode.',
        default=0.025,
    )
    line_density: float = Field(
        description='Density of the point along the path, used by the sumo interface.',
        default=20,
    )
    dos_kpoints_distance: float = Field(
        description=(
            'Kpoints for running DOS calculations in A^-1 * 2pi. Will perform non-SCF DOS calculation is supplied.'
        ),
        default=0.03,
    )
    only_dos: bool = Field(
        description='Flag for running only DOS calculations',
        default=False,
    )
    run_dos: bool = Field(
        description='Flag for running DOS calculations',
        default=False,
    )
    additional_band_analysis_parameters: dict = Field(
        description='Additional keyword arguments for the seekpath/ interface, available keys are:'
        "  ['with_time_reversal', 'reference_distance', 'recipe', 'threshold', 'symprec', 'angle_tolerance']",
        default={},
    )
    kpoints_per_split: int = Field(
        description='Number of kpoints per split for the band structure calculation',
        default=80,
    )
    hybrid_reuse_wavecar: bool = Field(
        description='Whether to reuse the WAVECAR from the previous relax/singlepoint calculation', default=False
    )
