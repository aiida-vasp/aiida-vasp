"""
Module for storing protocols for AiiDA VASP workflows.
"""

from __future__ import annotations

import pathlib

import yaml

from aiida_vasp.utils.dict_merge import recursive_merge

from .generator import *


class ProtocolMixin:
    """Utility class for processes to build input mappings for a given protocol based on a YAML configuration file."""

    _protocol_tag: str = 'NULL'
    _load_root: str = '~/.aiida-vasp/protocols'

    @staticmethod
    def _split_protocol_file_name(name):
        """
        Split the protocol name into its components.
        For example, "balance@my_protocol" becomes ("balance", "my_protocol").
        This allow the protocol to be loaded from a user define file, e.g ~/.aiida_vasp/relax/my_protocol.yaml
        """
        parts = name.split('@', maxsplit=1)
        if len(parts) == 1:
            return name, None
        return parts

    @classmethod
    def list_protocol_files(cls, protocol_tag=None) -> list[tuple[str | None, str, pathlib.Path]]:
        """List avaliable protocols"""

        protocol_tag = protocol_tag or '*'
        user_path = pathlib.Path(f'{cls._load_root}/{protocol_tag}').expanduser()
        system_path = pathlib.Path(__file__).parent.parent / 'protocols'

        user_files = []
        system_files = []
        for user_file in user_path.glob('*.yaml'):
            alias = user_file.stem
            tag = user_file.parent.stem
            user_files.append((alias, tag, user_file))

        for system_file in system_path.glob(f'{protocol_tag}.yaml'):
            alias = None
            tag = system_file.stem
            system_files.append((alias, tag, system_file))

        return user_files + system_files

    @classmethod
    def get_protocol_filepath(cls, file_alias: str | None = None) -> pathlib.Path:
        """Return the ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        # If user has custom defined protocols, use them as default
        assert cls._protocol_tag != 'NULL', 'Protocol tag must be set before calling this method.'
        # Use the default name
        if file_alias is None:
            file_alias = cls._protocol_tag
        else:
            file_alias = str(file_alias)
        # Return the path if it points to a file
        if (file_alias.endswith('.yaml') or file_alias.endswith('.yml')) and pathlib.Path(file_alias).is_file():
            return pathlib.Path(file_alias).absolute()
        # Check if the alias refers to a custom defined protocol file
        user_path = pathlib.Path(f'{cls._load_root}/{cls._protocol_tag}/{file_alias}.yaml').expanduser()
        if user_path.is_file():
            return user_path
        # Load the default protocol
        default_path = pathlib.Path(__file__).parent.parent / f'protocols/{cls._protocol_tag}.yaml'
        if not default_path.exists():
            raise FileNotFoundError(f'Protocol file not found at {default_path}. Please ensure it exists.')
        return default_path

    @classmethod
    def get_default_protocol(cls) -> str:
        """Return the default protocol for a given workflow class.

        :param cls: the workflow class.
        :return: the default protocol.
        """
        return cls._load_protocol_file()['default_protocol']

    @classmethod
    def get_available_protocols(cls, file_alias=None) -> dict:
        """Return the available protocols for a given workflow class.

        :param cls: the workflow class.
        :return: dictionary of available protocols, where each key is a protocol and value is another dictionary that
            contains at least the key `description` and optionally other keys with supplementary information.
        """
        data = cls._load_protocol_file(file_alias)
        return {protocol: {'description': values['description']} for protocol, values in data['protocols'].items()}

    @classmethod
    def get_protocol_inputs(
        cls,
        protocol: str | None = None,
        overrides: dict | pathlib.Path | None = None,
    ) -> dict:
        """Return the inputs for the given workflow class and protocol.

        :param cls: the workflow class.
        :param protocol: optional specific protocol, if not specified, the default will be used. An '@' symbol can be
          added to indicate which file to load the protocol from. For example, 'balanced@phonon' will load the protocol
          from '~/.aiida-vasp/cls._protocol_tag/phonon.yaml'
        :param overrides: dictionary of inputs that should override those specified by the protocol. The mapping should
            maintain the exact same nesting structure as the input port namespace of the corresponding workflow class.
        :return: mapping of inputs to be used for the workflow class.
        """
        if protocol is None:
            data = cls._load_protocol_file()
            protocol = data['default_protocol']
        else:
            protocol_name, file_alias = cls._split_protocol_file_name(protocol)
            data = cls._load_protocol_file(file_alias)
            protocol = protocol_name or data['default_protocol']

        try:
            protocol_inputs = data['protocols'][protocol]
        except KeyError as exception:
            alias_protocol = cls._check_if_alias(protocol)
            if alias_protocol is not None:
                protocol_inputs = data['protocols'][alias_protocol]
            else:
                raise ValueError(
                    f'`{protocol}` is not a valid protocol. Call ``get_available_protocols`` to show available '
                    'protocols.'
                ) from exception
        inputs = recursive_merge(data['default_inputs'], protocol_inputs)
        inputs.pop('description')

        if isinstance(overrides, pathlib.Path):
            with overrides.open() as file:
                overrides = yaml.safe_load(file)

        if overrides:
            return recursive_merge(inputs, overrides)

        return inputs

    @classmethod
    def _load_protocol_file(cls, file_alias=None) -> dict:
        """Return the contents of the protocol file for workflow class."""
        with cls.get_protocol_filepath(file_alias).open() as file:
            return yaml.safe_load(file)

    @staticmethod
    def _check_if_alias(alias: str):
        """Check if a given alias corresponds to a valid protocol."""
        aliases_dict = {
            'moderate': 'balanced',
            'precise': 'stringent',
        }
        return aliases_dict.get(alias, None)
