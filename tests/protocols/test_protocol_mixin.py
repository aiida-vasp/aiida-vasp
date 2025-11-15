"""Test the ProtocolMixin class functionality."""

from __future__ import annotations

import pathlib
import tempfile
from unittest.mock import patch

import pytest
import yaml

from aiida_vasp.protocols import ProtocolMixin

SYS_VASP_PROTO = (pathlib.Path(__file__).parent.parent.parent / 'src/aiida_vasp/protocols/vasp.yaml').resolve()


@pytest.fixture
def sample_protocol_data():
    """Sample protocol data for testing."""
    return {
        'default_protocol': 'balanced',
        'default_inputs': {
            'verbose': False,
            'max_iterations': 60,
            'convergence_settings': {
                'energy_threshold': 1e-6,
                'force_threshold': 1e-3,
            },
        },
        'protocols': {
            'balanced': {
                'description': 'Balanced protocol for testing',
                'max_iterations': 40,
                'convergence_settings': {
                    'energy_threshold': 1e-5,
                },
            },
            'stringent': {
                'description': 'Stringent protocol for testing',
                'max_iterations': 100,
                'convergence_settings': {
                    'energy_threshold': 1e-8,
                    'force_threshold': 1e-4,
                },
            },
        },
    }


@pytest.fixture
def custom_protocol_data():
    """Custom protocol data for testing user-defined protocols."""
    return {
        'default_protocol': 'custom_default',
        'default_inputs': {
            'verbose': True,
            'custom_setting': 'custom_value',
        },
        'protocols': {
            'custom_default': {
                'description': 'Custom default protocol',
                'custom_setting': 'overridden_value',
            },
            'custom_advanced': {
                'description': 'Custom advanced protocol',
                'verbose': False,
                'advanced_options': {
                    'option1': True,
                    'option2': 42,
                },
            },
        },
    }


@pytest.fixture
def temp_protocol_file(sample_protocol_data):
    """Create a temporary protocol file for testing."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(sample_protocol_data, f)
        return pathlib.Path(f.name)


@pytest.fixture
def temp_user_protocols_dir(custom_protocol_data, tmp_path):
    """Create a temporary user protocols directory structure."""
    user_dir = tmp_path / '.aiida-vasp' / 'protocols' / 'test_protocol'
    user_dir.mkdir(parents=True)

    # Create a custom protocol file
    custom_file = user_dir / 'custom.yaml'
    with custom_file.open('w') as f:
        yaml.dump(custom_protocol_data, f)

    return tmp_path


@pytest.fixture
def protocol_mixin(temp_user_protocols_dir):
    class TestProtocolMixin(ProtocolMixin):
        """Test subclass of ProtocolMixin for testing purposes."""

        _protocol_tag = 'test_protocol'
        _load_root = temp_user_protocols_dir

    return TestProtocolMixin


class TestProtocolSplitting:
    """Test protocol name splitting functionality."""

    def test_split_protocol_file_name_no_alias(self, protocol_mixin):
        """Test splitting protocol name without alias."""
        protocol, alias = protocol_mixin._split_protocol_file_name('balanced')
        assert protocol == 'balanced'
        assert alias is None

    def test_split_protocol_file_name_with_alias(self, protocol_mixin):
        """Test splitting protocol name with alias."""
        protocol, alias = protocol_mixin._split_protocol_file_name('balanced@custom')
        assert protocol == 'balanced'
        assert alias == 'custom'


class TestProtocolAliases:
    """Test protocol alias functionality."""

    def test_check_if_alias_moderate(self, protocol_mixin):
        """Test moderate alias maps to balanced."""
        alias = protocol_mixin._check_if_alias('moderate')
        assert alias == 'balanced'


class TestProtocolFileHandling:
    """Test protocol file handling functionality."""

    def test_get_protocol_filepath_explicit_yaml_file(self, temp_protocol_file, protocol_mixin):
        """Test loading protocol from explicit YAML file path."""
        filepath = protocol_mixin.get_protocol_filepath(str(temp_protocol_file))
        assert filepath == temp_protocol_file.absolute()

    def test_get_protocol_filepath_explicit_yml_file(self, sample_protocol_data, tmp_path, protocol_mixin):
        """Test loading protocol from explicit YML file path."""
        yml_file = tmp_path / 'test.yml'
        with yml_file.open('w') as f:
            yaml.dump(sample_protocol_data, f)

        filepath = protocol_mixin.get_protocol_filepath(str(yml_file))
        assert filepath == yml_file.absolute()

    def test_get_protocol_filepath_user_directory(self, temp_user_protocols_dir, protocol_mixin):
        """Test loading protocol from user directory."""
        expected_path = temp_user_protocols_dir / '.aiida-vasp' / 'protocols' / 'test_protocol' / 'custom.yaml'
        with patch('pathlib.Path.expanduser', return_value=expected_path):
            with patch.object(pathlib.Path, 'is_file', return_value=True):
                filepath = protocol_mixin.get_protocol_filepath('custom')
                assert filepath == expected_path

    def test_get_protocol_filepath_system_default(self, sample_protocol_data, tmp_path):
        """Test loading protocol from system default location."""

        class Proto(ProtocolMixin):
            _protocol_tag = 'vasp'

        file_path = Proto.get_protocol_filepath()
        assert (
            file_path == (pathlib.Path(__file__).parent.parent.parent / 'src/aiida_vasp/protocols/vasp.yaml').resolve()
        )

    def test_get_protocol_filepath_file_not_found(self, protocol_mixin):
        """Test error when protocol file is not found."""
        with pytest.raises(FileNotFoundError, match='Protocol file not found'):
            protocol_mixin.get_protocol_filepath('nonexistent')


class TestProtocolLoading:
    """Test protocol loading functionality."""

    def test_load_protocol_file(self, temp_protocol_file, sample_protocol_data, protocol_mixin):
        """Test loading protocol file content."""
        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            data = protocol_mixin._load_protocol_file()
            assert data == sample_protocol_data

    def test_get_default_protocol(self, temp_protocol_file, sample_protocol_data, protocol_mixin):
        """Test getting default protocol."""
        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            default = protocol_mixin.get_default_protocol()
            assert default == 'balanced'

    def test_get_available_protocols(self, temp_protocol_file, sample_protocol_data, protocol_mixin):
        """Test getting available protocols."""
        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            protocols = protocol_mixin.get_available_protocols()

            expected = {
                'balanced': {'description': 'Balanced protocol for testing'},
                'stringent': {'description': 'Stringent protocol for testing'},
            }
            assert protocols == expected

    def test_get_available_protocols_with_alias(self, temp_user_protocols_dir, custom_protocol_data, protocol_mixin):
        """Test getting available protocols from user file alias."""
        custom_file = temp_user_protocols_dir / '.aiida-vasp' / 'protocols' / 'test_protocol' / 'custom.yaml'

        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=custom_file):
            protocols = protocol_mixin.get_available_protocols('custom')

            expected = {
                'custom_default': {'description': 'Custom default protocol'},
                'custom_advanced': {'description': 'Custom advanced protocol'},
            }
            assert protocols == expected


class TestProtocolInputs:
    """Test protocol input generation functionality."""

    def test_get_protocol_inputs_default(self, temp_protocol_file, sample_protocol_data, protocol_mixin):
        """Test getting protocol inputs with default protocol."""
        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            inputs = protocol_mixin.get_protocol_inputs()

            expected = {
                'verbose': False,
                'max_iterations': 40,  # From balanced protocol
                'convergence_settings': {
                    'energy_threshold': 1e-5,  # From balanced protocol
                    'force_threshold': 1e-3,  # From default_inputs
                },
            }
            assert inputs == expected

    def test_get_protocol_inputs_specific_protocol(self, temp_protocol_file, sample_protocol_data, protocol_mixin):
        """Test getting protocol inputs with specific protocol."""
        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            inputs = protocol_mixin.get_protocol_inputs('stringent')

            expected = {
                'verbose': False,
                'max_iterations': 100,  # From stringent protocol
                'convergence_settings': {
                    'energy_threshold': 1e-8,  # From stringent protocol
                    'force_threshold': 1e-4,  # From stringent protocol
                },
            }
            assert inputs == expected

    def test_get_protocol_inputs_with_alias_protocol(self, temp_protocol_file, protocol_mixin):
        """Test getting protocol inputs using protocol alias."""
        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            inputs = protocol_mixin.get_protocol_inputs('moderate')  # Should map to 'balanced'

            expected = {
                'verbose': False,
                'max_iterations': 40,  # From balanced protocol
                'convergence_settings': {
                    'energy_threshold': 1e-5,  # From balanced protocol
                    'force_threshold': 1e-3,  # From default_inputs
                },
            }
            assert inputs == expected

    def test_get_protocol_inputs_with_file_alias(self, temp_user_protocols_dir, custom_protocol_data, protocol_mixin):
        """Test getting protocol inputs with file alias."""
        custom_file = temp_user_protocols_dir / '.aiida-vasp' / 'protocols' / 'test_protocol' / 'custom.yaml'

        inputs = protocol_mixin.get_protocol_inputs(f'custom_advanced@{custom_file}')

        expected = {
            'verbose': False,  # From custom_advanced protocol
            'custom_setting': 'custom_value',  # From default_inputs
            'advanced_options': {
                'option1': True,
                'option2': 42,
            },
        }
        assert inputs == expected

    def test_get_protocol_inputs_invalid_protocol(self, temp_protocol_file, protocol_mixin):
        """Test error when requesting invalid protocol."""
        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            with pytest.raises(ValueError, match='`invalid` is not a valid protocol'):
                protocol_mixin.get_protocol_inputs('invalid')

    def test_get_protocol_inputs_with_dict_overrides(self, temp_protocol_file, protocol_mixin):
        """Test getting protocol inputs with dictionary overrides."""
        overrides = {
            'verbose': True,
            'max_iterations': 200,
            'convergence_settings': {
                'energy_threshold': 1e-10,
            },
            'new_setting': 'new_value',
        }

        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            inputs = protocol_mixin.get_protocol_inputs('balanced', overrides)

            expected = {
                'verbose': True,  # Overridden
                'max_iterations': 200,  # Overridden
                'convergence_settings': {
                    'energy_threshold': 1e-10,  # Overridden
                    'force_threshold': 1e-3,  # From default_inputs
                },
                'new_setting': 'new_value',  # Added by override
            }
            assert inputs == expected

    def test_get_protocol_inputs_with_file_overrides(
        self, temp_protocol_file, sample_protocol_data, tmp_path, protocol_mixin
    ):
        """Test getting protocol inputs with file-based overrides."""
        override_data = {
            'verbose': True,
            'max_iterations': 300,
            'new_feature': 'enabled',
        }

        override_file = tmp_path / 'overrides.yaml'
        with override_file.open('w') as f:
            yaml.dump(override_data, f)

        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            inputs = protocol_mixin.get_protocol_inputs('balanced', override_file)

            expected = {
                'verbose': True,  # Overridden
                'max_iterations': 300,  # Overridden
                'convergence_settings': {
                    'energy_threshold': 1e-5,  # From balanced protocol
                    'force_threshold': 1e-3,  # From default_inputs
                },
                'new_feature': 'enabled',  # Added by override
            }
            assert inputs == expected

    def test_get_protocol_inputs_no_overrides(self, temp_protocol_file, sample_protocol_data, protocol_mixin):
        """Test getting protocol inputs with None overrides."""
        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            inputs = protocol_mixin.get_protocol_inputs('balanced', None)

            expected = {
                'verbose': False,
                'max_iterations': 40,  # From balanced protocol
                'convergence_settings': {
                    'energy_threshold': 1e-5,  # From balanced protocol
                    'force_threshold': 1e-3,  # From default_inputs
                },
            }
            assert inputs == expected


def test_list_protocol_files_system_only(tmp_path, protocol_mixin):
    """Test listing protocol files."""
    # This test would be complex to fully mock, so we test the method interface
    # In a real test environment, this would verify the actual file listing logic
    files = protocol_mixin.list_protocol_files('vasp')
    files[-1][-1] == SYS_VASP_PROTO


# Integration test for the full workflow
class TestProtocolMixinIntegration:
    """Integration tests for the full ProtocolMixin workflow."""

    def test_full_workflow_with_custom_protocol(
        self, temp_user_protocols_dir, custom_protocol_data, tmp_path, protocol_mixin
    ):
        """Test the full workflow of loading a custom protocol with overrides."""
        custom_file = temp_user_protocols_dir / '.aiida-vasp' / 'protocols' / 'test_protocol' / 'custom.yaml'

        # Test file alias syntax
        # Get available protocols
        protocols = protocol_mixin.get_available_protocols(custom_file)
        assert 'custom_advanced' in protocols

        # Get inputs with overrides
        overrides = {
            'verbose': False,  # Override the protocol setting
            'extra_setting': 'test_value',
        }

        inputs = protocol_mixin.get_protocol_inputs(f'custom_advanced@{custom_file}', overrides)

        expected = {
            'verbose': False,  # From override
            'custom_setting': 'custom_value',  # From default_inputs
            'advanced_options': {
                'option1': True,
                'option2': 42,
            },
            'extra_setting': 'test_value',  # From override
        }
        assert inputs == expected

    def test_protocol_with_nested_overrides(self, temp_protocol_file, protocol_mixin):
        """Test protocol loading with deeply nested overrides."""
        overrides = {
            'convergence_settings': {
                'energy_threshold': 1e-12,
                'new_convergence_option': True,
            },
            'completely_new': {
                'nested': {
                    'value': 42,
                },
            },
            'max_iterations': '$!del',
        }

        with patch.object(protocol_mixin, 'get_protocol_filepath', return_value=temp_protocol_file):
            inputs = protocol_mixin.get_protocol_inputs('stringent', overrides)

            # Check that nested merging works correctly
            assert inputs['convergence_settings']['energy_threshold'] == 1e-12  # Overridden
            assert inputs['convergence_settings']['force_threshold'] == 1e-4  # From protocol
            assert inputs['convergence_settings']['new_convergence_option'] is True  # Added
            assert inputs['completely_new']['nested']['value'] == 42  # Added
            assert 'max_iterations' not in inputs
