"""Test the generator.py utility functions."""

from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

from aiida import orm
from aiida.engine.processes.builder import ProcessBuilderNamespace

from aiida_vasp.protocols.generator import (
    get_library_path,
    has_content,
    incar_dict_to_relax_settings,
    list_protocol_presets,
    recursive_search_dict_with_key,
    recursive_search_port_basename,
    update_dict_node,
)


class TestGetLibraryPath:
    """Test the get_library_path function."""

    def test_get_library_path(self):
        """Test that get_library_path returns the correct path."""
        path = get_library_path()
        assert isinstance(path, Path)
        assert path.name == 'presets'
        assert 'protocols' in str(path)


class TestListProtocolPresets:
    """Test the list_protocol_presets function."""

    def test_list_protocol_presets_empty_directories(self):
        """Test list_protocol_presets with no preset files."""
        with patch('aiida_vasp.protocols.generator.get_library_path') as mock_get_path:
            # Create temporary directories with no files
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                mock_get_path.return_value = temp_path

                with patch('pathlib.Path.expanduser') as mock_expanduser:
                    mock_expanduser.return_value = temp_path
                    presets = list_protocol_presets()
                    assert isinstance(presets, list)
                    assert len(presets) == 0

    def test_list_protocol_presets_with_files(self):
        """Test list_protocol_presets with yaml files."""
        with patch('aiida_vasp.protocols.generator.get_library_path') as mock_get_path:
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                mock_get_path.return_value = temp_path

                # Create test yaml files
                (temp_path / 'test1.yaml').touch()
                (temp_path / 'test2.yml').touch()
                (temp_path / 'not_yaml.txt').touch()

                with patch('pathlib.Path.expanduser') as mock_expanduser:
                    mock_expanduser.return_value = Path(temp_dir) / 'nonexistent'
                    presets = list_protocol_presets()

                    assert isinstance(presets, list)
                    assert len(presets) == 2
                    assert all(isinstance(p, Path) for p in presets)
                    assert any('test1.yaml' in str(p) for p in presets)
                    assert any('test2.yml' in str(p) for p in presets)


class TestUpdateDictNode:
    """Test the update_dict_node function."""

    def test_update_dict_node_empty_content(self, aiida_profile):
        """Test update_dict_node with empty update content."""
        # Create an unstored Dict node
        original_dict = {'incar': {'nsw': 50, 'ibrion': 2}}
        node = orm.Dict(dict=original_dict)

        # Update with empty content
        update_content = {}

        # Update the node
        updated_node = update_dict_node(node, update_content, namespace='incar')

        # Should return the same node, content unchanged due to empty update
        assert updated_node is node
        expected = {'incar': {'nsw': 50, 'ibrion': 2}}
        assert updated_node.get_dict() == expected

    def test_update_dict_node_unstored_node(self, aiida_profile):
        """Test update_dict_node with unstored node."""
        # Create an unstored Dict node
        original_dict = {'incar': {'nsw': 50, 'ibrion': 2}}
        node = orm.Dict(dict=original_dict)

        # Update content
        update_content = {'nsw': 100, 'ediffg': -0.01}

        # Update the node
        updated_node = update_dict_node(node, update_content, namespace='incar')

        # Should return the same node (modified in place)
        assert updated_node is node
        expected = {'incar': {'nsw': 100, 'ibrion': 2, 'ediffg': -0.01}}
        assert updated_node.get_dict() == expected

    def test_update_dict_node_stored_node_different_content(self, aiida_profile):
        """Test update_dict_node with stored node and different content."""
        # Create and store a Dict node
        original_dict = {'incar': {'nsw': 50, 'ibrion': 2}}
        node = orm.Dict(dict=original_dict)
        node.store()

        # Update with different content
        update_content = {'nsw': 100, 'ediffg': -0.01}

        # Update the node
        updated_node = update_dict_node(node, update_content, namespace='incar')
        expected = {'incar': {'nsw': 100, 'ibrion': 2, 'ediffg': -0.01}}

        assert updated_node is not node
        assert updated_node.get_dict() == expected

    def test_update_dict_node_stored_node_same_content(self, aiida_profile):
        """Test update_dict_node with stored node and same content (reuse_if_possible=True)."""
        # Create and store a Dict node
        original_dict = {'incar': {'nsw': 50, 'ibrion': 2}}
        node = orm.Dict(dict=original_dict)
        node.store()

        # Update with same content
        update_content = {}

        # Update the node with reuse_if_possible=True
        updated_node = update_dict_node(node, update_content, namespace='incar', reuse_if_possible=True)

        # Should return the same node
        assert updated_node is node

    def test_update_dict_node_stored_node_no_reuse(self, aiida_profile):
        """Test update_dict_node with stored node and reuse_if_possible=False."""
        # Create and store a Dict node
        original_dict = {'incar': {'nsw': 50, 'ibrion': 2}}
        node = orm.Dict(dict=original_dict)
        node.store()

        # Update with content
        update_content = {'nsw': 100}

        # Update the node with reuse_if_possible=False
        updated_node = update_dict_node(node, update_content, namespace='incar', reuse_if_possible=False)

        # Should return a new node even if content might be similar
        assert updated_node is not node
        assert isinstance(updated_node, orm.Dict)

    def test_update_dict_node_no_namespace(self, aiida_profile):
        """Test update_dict_node without namespace."""
        # Create an unstored Dict node
        original_dict = {'nsw': 50, 'ibrion': 2}
        node = orm.Dict(dict=original_dict)

        # Update content without namespace
        update_content = {'nsw': 100, 'ediffg': -0.01}

        # Update the node
        updated_node = update_dict_node(node, update_content)

        # Should return the same node (modified in place)
        assert updated_node is node
        assert updated_node.get_dict()['nsw'] == 100
        assert updated_node.get_dict()['ibrion'] == 2
        assert updated_node.get_dict()['ediffg'] == -0.01

    def test_update_dict_node_nonexistent_namespace(self, aiida_profile):
        """Test update_dict_node with non-existent namespace."""
        # Create an unstored Dict node
        original_dict = {'incar': {'nsw': 50}}
        node = orm.Dict(dict=original_dict)

        # Update content with non-existent namespace
        update_content = {'new_param': 123}

        # Update the node
        updated_node = update_dict_node(node, update_content, namespace='kpoints')

        # Should return the same node (modified in place)
        assert updated_node is node
        # The non-existent namespace gets created
        expected = {'incar': {'nsw': 50}, 'kpoints': {'new_param': 123}}
        assert updated_node.get_dict() == expected


class TestIncarDictToRelaxSettings:
    """Test the incar_dict_to_relax_settings function."""

    def test_incar_dict_to_relax_settings_all_params(self):
        """Test converting all relaxation parameters."""
        incar_dict = {
            'incar': {
                'nsw': 100,
                'ibrion': 2,
                'ediffg': -0.01,
                'encut': 400,  # Should remain in incar
            }
        }

        updated_incar, relax_settings = incar_dict_to_relax_settings(incar_dict)

        # Check relax_settings
        assert relax_settings['steps'] == 100
        assert relax_settings['algo'] == 'cg'
        assert relax_settings['force_cutoff'] == -0.01

        # Check updated incar (relaxation params removed)
        assert 'nsw' not in updated_incar['incar']
        assert 'ibrion' not in updated_incar['incar']
        assert 'ediffg' not in updated_incar['incar']
        assert updated_incar['incar']['encut'] == 400

    def test_incar_dict_to_relax_settings_ibrion_rd(self):
        """Test converting ibrion=1 to RMM-DIIS algorithm."""
        incar_dict = {
            'incar': {
                'ibrion': 1,
                'encut': 400,
            }
        }

        _, relax_settings = incar_dict_to_relax_settings(incar_dict)

        # Check relax_settings
        assert relax_settings['algo'] == 'rd'
        assert 'steps' not in relax_settings
        assert 'force_cutoff' not in relax_settings

    def test_incar_dict_to_relax_settings_partial_params(self):
        """Test converting with only some relaxation parameters present."""
        incar_dict = {
            'incar': {
                'nsw': 50,
                'encut': 400,
            }
        }

        _, relax_settings = incar_dict_to_relax_settings(incar_dict)

        # Check relax_settings
        assert relax_settings['steps'] == 50
        assert 'algo' not in relax_settings
        assert 'force_cutoff' not in relax_settings

    def test_incar_dict_to_relax_settings_no_relax_params(self):
        """Test with no relaxation parameters."""
        incar_dict = {
            'incar': {
                'encut': 400,
                'ismear': 0,
            }
        }

        updated_incar, relax_settings = incar_dict_to_relax_settings(incar_dict)

        # Check relax_settings (should be empty)
        assert relax_settings == {}

        # Check updated incar (should be unchanged)
        assert updated_incar['incar']['encut'] == 400
        assert updated_incar['incar']['ismear'] == 0

    def test_incar_dict_to_relax_settings_ibrion_other_values(self):
        """Test with ibrion values other than 1 or 2."""
        incar_dict = {
            'incar': {
                'ibrion': 3,  # Not 1 or 2
                'nsw': 100,
            }
        }

        _, relax_settings = incar_dict_to_relax_settings(incar_dict)

        # Check relax_settings
        assert relax_settings['steps'] == 100
        assert 'algo' not in relax_settings  # Should not be set for ibrion != 1 or 2


class TestRecursiveSearchDictWithKey:
    """Test the recursive_search_dict_with_key function."""

    def test_recursive_search_dict_with_key_found(self, aiida_profile):
        """Test finding Dict nodes with specific key."""
        # Create mock namespace with Dict nodes
        namespace = MagicMock(spec=ProcessBuilderNamespace)

        # Create Dict nodes
        dict_with_key = orm.Dict(dict={'incar': {'nsw': 50}})
        dict_without_key = orm.Dict(dict={'other': {'value': 100}})

        # Mock the namespace structure
        namespace._valid_fields = ['parameters', 'settings', 'other']
        namespace.get.side_effect = lambda key: {
            'parameters': dict_with_key,
            'settings': dict_without_key,
            'other': 'not_a_dict',
        }[key]

        results = recursive_search_dict_with_key(namespace, 'incar')

        assert len(results) == 1
        assert results[0][0] == 'parameters'
        assert results[0][1] is dict_with_key

    def test_recursive_search_dict_with_key_nested(self, aiida_profile):
        """Test recursive search in nested namespaces."""
        # Create nested namespace structure
        sub_namespace = MagicMock(spec=ProcessBuilderNamespace)
        main_namespace = MagicMock(spec=ProcessBuilderNamespace)

        # Create Dict nodes
        dict_with_key = orm.Dict(dict={'incar': {'nsw': 50}})

        # Set up sub namespace
        sub_namespace._valid_fields = ['parameters']
        sub_namespace.get.side_effect = lambda key: {'parameters': dict_with_key}[key]

        # Set up main namespace
        main_namespace._valid_fields = ['calc']
        main_namespace.get.side_effect = lambda key: {'calc': sub_namespace}[key]

        results = recursive_search_dict_with_key(main_namespace, 'incar')

        assert len(results) == 1
        assert results[0][0] == 'calc.parameters'
        assert results[0][1] is dict_with_key

    def test_recursive_search_dict_with_key_no_matches(self, aiida_profile):
        """Test when no Dict nodes contain the search key."""
        namespace = MagicMock(spec=ProcessBuilderNamespace)

        dict_without_key = orm.Dict(dict={'other': {'value': 100}})

        namespace._valid_fields = ['settings']
        namespace.get.side_effect = lambda key: {'settings': dict_without_key}[key]

        results = recursive_search_dict_with_key(namespace, 'incar')

        assert len(results) == 0


class TestRecursiveSearchPortBasename:
    """Test the recursive_search_port_basename function."""

    def test_recursive_search_port_basename_direct_match(self):
        """Test finding ports with matching basename."""
        namespace = MagicMock(spec=ProcessBuilderNamespace)

        test_value = 'test_calc'

        namespace._valid_fields = ['calc', 'settings', 'other']
        namespace.get.side_effect = lambda key: {
            'calc': test_value,
            'settings': 'other_value',
            'other': 'another_value',
        }[key]

        results = recursive_search_port_basename(namespace, 'calc')

        assert len(results) == 1
        assert results[0][0] == 'calc'
        assert results[0][1] == test_value

    def test_recursive_search_port_basename_nested(self):
        """Test recursive search for ports in nested namespaces."""
        sub_namespace = MagicMock(spec=ProcessBuilderNamespace)
        main_namespace = MagicMock(spec=ProcessBuilderNamespace)

        test_value = 'nested_calc'

        # Set up sub namespace
        sub_namespace._valid_fields = ['calc']
        sub_namespace.get.side_effect = lambda key: {'calc': test_value}[key]

        # Set up main namespace
        main_namespace._valid_fields = ['vasp']
        main_namespace.get.side_effect = lambda key: {'vasp': sub_namespace}[key]

        results = recursive_search_port_basename(main_namespace, 'calc')

        assert len(results) == 1
        assert results[0][0] == 'vasp.calc'
        assert results[0][1] == test_value

    def test_recursive_search_port_basename_multiple_matches(self):
        """Test finding multiple ports with same basename."""
        sub_namespace1 = MagicMock(spec=ProcessBuilderNamespace)
        sub_namespace2 = MagicMock(spec=ProcessBuilderNamespace)
        main_namespace = MagicMock(spec=ProcessBuilderNamespace)

        test_value1 = 'calc1'
        test_value2 = 'calc2'

        # Set up sub namespaces
        sub_namespace1._valid_fields = ['calc']
        sub_namespace1.get.side_effect = lambda key: {'calc': test_value1}[key]

        sub_namespace2._valid_fields = ['calc']
        sub_namespace2.get.side_effect = lambda key: {'calc': test_value2}[key]

        # Set up main namespace
        main_namespace._valid_fields = ['scf', 'bands']
        main_namespace.get.side_effect = lambda key: {'scf': sub_namespace1, 'bands': sub_namespace2}[key]

        results = recursive_search_port_basename(main_namespace, 'calc')

        assert len(results) == 2
        # Results could be in any order
        result_keys = [r[0] for r in results]
        result_values = [r[1] for r in results]
        assert 'scf.calc' in result_keys
        assert 'bands.calc' in result_keys
        assert test_value1 in result_values
        assert test_value2 in result_values

    def test_recursive_search_port_basename_no_matches(self):
        """Test when no ports match the basename."""
        namespace = MagicMock(spec=ProcessBuilderNamespace)

        namespace._valid_fields = ['settings', 'parameters']
        namespace.get.side_effect = lambda key: {'settings': 'some_value', 'parameters': 'other_value'}[key]

        results = recursive_search_port_basename(namespace, 'calc')

        assert len(results) == 0


class TestHasContent:
    """Test the has_content function."""

    def test_has_content_empty_dict(self):
        """Test has_content with empty dictionary."""
        mapping = {}
        assert has_content(mapping) is False

    def test_has_content_dict(self):
        """Test has_content with simple non-empty dictionary."""
        mapping = {'key': 'value'}
        assert has_content(mapping) is True

        # Test has_content with nested empty dictionaries.
        mapping = {'level1': {'level2': {}}}
        assert has_content(mapping) is False

        # Test has_content with nested dictionaries containing content.
        mapping = {'level1': {'level2': {'key': 'value'}}}
        assert has_content(mapping) is True

    def test_has_content_mixed_empty_and_content(self):
        """Test has_content with mix of empty and content dictionaries."""
        mapping = {'empty': {}, 'with_content': {'key': 'value'}}
        assert has_content(mapping) is True

    def test_has_content_deep_nested(self):
        """Test has_content with deeply nested structure."""
        mapping = {'level1': {'level2': {'level3': {'level4': {'key': 'value'}}}}}
        assert has_content(mapping) is True

    def test_has_content_all_nested_empty(self):
        """Test has_content with all nested dictionaries empty."""
        mapping = {'level1': {'level2a': {}, 'level2b': {'level3': {}}}, 'other': {}}
        assert has_content(mapping) is False

    def test_has_content_non_dict_values(self):
        """Test has_content with non-dictionary values."""
        mapping = {'string': 'value', 'number': 42, 'list': [1, 2, 3], 'none': None}
        # Should return True on first non-dict value found
        assert has_content(mapping) is True
