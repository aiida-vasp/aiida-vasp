import pytest

from aiida_vasp.utils.dict_merge import recursive_merge


@pytest.mark.parametrize(
    'left,right,expected',
    [
        # Simple merge
        ({'a': 1, 'b': {'c': 2}}, {'b': {'c': 3}, 'd': 4}, {'a': 1, 'b': {'c': 3}, 'd': 4}),
        # Nested merge and list extend
        (
            {'a': {'x': 1, 'y': {'z': 2}}, 'b': [1, 2]},
            {'a': {'y': {'z': 3}}, 'b': {'$!extend': [3, 4]}},
            {'a': {'x': 1, 'y': {'z': 3}}, 'b': [1, 2, 3, 4]},
        ),
        # Special: append and delete
        ({'a': [1, 2], 'b': {'c': 5}}, {'a': {'$!append': 3}, 'b': {'$!del': True}}, {'a': [1, 2, 3]}),
        ({'a': [1, 2], 'b': {'c': 5}}, {'a': {'$!append': 3}, 'b': '$!del'}, {'a': [1, 2, 3]}),
        # Replace
        ({'a': {'x': 1}}, {'a': {'$!replace': {'y': 2}}}, {'a': {'y': 2}}),
    ],
)
def test_recursive_merge(left, right, expected):
    assert recursive_merge(left, right) == expected
