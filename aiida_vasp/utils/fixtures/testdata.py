"""
Test data retrieval.

--------------------
Here we provide pytest fixtures to make it more convenient to retrieving test data.
"""
import os


def data_path(*args):
    """Give the path to test data."""
    path = os.path.realpath(os.path.join(__file__, '../../../test_data', *args))
    assert os.path.exists(path)
    assert os.path.isabs(path)
    return path


def read_file(*args, **kwargs):
    """Give the content (string) of a test data file."""
    path = kwargs.get('path', None)
    mode = kwargs.pop('mode', None)
    encoding = kwargs.pop('encoding', None)
    if not mode:
        mode = 'r'
    if not path:
        path = data_path(*args)
    if encoding is not None and 'b' not in mode:
        encoding = 'utf8'
    with open(path, mode, encoding=encoding) as testdata_fo:
        testdata_content = testdata_fo.read()
    return testdata_content
