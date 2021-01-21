"""
A general archive class.

------------------------
Archive data class: store multiple files together in a compressed archive in the repository.
"""
# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
import tarfile
from contextlib import contextmanager
import os
import tempfile
from aiida.orm.nodes import Data


class ArchiveData(Data):
    """Compressed archive data node, contains a group of files that don't need to be readily accessible on their own."""

    def __init__(self, *args, **kwargs):
        self._filelist = []
        super(ArchiveData, self).__init__(*args, **kwargs)

    @contextmanager
    def get_archive(self):
        with self.open('archive.tar.gz', mode='rb') as fobj:
            with tarfile.open(fileobj=fobj, mode='r:gz') as tar:
                yield tar

    @contextmanager
    def archive(self):
        with self.open('archive.tar.gz', mode='rb') as fobj:
            with tarfile.open(fileobj=fobj, mode='r:gz') as tar:
                yield tar

    def get_archive_list(self):
        with self.get_archive() as archive:
            return archive.list()

    def add_file(self, src_abs, dst_filename=None):
        if not dst_filename:
            dst_filename = os.path.basename(src_abs)
        self._filelist.append((src_abs, dst_filename))

    def _make_archive(self):
        """Create the archive file on disk with all it's contents."""
        _, path = tempfile.mkstemp()
        try:
            with tarfile.open(path, mode='w:gz') as archive:
                for src, dstn in self._filelist:
                    archive.add(src, arcname=dstn)
            self.put_object_from_file(path, path='archive.tar.gz')
        finally:
            os.remove(path)

    # pylint: disable=arguments-differ, signature-differs
    def store(self, *args, **kwargs):
        self._make_archive()
        del self._filelist
        super(ArchiveData, self).store(*args, **kwargs)
