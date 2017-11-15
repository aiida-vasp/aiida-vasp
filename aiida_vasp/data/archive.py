# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""Archive data class: store multiple files together in a compressed archive in the repository"""
import tarfile
import os
import StringIO

from aiida.orm.data import Data


class ArchiveData(Data):
    """Compressed archive data node, contains a group of files that don't need to be readily accessible on their own"""

    def __init__(self, *args, **kwargs):
        super(ArchiveData, self).__init__(*args, **kwargs)
        self._filelist = []

    def get_archive_abs_path(self):
        return self.get_abs_path('archive.tar.gz')

    def get_archive(self):
        return tarfile.open(self.get_archive_abs_path(), mode='r')

    def get_archive_list(self):
        archive = self.get_archive()
        archive.list()

    def add_file(self, src_abs, dst_filename=None):
        if not dst_filename:
            dst_filename = os.path.basename(src_abs)
        self._filelist.append((src_abs, dst_filename))

    def _make_archive(self):
        """Create the archive file on disk with all it's contents"""
        self.folder.create_file_from_filelike(StringIO.StringIO(), 'path/archive.tar.gz')
        archive = tarfile.open(self.get_archive_abs_path(), mode='w:gz')
        for src, dstn in self._filelist:
            archive.add(src, arcname=dstn)
        archive.close()

    def store(self, with_transaction=True):
        self._make_archive()
        del self._filelist
        super(ArchiveData, self).store(with_transaction)

    @property
    def archive(self):
        return self.get_archive()
