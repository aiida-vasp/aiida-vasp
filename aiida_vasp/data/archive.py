# pylint: disable=abstract-method
# explanation: pylint wrongly complains about (aiida) Node not implementing query
"""Archive data class: store multiple files together in a compressed archive in the repository"""
import tarfile
import os
import six
from io import StringIO

from aiida.orm.nodes import Data


class ArchiveData(Data):
    """Compressed archive data node, contains a group of files that don't need to be readily accessible on their own"""

    def __init__(self, *args, **kwargs):
        self._filelist = []
        super(ArchiveData, self).__init__(*args, **kwargs)

    def get_archive(self):
        return tarfile.open(fileobj=self.open('archive.tar.gz', mode='rb'), mode='r:gz')

    def get_archive_list(self):
        archive = self.get_archive()
        archive.list()

    def add_file(self, src_abs, dst_filename=None):
        if not dst_filename:
            dst_filename = os.path.basename(src_abs)
        self._filelist.append((src_abs, dst_filename))

    def _make_archive(self):
        """Create the archive file on disk with all it's contents"""
        if six.PY2:
            self.put_object_from_filelike(StringIO.StringIO(), 'archive.tar.gz')
        else:
            self.put_object_from_filelike(StringIO(), 'archive.tar.gz')

        archive = tarfile.open(fileobj=self.open('archive.tar.gz', mode='wb'), mode='w:gz')

        for src, dstn in self._filelist:
            archive.add(src, arcname=dstn)

        archive.close()

    # pylint: disable=arguments-differ
    def store(self, *args, **kwargs):
        self._make_archive()
        del self._filelist
        super(ArchiveData, self).store(*args, **kwargs)

    @property
    def archive(self):
        return self.get_archive()
