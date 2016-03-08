# ~ import tempfile
import tarfile
from aiida.orm.data import Data
import os
import StringIO


class ArchiveData(Data):
    def __init__(self, *args, **kwargs):
        super(ArchiveData, self).__init__(*args, **kwargs)
        self._filelist = []

    def get_archive_abs_path(self):
        return self.get_abs_path('archive.tar.gz')

    def get_archive(self):
        return tarfile.open(self.get_archive_abs_path(), mode='r')

    def get_archive_list(self):
        ar = get_archive()
        ar.list()

    def add_file(self, src_abs, dst_filename=None):
        if not dst_filename:
            dst_filename = os.path.basename(src_abs)
        self._filelist.append((src_abs, dst_filename))

    def _make_archive(self):
        self.folder.create_file_from_filelike(
            StringIO.StringIO(), 'path/archive.tar.gz')
        ar = tarfile.open(self.get_archive_abs_path(), mode='w:gz')
        for src, dstn in self._filelist:
            ar.add(src, arcname=dstn)
        ar.close()

    def store(self, *args, **kwargs):
        self._make_archive()
        del(self._filelist)
        super(ArchiveData, self).store(*args, **kwargs)

    @property
    def archive(self):
        return self.get_archive()
