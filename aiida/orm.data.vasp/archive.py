import tempfile
import tarfile
from aiida.orm.data.singlefile import SinglefileData
import os
import shutil


class ArchiveData(SinglefileData):
    def __init__(self):
        super(ArchiveData, self).__init__()
        with tempfile.NamedTemporaryFile() as tmp:
            self.add_path(tmp.name, dst_filename='archive.tar.gz')
        self._filelist = []

    def add_file(self, src_abs, dst_filename=None):
        if not dst_filename:
            dst_filename = os.path.basename(src_abs)
        self._filelist.append((src_abs, dst_filename))

    def _make_archive(self):
        # ~ with tarfile.open(fileobj=tmp, mode='w:gz') as ar:
            # ~ for src, dstn in self._filelist:
                # ~ ar.add(src, arcname=dstn)
        # ~ mytemp = os.path.expanduser('~/tmp/archive.tar.gz')
        # ~ shutil.copy(tmp.name, mytemp)
        # ~ self.add_path(mytemp)
        ar = tarfile.open(self.get_file_abs_path(), mode='w:gz')
        for src, dstn in self._filelist:
            ar.add(src, arcname=dstn)
        ar.close()

    def store(self):
        self._make_archive()
        del(self._filelist)
        super(ArchiveData, self).store()
