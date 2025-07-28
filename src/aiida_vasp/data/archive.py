"""
A general archive class.

Archive data class: store multiple files together in a compressed archive in the repository.
"""

from __future__ import annotations

import os
import tarfile
import tempfile
from contextlib import contextmanager
from typing import Any, Generator

from aiida.orm.nodes import Data


class ArchiveData(Data):
    """Compressed archive data node, contains a group of files that don't need to be readily accessible on their own."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self._filelist: list[tuple[str, str]] = []
        super().__init__(*args, **kwargs)

    @contextmanager
    def get_archive(self) -> Generator[tarfile.TarFile, None, None]:
        with self.base.repository.open('archive.tar.gz', mode='rb') as fobj:  # pylint: disable=not-context-manager
            with tarfile.open(fileobj=fobj, mode='r:gz') as tar:
                yield tar

    @contextmanager
    def archive(self) -> Generator[tarfile.TarFile, None, None]:
        with self.base.repository.open('archive.tar.gz', mode='rb') as fobj:  # pylint: disable=not-context-manager
            with tarfile.open(fileobj=fobj, mode='r:gz') as tar:
                yield tar

    def get_archive_list(self) -> None:
        with self.get_archive() as archive:
            return archive.list()

    def add_file(self, src_abs: str, dst_filename: str | None = None) -> None:
        if not dst_filename:
            dst_filename = os.path.basename(src_abs)
        self._filelist.append((src_abs, dst_filename))

    def _make_archive(self) -> None:
        """Create the archive file on disk with all it's contents."""
        _, path = tempfile.mkstemp()
        try:
            with tarfile.open(path, mode='w:gz') as archive:
                for src, dstn in self._filelist:
                    archive.add(src, arcname=dstn)
            self.base.repository.put_object_from_file(path, path='archive.tar.gz')
        finally:
            os.remove(path)

    # pylint: disable=arguments-differ, signature-differs
    def store(self, *args: Any, **kwargs: Any) -> Data:
        self._make_archive()
        del self._filelist
        return super().store(*args, **kwargs)
