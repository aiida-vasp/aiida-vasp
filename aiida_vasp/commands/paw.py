"""AiiDA Command plugin to work with .paw pseudopotential files and families"""
import sys

from aiida.cmdline.baseclass import VerdiCommandWithSubcommands
from aiida.cmdline.commands.data import Importable


def _try_load_dbenv():
    """
    Run `load_dbenv` unless the dbenv has already been loaded.
    """
    from aiida import load_dbenv, is_dbenv_loaded
    if not is_dbenv_loaded():
        load_dbenv()
        return True
    return False


class _Paw(VerdiCommandWithSubcommands, Importable):
    """
    Setup and manage paw pseudopotentials and families
    """

    def __init__(self):
        """setup and register subcommands"""
        from aiida_vasp.data.paw import PawData

        self.dataclass = PawData
        self.valid_subcommands = {
            'uploadfamily': (self.uploadfamily, self.complete_auto),
            'listfamilies': (self.listfamilies, self.complete_none),
            'import': (self.importfile, self.complete_none),
            'exportfamily': (self.exportfamily, self.complete_auto)
        }

    def uploadfamily(self, *args):
        """Upload a new PAW pseudopotential family."""
        import os.path
        import argparse as arp

        parser = arp.ArgumentParser(prog=self.get_full_command_name(), description='Upload a new PAW pseudopotential family.')
        parser.add_argument(
            '--stop-if-existing',
            action='store_true',
            dest='stop_if_existing',
            help='do not create the family or upload any '
            ' files, if existing files with matching '
            ' md5 sum are found.')
        parser.add_argument('folder')
        parser.add_argument('group_name')
        parser.add_argument('group_description')
        parser.set_defaults(stop_if_existing=False)

        params = parser.parse_args(args)
        folder = os.path.abspath(os.path.expanduser(params.folder))
        group_name = params.group_name
        group_description = params.group_description
        stop_if_existing = params.stop_if_existing

        if not os.path.isdir(folder):
            print >> sys.stderr, 'Cannot find directory: ' + folder
            sys.exit(1)

        _try_load_dbenv()
        from aiida.orm import DataFactory
        paw_cls = DataFactory('vasp.paw')
        files_found, files_uploaded = paw_cls.import_family(
            folder, familyname=group_name, family_desc=group_description, stop_if_existing=stop_if_existing)

        print "POTCAR files found in subfolders: {}. New files uploaded from: {}".format(files_found, files_uploaded)

    def listfamilies(self, *args):
        """Subcommand to list AiiDA PAW families present in the DB"""
        import argparse

        parser = argparse.ArgumentParser(prog=self.get_full_command_name(), description='List AiiDA PAW families.')
        parser.add_argument(
            '-e',
            '--element',
            nargs='+',
            type=str,
            default=[],
            help="Filter the families only to those containing "
            "potentials for each of the specified elements")
        parser.add_argument(
            '-s',
            '--symbol',
            nargs='+',
            type=str,
            default=[],
            help="Filter the families only to those containing "
            "a potential for each of the specified symbols")
        parser.add_argument(
            '-d', '--with-description', dest='with_description', action='store_true', help="Show also the description for the PAW family")
        parser.set_defaults(with_description=False)

        params = parser.parse_args(args)

        _try_load_dbenv()
        from aiida.orm import DataFactory
        paw_cls = DataFactory('vasp.paw')
        groups = paw_cls.get_paw_groups(elements=list(params.element), symbols=list(params.symbol))

        if groups:
            for group in groups:
                paws = paw_cls.query(dbgroups=group.dbgroup).distinct()
                num_paws = paws.count()
                description_string = ''
                if params.with_description:
                    description_string = ': {}'.format(group.description)
                groupitem = '* {} [{} pseudos]{}'
                print groupitem.format(group.name, num_paws, description_string)
        else:
            print 'No PAW pseudopotential family found.'

    @staticmethod
    def _import_paw_parameters(parser):
        parser.add_argument(
            '--psctr',
            nargs=1,
            type=str,
            default=None,
            help='path to the psctr file, only necessary if you wish '
            'to store it and you are not passing a folder which contains it.')

    def _import_paw(self, filename, **kwargs):
        """Import a PAW potential."""
        from os.path import expanduser

        # ~ parser.add_argument('path', help='path to a file or a folder. '
        # ~ 'file: must be in POTCAR format'
        # ~ 'folder: must contain one POTCAR and optionally a PSCTR file')
        path = expanduser(filename)
        psctr = kwargs.get('psctr')
        psctr = expanduser(psctr[0]) if psctr else None
        try:
            node, created = self.dataclass.get_or_create(path, psctr=psctr, store=True)
            print repr(node)
            if not created:
                print 'part of the following families:'
                for group in node.dbnode.dbgroups.all():
                    print ' * {}'.format(group.name)
        except ValueError as err:
            print err

    # pylint: disable=too-many-locals
    def exportfamily(self, *args):
        """Export a PAW potential family into a folder"""
        import argparse
        import os
        from os.path import abspath, expanduser

        from aiida.orm import DataFactory
        from aiida.common.exceptions import NotExistent

        parser = argparse.ArgumentParser(prog=self.get_full_command_name(), description='Export a PAW potential family into a folder')
        parser.add_argument('--name', nargs=1, type=str, default=None, help='name of the folder, defaults to ' 'potpaw_<family_name>.')
        parser.add_argument('folder', help='path to where the potpaw ' 'folder should be created')
        parser.add_argument('family_name')

        params = parser.parse_args(args)
        folder = abspath(expanduser(params.folder))

        _try_load_dbenv()
        paw_cls = DataFactory('vasp.paw')
        try:
            group = paw_cls.get_famgroup(params.family_name)
        except NotExistent:
            print >> sys.stderr, ('paw family {} not found'.format(params.family_name))

        #create folder
        potpaw = params.name[0] if params.name else os.path.join(folder, 'potpaw_%s' % params.family_name)
        for paw in group.nodes:
            pawf = os.path.join(potpaw, paw.symbol)
            if not os.path.exists(pawf):
                os.makedirs(pawf)
            potcar_path = os.path.join(pawf, 'POTCAR')
            psctr_path = os.path.join(pawf, 'PSCTR')
            isfile_msg = 'file {} is already present in the destination folder.'
            if not os.path.isfile(potcar_path):
                with open(potcar_path, 'w') as dest:
                    with open(paw.potcar) as src:
                        dest.write(src.read())
            else:
                print isfile_msg.format(os.path.join(paw.symbol, 'POTCAR'))
            if not os.path.isfile(psctr_path):
                try:
                    with open(paw.psctr) as src:
                        with open(psctr_path, 'w') as dest:
                            dest.write(src.read())
                except OSError:
                    pass
            else:
                print isfile_msg.format(os.path.join(paw.symbol, 'PSCTR'))
