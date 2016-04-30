from aiida.cmdline.baseclass import VerdiCommandWithSubcommands
from aiida.cmdline.commands.data import Importable
import sys


class _Paw(VerdiCommandWithSubcommands, Importable):
    '''
    Setup and manage paw pseudopotentials and families
    '''

    def __init__(self):
        '''setup and register subcommands'''
        from aiida.orm.data.vasp.paw import PawData

        self.dataclass = PawData
        self.valid_subcommands = {
            'uploadfamily': (self.uploadfamily, self.complete_auto),
            'listfamilies': (self.listfamilies, self.complete_none),
            'import': (self.importfile, self.complete_none),
            'exportfamily': (self.exportfamily, self.complete_auto)
        }

    def uploadfamily(self, *args):
        '''Upload a new PAW pseudopotential family.'''
        from aiida import load_dbenv
        import os.path
        import argparse as arp

        parser = arp.ArgumentParser(
            prog=self.get_full_command_name(),
            description='Upload a new PAW pseudopotential family.')
        parser.add_argument('--stop-if-existing', action='store_true',
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
            print >> sys.stderrm, 'Cannot find directory: ' + folder
            sys.exit(1)

        load_dbenv()
        from aiida.orm import DataFactory
        Paw = DataFactory('vasp.paw')
        files_found, files_uploaded = Paw.import_family(folder,
                                                        familyname=group_name,
                                                        family_desc=group_description,
                                                        stop_if_existing=stop_if_existing)

        print "POTCAR files found in subfolders: {}. New files uploaded from: {}".format(files_found, files_uploaded)

    def listfamilies(self, *args):
        from aiida import load_dbenv
        import argparse

        parser = argparse.ArgumentParser(
            prog=self.get_full_command_name(),
            description='List AiiDA PAW families.')
        parser.add_argument('-e', '--element', nargs='+', type=str, default=[],
                            help="Filter the families only to those containing "
                                 "potentials for each of the specified elements")
        parser.add_argument('-s', '--symbol', nargs='+', type=str, default=[],
                            help="Filter the families only to those containing "
                                 "a potential for each of the specified symbols")
        parser.add_argument('-d', '--with-description',
                            dest='with_description', action='store_true',
                            help="Show also the description for the PAW family")
        parser.set_defaults(with_description=False)

        params = parser.parse_args(args)

        load_dbenv()
        from aiida.orm import DataFactory
        Paw = DataFactory('vasp.paw')
        groups = Paw.get_paw_groups(elements=list(params.element),
                                    symbols=list(params.symbol))

        if groups:
            for g in groups:
                paws = Paw.query(dbgroups=g.dbgroup).distinct()
                num_paws = paws.count()
                description_string = ''
                if params.with_description:
                    description_string = ': {}'.format(g.description)
                groupitem = '* {} [{} pseudos]{}'
                print groupitem.format(g.name, num_paws, description_string)
        else:
            print 'No PAW pseudopotential family found.'

    def _import_paw_parameters(self, parser):
        parser.add_argument('--psctr', nargs=1, type=str, default=None,
                            help='path to the psctr file, only necessary if you wish '
                            'to store it and you are not passing a folder which contains it.')

    def _import_paw(self, filename, **kwargs):
        '''Import a PAW potential.'''
        from os.path import expanduser

        # ~ parser.add_argument('path', help='path to a file or a folder. '
                            # ~ 'file: must be in POTCAR format'
                            # ~ 'folder: must contain one POTCAR and optionally a PSCTR file')
        path = expanduser(filename)
        psctr = kwargs.get('psctr')
        psctr = psctr and expanduser(psctr[0]) or None
        try:
            node, created = self.dataclass.get_or_create(path, psctr=psctr, store=True)
            print repr(node)
            if not created:
                print 'part of the following families:'
                for g in node.dbnode.dbgroups.all():
                    print ' * {}'.format(g.name)
        except ValueError as e:
            print e

    def exportfamily(self, *args):
        '''Export a PAW potential family into a folder'''
        from aiida import load_dbenv
        from aiida.orm import DataFactory
        from aiida.common.exceptions import NotExistent
        import argparse
        import os
        from os.path import abspath, expanduser

        parser = argparse.ArgumentParser(
            prog=self.get_full_command_name(),
            description='Export a PAW potential family into a folder')
        parser.add_argument('--name', nargs=1, type=str, default=None,
                            help='name of the folder, defaults to '
                            'potpaw_<family_name>.')
        parser.add_argument('folder', help='path to where the potpaw '
                            'folder should be created')
        parser.add_argument('family_name')

        params = parser.parse_args(args)
        folder = abspath(expanduser(params.folder))

        load_dbenv()
        Paw = DataFactory('vasp.paw')
        try:
            group = Paw.get_famgroup(params.family_name)
        except NotExistent:
            print >> sys.stderr, ('paw family {} not found'.format(params.family_name))

        #create folder
        potpaw = params.name and params.name[0] or os.path.join(folder, 'potpaw_%s' % params.family_name)
        for paw in group.nodes:
            pawf = os.path.join(potpaw, paw.symbol)
            if not os.path.exists(pawf):
                os.makedirs(pawf)
            pp = os.path.join(pawf, 'POTCAR')
            cp = os.path.join(pawf, 'PSCTR')
            isfile_msg = 'file {} is already present in the destination folder.'
            if not os.path.isfile(pp):
                with open(pp, 'w') as dest:
                    with open(paw.potcar) as src:
                        dest.write(src.read())
            else:
                print isfile_msg.format(os.path.join(paw.symbol, 'POTCAR'))
            if not os.path.isfile(cp):
                try:
                    with open(paw.psctr) as src:
                        with open(cp, 'w') as dest:
                            dest.write(src.read())
                except OSError:
                    pass
            else:
                print isfile_msg.format(os.path.join(paw.symbol, 'PSCTR'))

