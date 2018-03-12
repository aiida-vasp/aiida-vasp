"""Run a vasp workflow from json input."""

import argparse
from aiida import load_dbenv
load_dbenv()


def get_parser():
    """Create a cmdline parser for the tool"""
    parser = argparse.ArgumentParser(description=('Run a vasp workflow reading parameters ' 'from a json file'))
    parser.add_argument(
        '--store-template', action='store_true', help=('store an input template '
                                                       'in <input file> instead of running the workflow'))
    # ~ parser.add_argument('-g', '--group', required=False,
    # ~ help=('aiida group to store the wf in'))
    parser.add_argument('workflow', help=('a valid string to load ' 'a workflow using WorkflowFactory'))
    parser.add_argument('input_file', help=('a .json file containing ' 'all necessary parameters for the workflow'))

    return parser


def main():
    """Write or read a workflow input file and start a workflow if requested"""
    parser = get_parser()
    args = parser.parse_args()
    from aiida.orm import WorkflowFactory
    from os.path import expanduser, abspath
    import json
    workflow_cls = WorkflowFactory(args.workflow)
    if args.store_template:
        workflow_cls().get_template(path=args.input_file)
    else:
        path = abspath(expanduser(args.input_file))
        with open(path) as inputf:
            params = json.load(inputf)
            workflow = workflow_cls(params=params)
            workflow.label = params.get('label')
            valid, log = workflow.helper._verify_params(  # pylint: disable=protected-access
                workflow.get_parameters(), silent=True)
            if not valid:
                raise IOError('invalid input:\n' + log)
            # ~ wf.store()
            workflow.start()
            print '\n'.join(workflow.get_report())


if __name__ == '__main__':
    main()
