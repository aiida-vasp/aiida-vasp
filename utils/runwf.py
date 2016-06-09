from aiida import load_dbenv
load_dbenv()

import argparse

parser = argparse.ArgumentParser(
    description=('Run a vasp workflow reading parameters '
                 'from a json file'))
parser.add_argument('--store-template', action='store_true', help=('store an input template '
                    'in <input file> instead of running the workflow'))
parser.add_argument('workflow', help=('a valid string to load '
                    'a workflow using WorkflowFactory'))
parser.add_argument('input_file', help=('a .json file containing '
                    'all necessary parameters for the workflow'))

if __name__ == '__main__':
    args = parser.parse_args()
    from aiida.orm import WorkflowFactory
    from os.path import expanduser, abspath
    import json
    WorkflowClass = WorkflowFactory(args.workflow)
    if args.store_template:
        WorkflowClass().get_template(path=args.input_file)
    else:
        path = abspath(expanduser(args.input_file))
        with open(path) as inputf:
            params = json.load(inputf)
            wf = WorkflowClass(params=params)
            wf.label = params.get('label')
            valid, log = wf.helper._verify_params(wf.get_parameters(), silent=True)
            if not valid:
                raise IOError('invalid input:\n' + log)
            # ~ wf.store()
            wf.start()
