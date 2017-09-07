"""
Utils to query for Vasp calculations

THIS IS LIKELY OUTDATED - It was written before QueryBuilder was added to aiida-core
"""
from aiida.orm.querytool import QueryTool
from aiida.orm.calculation.job import JobCalculation


class VaspFinder(object):
    """
    Utility class to query for Vasp calculations

    THIS IS LIKELY OUTDATED - It was written before QueryBuilder was added to aiida-core
    """

    _vaspclass = [
        'vasp.vasp', 'vasp.asevasp', 'vasp.vasp5', 'vasp.vasp2w90',
        'vasp.vasp2w90.Vasp2W90Calculation'
    ]

    @classmethod
    def cmp_ctime(cls, calc1, calc2):
        """Compare calculations by creation time"""
        time_1 = calc1.ctime
        time_2 = calc2.ctime
        if time_1 < time_2:
            return -1
        elif time_1 > time_2:
            return 1
        return 0

    @classmethod
    def table_element(cls, calc):
        """Create a table element from a calculation"""
        ex = calc.get_extras()
        success = ex.get('success')
        if isinstance(success, (bool, int)):
            success_str = 'yes' if success else 'no'
        else:
            success_str = success
        element = {
            'creation_time': calc.ctime.strftime(format='%Y-%m-%d %H:%M'),
            'ctime': calc.ctime,
            'class': calc.__class__.__name__,
            'successful': success_str and success_str or 'N/A',
            'experiment': ex.get('experiment', 'N/A'),
            'exp_run': ex.get('run', 'N/A'),
            'type': ex.get('type', 'N/A'),
            'state': calc.get_state(),
            'pk': calc.pk
        }
        return element

    @classmethod
    def str_table(cls, tab):
        """Convert a table of calculations to a string"""
        header = {
            'creation_time': 'Creation Time',
            'class': 'Class',
            'successful': 'Success',
            'experiment': 'Experiment',
            'exp_run': 'Run Nr.',
            'type': 'Tags',
            'state': 'AiiDA-State',
            'pk': 'PK'
        }
        line = '{pk:>5} {creation_time:18} {state:20} {successful:>6} '
        line += '{experiment:20} {exp_run:>7} {class} {type}'
        tab.insert(0, header)
        return '\n'.join([line.format(**c) for c in tab])

    @classmethod
    def cstate(cls, calc):
        if hasattr(calc, 'get_state'):
            return calc.get_state()
        return 'N/A'

    @classmethod
    def history(cls, last=0):
        """
        Print a history of the most recent calculations

        :param last: int, show this many calculations
        """
        query_tool = QueryTool()
        query_tool.set_class(JobCalculation)
        calc_list = [
            calc for calc in query_tool.run_query()
            if 'vasp' in str(calc.__class__)
        ]
        calc_list.sort(cls.cmp_ctime)
        res = calc_list[-last:]
        res.reverse()
        tab = [cls.table_element(c) for c in res]
        print cls.str_table(tab)

    @classmethod
    def status(cls):
        """Print a table of calculations with the state they are in"""
        query_tool = QueryTool()
        query_tool.set_class(JobCalculation)
        states = [
            'TOSOBMIT', 'SUBMITTING', 'WITHSCHEDULER', 'COMPUTED', 'PARSING',
            'RETRIEVING'
        ]
        res_list = [
            calc for calc in query_tool.run_query()
            if cls.cstate(calc) in states
        ]
        res_list.sort(cls.cmp_ctime)
        res_list.reverse()
        tab = [cls.table_element(c) for c in res_list]
        print cls.str_table(tab)
