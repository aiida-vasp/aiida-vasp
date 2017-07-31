from aiida.orm.querytool import QueryTool
from aiida.orm.calculation.job import JobCalculation


class VaspFinder(object):
    _vaspclass = ['vasp.vasp', 'vasp.vasp2w90']

    @classmethod
    def cmp_ctime(cls, calc1, calc2):
        t1 = calc1.ctime
        t2 = calc2.ctime
        if t1 < t2:
            return -1
        elif t1 > t2:
            return 1
        else:
            return 0

    @classmethod
    def table_element(cls, calc):
        ex = calc.get_extras()
        s = ex.get('success')
        if isinstance(s, (bool, int)):
            ok = s and 'yes' or 'no'
        else:
            ok = s
        element = {
            'creation_time': calc.ctime.strftime(format='%Y-%m-%d %H:%M'),
            'ctime': calc.ctime,
            'class': calc.__class__.__name__,
            'successful': ok and ok or 'N/A',
            'experiment': ex.get('experiment', 'N/A'),
            'exp_run': ex.get('run', 'N/A'),
            'type': ex.get('type', 'N/A'),
            'state': calc.get_state(),
            'pk': calc.pk
        }
        return element

    @classmethod
    def str_table(cls, tab):
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
        else:
            return 'N/A'

    @classmethod
    def history(cls, last=0, vaspclass=None):
        q = QueryTool()
        q.set_class(JobCalculation)
        l = filter(lambda c: 'vasp' in str(c.__class__), q.run_query())
        l.sort(cls.cmp_ctime)
        res = l[-last:]
        res.reverse()
        tab = [cls.table_element(c) for c in res]
        print cls.str_table(tab)

    @classmethod
    def status(cls, vaspclass=None):
        q = QueryTool()
        q.set_class(JobCalculation)
        st = ['TOSOBMIT', 'SUBMITTING', 'WITHSCHEDULER',
              'COMPUTED', 'PARSING', 'RETRIEVING']
        l = filter(lambda c: cls.cstate(c) in st, q.run_query())
        l.sort(cls.cmp_ctime)
        l.reverse()
        tab = [cls.table_element(c) for c in l]
        print cls.str_table(tab)
