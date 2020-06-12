"""Utils to work with Wannier90 .win format."""
import re

from aiida_vasp.parsers.file_parsers.parser import KeyValueParser


class WinParser(KeyValueParser):
    """Parses wannier90.win files."""
    block = re.compile(r'begin (?P<name>\w*)\s*\n\s*(?P<content>[\w\W]*)\s*\n\s*end \1')
    comment = re.compile(r'(!.*)\n?')

    def __init__(self, file_path):
        super(WinParser, self).__init__()
        self.result = {}
        with open(file_path) as winf:
            self.keywords, self.blocks, self.comments = WinParser.parse_win(winf)
        self.result.update(self.keywords)
        self.result.update(self.blocks)

    @classmethod
    def parse_win(cls, fobj_or_str):
        """Parse a wannier90 .win file."""
        if isinstance(fobj_or_str, str):
            content = fobj_or_str
        else:
            content = fobj_or_str.read()
        comments = re.findall(cls.comment, content)
        content = re.sub(cls.comment, '', content)
        blocks = re.findall(cls.block, content)
        content = re.sub(cls.block, '', content)
        kvd = dict(re.findall(cls.assignment, content))
        bld = {}
        for keyword, value in blocks:
            # do not split individual lines
            bld[keyword] = [line.strip() for line in value.split('\n')]
        return kvd, bld, comments
