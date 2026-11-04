"""Fortran list-directed text reading (``read(unit, *)``): numbers with d exponents, quoted strings."""


def fortran_float(s):
    """Float from a Fortran number string (d/D exponent allowed)."""
    return float(s.replace('d', 'e').replace('D', 'e'))


def parse_floats(s):
    """Whitespace-separated floats, Fortran d exponent allowed."""
    return [fortran_float(x) for x in s.split()]


def parse_ints(s):
    """Whitespace-separated ints, read through fortran_float."""
    return [int(fortran_float(x)) for x in s.split()]


class FortranFile:
    """Fortran list-directed reading (read(unit,*)).

    - skips blank lines
    - reads N tokens across several lines
    - handles single- or double-quoted strings
    """

    def __init__(self, f):
        self.f = f
        self._buffer = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    @staticmethod
    def _tokenize(line):
        """Split a line on blanks and commas, keeping quoted strings whole."""
        tokens = []
        i = 0
        while i < len(line):
            while i < len(line) and line[i] in ' \t,':
                i += 1
            if i >= len(line):
                break
            if line[i] in "'\"":
                quote = line[i]
                j = i + 1
                while j < len(line) and line[j] != quote:
                    j += 1
                tokens.append(line[i + 1 : j])
                i = j + 1
            else:
                j = i
                while j < len(line) and line[j] not in ' \t,':
                    j += 1
                tokens.append(line[i:j])
                i = j
        return tokens

    def readline(self):
        """Next non-blank line, '' at end of file."""
        while True:
            line = self.f.readline()
            if not line:
                return ''
            if line.strip():
                return line

    def read_tokens(self, n):
        """Exactly n tokens, read across lines if needed; fewer only at end of file."""
        result = list(self._buffer)
        self._buffer = []
        while len(result) < n:
            line = self.readline().strip()
            if not line:
                break
            result.extend(self._tokenize(line))
        if len(result) > n:
            self._buffer = result[n:]
            result = result[:n]
        return result
