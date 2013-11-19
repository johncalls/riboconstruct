import math
import os.path
import shutil
import subprocess
import sys
import tempfile
import time


CMP_THR = 0.00000001
MAX_FLOAT = 10000000.0
MIN_FLOAT = -MAX_FLOAT


# RNAf
BB_FILE_EXT = "bb"
BB_TERMINATOR = "<<<"
PATTERN_FILE_NAME = "pattern.p"


def enum(*sequential, **named):
    """
    Helper function to construct C#-like enums.

    Example: ::

        >>> Numbers = enum('ZERO', 'One', five=5)
        >>> Numbers.ZERO
        0
        >>> Numbers.five
        5
        >>> Numbers.ident[1]  # reverse mapping
        'One'
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    reverse = dict((value, key) for key, value in enums.iteritems())
    enums['ident'] = reverse
    return type('Enum', (), enums)


def max_sum(*l):
    """
    Sums up the committed parameters. Returns *MAX_FLOAT* if one of the
    values exceeds this value.
    """
    if max(l) > MAX_FLOAT:
        return MAX_FLOAT
    try:
        return sum(l)
    except TypeError:
        return sum(filter(None, l))


def min_sum(a, b):
    """
    Subtracts *b* from *a*. Returns *MIN_FLOAT* if one of the values
    exceeds *MAX_FLOAT*.
    """
    if a > MAX_FLOAT or b > MAX_FLOAT:
        return MIN_FLOAT
    return a - b


def float_equal(a, b):
    """
    Returns whether floats *a* and *b* are equal up to some *threshold*.
    """
    try:
        return math.fabs(a - b) < CMP_THR
    except TypeError:
        return False


class RedirectStdStreams(object):
    """
    Can be used in combination with the ``with`` statement to redirect
    output streams.
    """

    def __init__(self, stdout=sys.stdout, stderr=sys.stderr):
        self._stdout = stdout
        self._stderr = stderr

    def __enter__(self):
        self.old_stdout = sys.stdout
        self.old_stderr = sys.stderr
        self.old_stdout.flush()
        self.old_stderr.flush()
        sys.stdout = self._stdout
        sys.stderr = self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush()
        self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr


def redirect_std_streams(stderr=sys.stderr, stdout=sys.stdout):
    """Decorator function to redirect output streams."""
    def wrap(f):
        def newf(*args, **kwargs):
            old_stderr, old_stdout = sys.stderr, sys.stdout
            sys.stderr = stderr
            sys.stdout = stdout
            try:
                return f(*args, **kwargs)
            finally:
                sys.stderr, sys.stdout = old_stderr, old_stdout
        return newf
    return wrap


def boltzmann(s, T=None, tries=0):
    """
    Computes basepair probabilities and returns as a list of
    dictionaries with the None entry giving the unpaired probability.
    `RNAfold <http://www.tbi.univie.ac.at/~ronny/RNA/RNAfold.html>`_
    is used to calculate the pairing probabilities and has to be
    installed on the machine :func:`boltzmann` is called from.
    """
    # Class for maintaining probabilities, if none is present return 0
    class _ProbVec(dict):
        def __getitem__(self, key):
            if key in self:
                return super(_ProbVec, self).__getitem__(key)
            else:
                return 0.0

    # Function defining behaviour if a folding program fails
    def onfail(f, tries, *args):
        if tries < 10:
            # Don't retry more than 10 times
            time.sleep(10)  # Take ten
            return f(*args, tries=tries + 1)
        else:
            raise RuntimeError("Function " + f.__name__ +
                               " did not receive result from external program")

    cmd = ["RNAfold", "-p"]
    tmpdir = tempfile.mkdtemp()
    if T is not None:
        cmd.extend(["-T", str(T)])
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         cwd=tmpdir)
    print >> p.stdin, s
    p.stdin.close()
    t = p.stdout.readlines()
    while p.stdout.readline() != "":
        # Make sure RNAfold has properly finished before reading probability plot
        pass
    p.stdout.close()
    try:
        mfe = t[1].strip().split(None, 1)
        ensemble = float(t[2].strip().split(None, 1)[1][1:-1])
        if mfe == [] or len(mfe[0]) != len(s):
            # Did not receive expected output from RNAfold
            raise IndexError
    except IndexError:
        return onfail(boltzmann, tries, s, T)
    p = [_ProbVec() for i in xrange(len(s))]
    for i in p:
        i[None] = 1.0
    with open(os.path.join(tmpdir, "dot.ps")) as fh:
        t = fh.readline()
        while "data starts here" not in t:
            t = fh.readline()
        t = fh.readline()
        while t != "" and "showpage" not in t:
            if "ubox" in t:
                t = t.split()
                i = int(t[0]) - 1
                j = int(t[1]) - 1
                q = pow(float(t[2]), 2)
                p[i][j] = q
                p[j][i] = q
                p[i][None] -= q
                p[j][None] -= q
            t = fh.readline()
        shutil.rmtree(tmpdir, ignore_errors=True)
        return (p, mfe[0], float(mfe[1][1:-1]), ensemble)


def fold(seq, struct=None):
    """
    Folds *seq* using
    `RNAfold <http://www.tbi.univie.ac.at/~ronny/RNA/RNAfold.html>`_
    which has to be installed on the machine :func:`fold` is called from.
    Returns the structure and the minimum free energy with which *seq*
    folds into this structure.
    """
    cmd = ["RNAfold", "--noPS"]
    tmp_dir = tempfile.mkdtemp(dir=os.getcwd())
    if struct is not None:
        cmd.extend(["-C %s" % struct])
    p = subprocess.Popen(cmd,
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         cwd=tmp_dir)
    print >> p.stdin, str(seq)
    p.stdin.close()
    t = p.stdout.readlines()[-1].strip().split(None, 1)
    p.stdout.close()
    shutil.rmtree(tmp_dir, ignore_errors=True)
    return t[0], float(t[1][1:-1])


def rnaplfold(seq, positions, lengths=None, L=70):
    """
    Folds *seq* using
    `RNAplfold <http://www.tbi.univie.ac.at/~ronny/RNA/RNAplfold.html>`_
    which has to be installed on the machine
    :func:`rnaplfold` is called from and returns
    the local accessibilities for each stretch starting at a position in
    *positions*. The length of each stretch is defined in *lengths*. If
    *lengths* is ``None`` than each stretch is of length 5.
    """
    ordered_positions = sorted(((p, l), i) for i, (p, l) in
                               enumerate(zip(positions, lengths)))
    if lengths is None:
        u = 5
    else:
        u = max(lengths)

    cmd = ["RNAplfold", "-u %i" % u, "-W %i" % L]
    tmp_dir = tempfile.mkdtemp(dir=os.getcwd())
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, cwd=tmp_dir)
    print >> p.stdin, str(seq)
    p.stdin.close()
    p.wait()

    pos_iter = iter(ordered_positions)
    accessibilities = dict()
    (pos, u), j = pos_iter.next()
    pos += u
    with open(os.path.join(tmp_dir, "plfold_lunp")) as fh:
        # skip the first two lines
        fh.next()
        fh.next()
        for i, line in enumerate(fh, 1):
            # wait for the line that has been asked for
            if i == pos:
                if i + 1 < u:
                    accessibilities[j] = float(line.rstrip().split('\t')[i + 1])
                else:
                    accessibilities[j] = float(line.rstrip().split('\t')[u])
                try:
                    (pos, u), j = pos_iter.next()
                    pos += u
                except StopIteration:
                    break
    shutil.rmtree(tmp_dir)

    return accessibilities.values()
