import copy
import sys

from . import element as rs_e
from .. import rna


def iter_riboswitch(config):
    """
    Yields the parsed riboswitch elements that are encoded in *config*
    in a FASTA-like format.

    Example riboswitch in FASTA-like format: ::

        >Aptamer|38,103|bound
        ((((((((((........(((((....)))))..(((((...........)))))))))))))))
        NNNNNNNCCUAAAACAUACCAGAGAAAUCUGGAGAGGUGAAGAAUACGACCACCUAGGNNNNNNN
        >Aptamer|50,103|unbound
        ......(((((....)))))..(((((...........)))))..........
        AACAUACCAGAGAAAUCUGGAGAGGUGAAGAAUACGACCACCUAGGNNNNNNN
        >Access_constraint|50,56
        AACAUA
        >Hairpin|15,32|unbound
        (((((((...)))))))
        >Hairpin|0,17|bound
        (((((((...)))))))
        >Restriction_site|0,10
        AUGAGUAUGU
        >Context_front|-45,0
        AAGCUAUACCAAGCAUACAAUCAACUCCAAGCUAGAUCUCUUAAG
        >Context_back|103,148
        AUCUAGCGCUGGUACCAUCCCAUUUAACUGUAAGAAGAAUUGCAC

    Read file: ::

        >>> from riboconstruct.riboswitch import iter_riboswitch
        >>> with open("example_file") as config:
        ...     for element in iter_riboswitch(config):
        ...             print element

    :func:`get_riboswitch_from_config_file` can be used to read a
    FASTA-like riboswitch file and return a :class:`Riboswitch`.
    """
    def next():
        seq = config_iter.next().strip()
        while seq[0] == '#':
            seq = config_iter.next().strip()
        if seq[0] == '>':
            raise ValueError("Config has wrong format.")
        return seq

    config_iter = iter(config)
    re_type = None
    for line in config_iter:
        line = line.strip()
        if not len(line):
            continue
        if line[0] == '>' and not re_type:
            line = line[1:].split('|')
            if len(line) < 2:
                raise ValueError("Config has wrong format.")
            re_type = getattr(rs_e.Type, line[0].strip().lower())
            pos = line[1].split(',')
            pos = int(pos[0]), int(pos[1])
            if (re_type == rs_e.Type.aptamer or
                re_type == rs_e.Type.hairpin):
                # catch these two exception to raise a ValueError which is
                # better suited
                try:
                    state = getattr(rs_e.State, line[2].strip().lower())
                except (AttributeError, IndexError):
                    raise ValueError("Config has wrong format.")
        elif line[0] == '#' or line[0] == '':
            continue
        elif re_type is not None:
            if line[0] == '>':
                raise ValueError('Config has wrong format.')
            if re_type == rs_e.Type.aptamer:
                struct = rna.Structure(line)
                seq = rna.IUPACSequence(next().upper())
                if pos[1] - pos[0] != len(struct) != len(seq):
                    raise ValueError("Aptamer element has wrong size.")
                yield rs_e.Aptamer(state, pos, struct, seq)
            elif re_type == rs_e.Type.hairpin:
                struct = rna.Structure(line)
                if pos[1] - pos[0] != len(struct):
                    raise ValueError("Hairpin element has wrong size.")
                yield rs_e.Hairpin(state, pos, struct)
            elif re_type == rs_e.Type.target_site:
                seq = rna.IUPACSequence(line.upper())
                if pos[1] - pos[0] != len(seq):
                    raise ValueError("Restriction site element has wrong size.")
                yield rs_e.TargetSite(pos, seq)
            elif re_type == rs_e.Type.context_front:
                seq = rna.IUPACSequence(line)
                if pos[1] - pos[0] != len(seq):
                    raise ValueError("Context front element has wrong size.")
                yield rs_e.ContextFront(pos, seq)
            elif re_type == rs_e.Type.context_back:
                seq = rna.IUPACSequence(line)
                if pos[1] - pos[0] != len(seq):
                    raise ValueError("Context back element has wrong size.")
                yield rs_e.ContextBack(pos, seq)
            elif re_type == rs_e.Type.seq_constraint:
                seq = rna.IUPACSequence(line.upper())
                if pos[1] - pos[0] != len(seq):
                    raise ValueError("Sequence constraint element has wrong "
                                     "size.")
                yield rs_e.Sequence(pos, seq)
            elif re_type == rs_e.Type.access_constraint:
                yield rs_e.AccessConstraint(pos)
            else:
                raise TypeError("Riboswitch config file contains invalid "
                                "riboswitch element.")
            re_type = None
        else:
            raise ValueError("Config has wrong format.")


def get_riboswitch_from_config_file(config_file, validator=None):
    """
    Returns a riboswitch instance by parsing the elements from a
    FASTA-like format using :func:`iter_riboswitch`.
    """
    riboswitch = Riboswitch()
    # read initial settings from config
    with open(config_file) as config:
        map(riboswitch.add, iter_riboswitch(config))
    if validator is not None:
        validator.validate(riboswitch)
    return riboswitch


EXPRESSION_PLATFORM_MAX_LEN = 50

PATTERNS = [None, None]

# in the bound case:
# due to ligand binding the aptamer is in a certain conformation forcing the
# hairpin to fold into a certain structure, i.e. use the aptamer as structural
# constraint and check if the hairpin folds into the modelled structure
PATTERNS[rs_e.State.bound] = (
    "({}|%%s%%s#%%s%s&%%s%%s|1)" % unichr(167).encode("utf-8"))

# in the unbound case:
# both hairpin and aptamer should fold into the modelled structure, i.e. use NO
# structural constraint and check if the aptamer and the hairpin fold into the
# modelled structures
PATTERNS[rs_e.State.unbound] = (
    "({}|%%s%%s#%%s%s%%s#%%s%%s|1)" % unichr(167).encode("utf-8"))

PATTERNS = tuple(PATTERNS)


def get_riboswitch_from_str(riboswitch_str):
    """
    Parse a riboswitch string representation and return the
    corresponding riboswitch instance.
    """
    riboswitch = Riboswitch()
    elements_iter = iter(riboswitch_str.split(';'))

    # first two elements are the positions
    positions = elements_iter.next().split(', ')
    riboswitch.pos = [int(positions[0][1:]), int(positions[1][:-1])]
    positions = elements_iter.next().split(', ')
    riboswitch.pos_riboswitch = [int(positions[0][1:]), int(positions[1][:-1])]

    for element in elements_iter:
        element = element[1:-1].split(' ')  # get rid of '(' and ')' and split
        type_ = element[0].rsplit('_', 1)  # split type from state
        try:
            type_, state = type_
        except ValueError:  # element has no state
            type_ = type_[0]
            pos = int(element[1][1:-1]), int(element[2][:-1])
            seq = rna.Sequence(element[3][1:-1])
            if type_ == rs_e.TargetSite.ident:
                riboswitch.add(rs_e.TargetSite(pos, seq))
            elif type_ == rs_e.ContextFront.ident:
                riboswitch.add(rs_e.ContextFront(pos, seq))
            elif type_ == rs_e.ContextBack.ident:
                riboswitch.add(rs_e.ContextBack(pos, seq))
            elif type_ == rs_e.AccessConstraint.ident:
                riboswitch.add(rs_e.AccessConstraint(pos, seq))
        else:  # element has a state
            state = (rs_e.State.bound
                     if state == rs_e.State.get_str(rs_e.State.bound)
                     else rs_e.State.unbound)
            pos = int(element[1][1:-1]), int(element[2][:-1])
            struct = rna.Structure(element[3][1:-1])
            if type_ == rs_e.Aptamer.ident:
                seq = rna.Sequence(element[4][1:-1])
                riboswitch.add(rs_e.Aptamer(state, pos, struct, seq))
            elif type_ == rs_e.Hairpin.ident:
                riboswitch.add(rs_e.Hairpin(state, pos, struct))
    return riboswitch


def get_constraints((start, end), riboswitch_elements):
    def set_base(pos, new_base):
        if new_base == rna.BASES[rna.BaseId.UNSPEC]:
            return
        if seq[pos] == rna.BASES[rna.BaseId.UNSPEC]:
            seq[pos] = new_base
        # TODO: remove later since this should not happen
        elif seq[pos] != new_base:
            raise ValueError("Overlapping bases at position %i not valid %s "
                             "vs. %s" % (pos, seq[pos], new_base))

    # set new_struct_elem if the elem has not been specified before
    def set_struct_elem(state, pos, new_struct_elem):
        if structs[state][pos] == rna.STRUCT_ELEM_UNSPEC:
            structs[state][pos] = new_struct_elem

    length = end - start
    seq = [rna.BASES[rna.BaseId.UNSPEC]] * length
    structs = ([rna.STRUCT_ELEM_UNSPEC] * length,
               [rna.STRUCT_ELEM_UNSPEC] * length)

    for rs_elem in riboswitch_elements:
        # clip rs_elem position into the start-end range
        rs_elem_start = rs_elem.pos[0] if rs_elem.pos[0] > start else start
        rs_elem_end = rs_elem.pos[1] if rs_elem.pos[1] < end else end
        i, j = rs_elem_start - start, rs_elem_end - start

        if (rs_elem.type == rs_e.Type.target_site or
            rs_elem.type == rs_e.Type.seq_constraint):
            for pos in xrange(i, j):
                set_base(pos, rs_elem.seq[pos - i])
        elif rs_elem.type == rs_e.Type.hairpin:
            for pos in xrange(i, j):
                set_struct_elem(rs_elem.state, pos, rs_elem.struct[pos - i])
        elif rs_elem.type == rs_e.Type.aptamer:
            for pos in xrange(i, j):
                set_base(pos, rs_elem.seq[pos - i])
                set_struct_elem(rs_elem.state, pos, rs_elem.struct[pos - i])
        elif rs_elem.type == rs_e.Type.access_constraint:
            # The sequence of the access constraint should fit to the aptamer
            # sequence it is overlapping with (has been checked before), i.e.
            # do nothing here
            pass
        else:  # context_front or context_back
            for pos in xrange(i, j):
                set_base(pos, rs_elem.seq[pos - i])

    return (''.join(structs[0]), ''.join(structs[1])), ''.join(seq)


# INFO: the currenty riboswitch model only supports ONE hairpin loop (for each
# state) although some functionalities are already implemented to support more
# than one hairpin

class Riboswitch(object):
    """
    A class to represent riboswitch instances which are defined by their
    riboswitch elements (:class:`~riboconstruct.riboswitch.element`).
    """

    __all__ = ["copy", "add", "remove", "replace", "reset_repr",
               "get_constraints", "get_constraints_riboswitch",
               "get_len", "get_len_riboswitch"]

    def __init__(self):
        self.elements = set()
        self._elements = [[] for _ in xrange(rs_e.Type.count)]

        maxint = sys.maxint
        minint = -maxint - 1
        self.pos = [maxint, minint]
        self.pos_riboswitch = [maxint, minint]

    def __hash__(self):
        """
        Returns the hash value of a riboswitch which is based on the
        riboswitch's string representation.
        """
        return hash(repr(self))

    def __eq__(self, other):
        """
        Compares two riboswitches based on their string representation.
        """
        return repr(self) == repr(other)

    def __repr__(self):
        """
        Returns the riboswitch instance string representation.

        The representation is computed only once and then stored in
        ``Instance._str``. ``Instance.reset_repr()`` is used to reset
        the representation whenever the riboswitch instance is changed.
        """
        try:
            return self._str
        except AttributeError:
            sorted_elements = sorted(self.elements,
                                     key=lambda e: (e.type, e.state, e.pos))
            self._str = (repr(self.pos) + ';' +
                         repr(self.pos_riboswitch) + ';' +
                         ';'.join(repr(e) for e in sorted_elements))
            return self._str

    def copy(self):
        """Returns a copy of the riboswitch."""
        new = Riboswitch()
        new.elements = set(self.elements)
        # TODO: remove later
        new._elements = (
            list(list(element)
                 for element in self._elements))
        new.pos = copy.deepcopy(self.pos)
        new.pos_riboswitch = copy.deepcopy(self.pos_riboswitch)
        return new

    def add(self, element):
        """Add a riboswitch element to the instance."""
        self.reset_repr()
        # TODO: remove later
        self._elements[element.type].append(element)
        self.elements.add(element)
        # check the positions of the riboswitch instance
        if element.type == rs_e.Type.context_front:
            self.pos[0] = element.pos[0]
        elif element.type == rs_e.Type.context_back:
            self.pos[1] = element.pos[1]
        elif element.type == rs_e.Type.aptamer:
            if element.pos[1] > self.pos[1]:
                self.pos[1] = element.pos[1]
                self.pos_riboswitch[1] = element.pos[1]
            elif element.pos[1] > self.pos_riboswitch[1]:
                self.pos_riboswitch[1] = element.pos[1]
        else:  # hairpin, seq_constraint, access_constraint, target_site
            if element.pos[0] < self.pos[0]:
                self.pos[0] = element.pos[0]
                self.pos_riboswitch[0] = element.pos[0]
            elif element.pos[0] < self.pos_riboswitch[0]:
                self.pos_riboswitch[0] = element.pos[0]

    def remove(self, element):
        """Remove a riboswitch *element* from the instance."""
        try:
            del self.elements[element]
        except KeyError:
            pass
        try:
            # remove old if it exists
            self._elements[rs_elem.type].remove(rs_elem)
        except ValueError:
            pass

    def replace(self, old, new):
        """Replace an *old* riboswitch element by a *new* one."""
        self.remove(old)
        self.add(new)

    def get_accessibility_positions(self):
        # TODO: HACK: this hole method is really hacky; replace/remove somehow!
        start = self.pos[0]
        rs = list(self._elements[rs_e.Type.target_site][0].pos)
        rs[0] += 3  # HACK!!!
        try:
            ac = self._elements[rs_e.Type.access_constraint][0].pos
        except IndexError:
            raise ValueError("Access constraint should be available!")
        a_b = sorted(self._elements[rs_e.Type.aptamer])[rs_e.State.bound]
        h_0, h_1 = self._elements[rs_e.Type.hairpin]
        l = max(a_b.pos[1] - a_b.pos[0],
                h_0.pos[1] - h_0.pos[0],
                h_1.pos[1] - h_1.pos[0]) + 2
        return (rs[0] - start, ac[0] - start), (rs[1] - rs[0], ac[1] - ac[0]), l

    def get_constraints(self):
        """
        Get the structural constraints for both conformations and the
        sequential constraint defined by the riboswitch's elements plus
        the sequential constraints of its contexts.
        """
        return get_constraints(self.pos, self.elements)

    def get_constraints_riboswitch(self):
        """
        Get the structural constraints for both conformations and the
        sequential constraint defined by the riboswitch's elements
        (without any context).
        """
        return get_constraints(self.pos_riboswitch,
                               (element for element in self.elements if
                                (element.type != rs_e.Type.context_front and
                                 element.type != rs_e.Type.context_back)))

    def reset_repr(self):
        """
        Reset the string representation stored in ``Instance._str``.
        """
        try:
            del self._str
        except AttributeError:
            pass
