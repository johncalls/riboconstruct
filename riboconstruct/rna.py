import re

from .helper import enum

# Boltzmann gas constant in kcal/mol
RT = 0.616
# Penalty for hairpin of size 3 consisting only of C's
C_TRI_LOOP = 1.4
# Penality for hairpin loop with terminal A-U or G-U
TERMINAL_AU = 0.5
# Minimum size of a hairpin
HAIRPIN_MIN_SIZE = 3
# Maximum size of a hairpin
HAIRPIN_MAX_SIZE = 30


BASES = ('A', 'C', 'G', 'U',
         'N')

IUPAC = set(('A', 'C', 'G', 'U',
             'N',
             'R', 'Y', 'S', 'W', 'K', 'M',
             'B', 'D', 'H', 'V'))

BASEPAIRS = ((0, 3), (1, 2), (2, 1), (2, 3), (3, 0), (3, 2))

STRUCT_ELEM_UNSPEC = '.'


BaseId = enum('A', 'C', 'G', 'U', 'N', UNSPEC=4, count=4)
"""
Struct-like representation of RNA bases as id, plus an additional
wildcard base matching the other bases. (See
:func:`riboconstruct.helper.enum` for details.)

**A** = 0,
**C** = 1,
**G** = 2,
**U** = 3,
**N** = 4

**UNSPEC** = 4

**count** = 4 (number of *real* bases)
"""


class IUPAC_Id:
    """Struct-like representation of RNA bases in IUPAC format as id."""
    count = len(IUPAC)

    A, C, G, U, N, R, Y, S, W, K, M, B, D, H, V = xrange(count)


class BasepairId:
    """
    Struct-like representation of RNA basepairs as id.

    (**A**-**U**) = 0,
    (**C**-**G**) = 1,
    (**G**-**C**) = 2,
    (**G**-**U**) = 3,
    (**U**-**A**) = 4,
    (**U**-**G**) = 5
    """
    count = len(BASEPAIRS)

    AU, CG, GC, GU, UA, UG = xrange(count)


class StructType:
    """
    Struct-like representation of the different structural elements of
    a RNA structure.

    **INTERIOR** = 0,
    **BULGE** = 1,
    **HAIRPIN** = 2,
    **STACKING** = 3,
    **MULTILOOP** = 4,
    **EXTERIOR** = 5
    """
    count = 6

    INTERIOR, BULGE, HAIRPIN, STACKING, MULTILOOP, EXTERIOR = xrange(count)


def valid_bp_b_ids(b_id_i, b_id_j):
    """
    Returns whether the base ids *b_id_i* and *b_id_j* can form a valid
    basepair.
    """
    if b_id_i == BaseId.A:
        return b_id_j == BaseId.U
    elif b_id_i == BaseId.C:
        return b_id_j == BaseId.G
    elif b_id_i == BaseId.G:
        return b_id_j == BaseId.C or b_id_j == BaseId.U
    elif b_id_i == BaseId.U:
        return b_id_j == BaseId.A or b_id_j == BaseId.G
    return False


def valid_bp(b_i, b_j):
    """
    Returns whether the bases *b_i* and *b_j* can form a valid basepair.
    """
    if b_i == BASES[BaseId.A]:
        return b_j == BASES[BaseId.U]
    elif b_i == BASES[BaseId.C]:
        return b_j == BASES[BaseId.G]
    elif b_i == BASES[BaseId.G]:
        return b_j == BASES[BaseId.C] or b_j == BASES[BaseId.U]
    elif b_i == BASES[BaseId.U]:
        return b_j == BASES[BaseId.A] or b_j == BASES[BaseId.G]
    return False


# IUPAC nucleotide code -> base(s)
#
# A -> Adenine
# C -> Cytosine
# G -> Guanine
# U -> Uracil
# R -> A or G
# Y -> C or U
# S -> G or C
# W -> A or U
# K -> G or U
# M -> A or C
# B -> C or G or U
# D -> A or G or U
# H -> A or C or U
# V -> A or C or G
# N -> any base
def valid_IUPAC(iu):
    """
    Returns whether base id *ui* is a valid IUPAC base id.

    When a IUPAC base *b* instead of an id is given, use something
    like: ::

        b_id = getattr(IUPAC_Id, b)
        valid_IUPAC(b_id)
    """
    return (iu == IUPAC_Id.A or iu == IUPAC_Id.C or iu == IUPAC_Id.G or
            iu == IUPAC_Id.U or
            iu == IUPAC_Id.R or iu == IUPAC_Id.Y or iu == IUPAC_Id.S or
            iu == IUPAC_Id.W or iu == IUPAC_Id.K or iu == IUPAC_Id.M or
            iu == IUPAC_Id.B or iu == IUPAC_Id.D or iu == IUPAC_Id.H or
            iu == IUPAC_Id.V or
            iu == IUPAC_Id.N)


def base_valid_IUPAC(iu, b_id):
    """
    Returns whether base id *b_id* suffices the IUPAC constraint *iu*.
    """
    if b_id == BaseId.A:
        if (iu == IUPAC_Id.A or iu == IUPAC_Id.M or iu == IUPAC_Id.R or
            iu == IUPAC_Id.W or iu == IUPAC_Id.V or iu == IUPAC_Id.H or
            iu == IUPAC_Id.D or iu == IUPAC_Id.N):
            return True
        else:
            return False
    elif b_id == BaseId.C:
        if (iu == IUPAC_Id.C or iu == IUPAC_Id.M or iu == IUPAC_Id.S or
            iu == IUPAC_Id.Y or iu == IUPAC_Id.V or iu == IUPAC_Id.H or
            iu == IUPAC_Id.B or iu == IUPAC_Id.N):
            return True
        else:
            return False
    elif b_id == BaseId.G:
        if (iu == IUPAC_Id.G or iu == IUPAC_Id.R or iu == IUPAC_Id.S or
            iu == IUPAC_Id.K or iu == IUPAC_Id.V or iu == IUPAC_Id.D or
            iu == IUPAC_Id.B or iu == IUPAC_Id.N):
            return True
        else:
            return False
    elif b_id == BaseId.U:
        if (iu == IUPAC_Id.U or iu == IUPAC_Id.W or iu == IUPAC_Id.Y or
            iu == IUPAC_Id.K or iu == IUPAC_Id.H or iu == IUPAC_Id.D or
            iu == IUPAC_Id.B or iu == IUPAC_Id.N):
            return True
        else:
            return False
    raise ValueError('Invalid IUPAC symbol %i' % b_id)


def check_struct_seq_match(struct, seq):
    """
    Returns whether IUPAC sequence *seq* can form structure *struct*.
    """
    for bp_pos_i, bp_pos_j in struct.bp_positions:
        if not valid_bp(seq[bp_pos_i], seq[bp_pos_j]):
            for bp_i, bp_j in BASEPAIRS:
                if (base_valid_IUPAC(getattr(IUPAC_Id, seq[bp_pos_i]),
                                     bp_i) and
                    base_valid_IUPAC(getattr(IUPAC_Id, seq[bp_pos_j]), bp_j)):
                    return
            raise ValueError("Structure and sequence do not match")


class Structure(object):
    """RNA structure representation."""

    __all__ = ["basepairs", "bp_positions"]

    def __init__(self, struct):
        self.struct = list(struct)
        self._basepairs = tuple(self._calc_basepairs())

    @property
    def basepairs(self):
        """
        Returns a tuple storing for each position whether it is forming
        a basepair with another position in the underlying structure or
        not.

        - ``basepairs[i] = j``, if *i* is forming a basepair
          with *j*
        - ``basepairs[i] = None``, if *i* is unpaired

        The information is calculated only once at initilization for
        each :class:`Structure`.
        """
        return self._basepairs

    @property
    def bp_positions(self):
        """
        Returns a tuple containing all basepairs in an ordered manner
        such that for every basepair *j* enclosed by another basepair
        *i*, ``i > j``.

        The position of a specific basepair within this order is the
        *basepair position id* (*bp_pos_id*).

        The order has to be calculated only once for each structure
        instance.
        """
        try:
            return self.__bp_positions
        except AttributeError:
            self.__bp_positions = tuple(self._calc_bp_positions())
            return self.__bp_positions

    def __getattr__(self, attr):
        """Everything else is delegated to the object"""
        return getattr(self.struct, attr)

    def __repr__(self):
        return ''.join(self.struct)

    def __iter__(self):
        return iter(self.struct)

    def __getitem__(self, pos):
        return self.struct[pos]

    def __len__(self):
        return len(self.struct)

    def __eq__(self, a):
        return repr(self) == repr(a)

    def _calc_bp_positions(self):
        """Returns an order over the basepairs."""
        bp_positions = ((i, j) for i, j in enumerate(self.basepairs)
                        if j is not None and i < j)
        # sort base pairs: (i1, j1) < (i2, j2) iff i1 > i2
        return sorted(bp_positions, key=lambda bp: bp[0], reverse=True)

    def _calc_basepairs(self):
        """
        Calculates for each position within the structure whether it is
        forming a basepair (and with which other position) or not.
        """
        basepairs = {}
        brackets = []
        within_hairpin = False
        unpaired_bases_counter = 0
        for i, s in enumerate(self):
            if s == '.':
                basepairs[i] = None
                if within_hairpin:
                    unpaired_bases_counter += 1
                    if unpaired_bases_counter > HAIRPIN_MAX_SIZE:
                        raise ValueError(
                            "Structure contains hairpin with size > %i." %
                            HAIRPIN_MAX_SIZE)
            elif s == '(':
                brackets.append(i)
                unpaired_bases_counter = 0
                within_hairpin = True
            elif s == ')':
                within_hairpin = False
                unpaired_bases_counter = 0
                try:
                    basepairs[i] = brackets.pop()
                except:
                    raise ValueError(
                        "No matching opening bracket in structure '%s'." % self)
                if i - basepairs[i] - 1 < HAIRPIN_MIN_SIZE:
                    raise ValueError(
                        "Structure contains hairpin with size < %i." %
                        HAIRPIN_MIN_SIZE)
                basepairs[basepairs[i]] = i
            else:
                raise ValueError("'%s' is not a valid structure element." %
                                 s)
        if len(brackets):
            raise ValueError("No matching closing bracket in structure "
                             "'%s'." % self)
        return [basepairs[k] for k in sorted(basepairs.iterkeys())]


class Sequence(object):
    """RNA sequence representation."""

    base_check = re.compile(r"[^%s]" % ''.join(BASES[:-1]), re.IGNORECASE)

    def __init__(self, seq):
        if type(seq) is not Sequence and self.base_check.search(seq):
            raise TypeError("Not a valid RNA sequence.")
        self.seq = list(seq)

    def __repr__(self):
        return ''.join(self.seq)

    def __eq__(self, seq):
        return all(i == j for i, j in zip(self, seq))

    def __hash__(self):
        return hash(repr(self.seq))

    def __len__(self):
        return len(self.seq)

    def __iter__(self):
        return iter(self.seq)

    def __getitem__(self, pos):
        return self.seq[pos]

    def __setitem__(self, pos, base):
        if base not in BASES:
            raise TypeError("Not a valid RNA base.")
        self.seq[pos] = base


class IUPACSequence(Sequence):
    """
    IUPAC sequence representation.

    Example: ::

        >>> s = IUPACSequence("GGVAUUCY")
        >>> print s
        GGVAUUCY
        >>> print s[2]
        V
        >>> IUPACSequence("GGIAUUCY")
        TypeError: No valid IUPAC sequence.
    """

    base_check = re.compile(r"[^%s]" % ''.join(IUPAC), re.IGNORECASE)

    def __init__(self, seq):
        if self.base_check.search(seq):
            raise TypeError("No valid IUPAC sequence.")
        self.seq = seq

    def __setitem__(self, pos, iupac_base):
        if iupac_base not in IUPAC:
            raise TypeError("No valid IUPAC base.")
        self.seq[pos] = iupac_base
