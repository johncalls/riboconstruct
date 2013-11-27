import types

from .. import rna
from ..helper import enum


State = enum("unbound", "bound", count=2)
"""
Struct-like representation of the riboswitch states. (See
:func:`riboconstruct.helper.enum` for details.)

**unbound** = 0, **bound** = 1

**count** = 2

.. function:: get_str(state)

    Returns the string representation (``ub``, ``b``) of a certain
    state.

Example: ::

    >>> from riboconstruct.riboswitch import element
    >>> state = element.State.unbound
    >>> element.State.get_str(state)
    'ub'
"""

State.get_str = types.MethodType(lambda _, state:
                                 state == State.unbound and "ub" or "b",
                                 State)


Type = enum("aptamer", "hairpin", "target_site", "context_front",
            "context_back", "seq_constraint", "access_constraint", count=7)
"""
Struct-like representation of the riboswitch element types. (See
:func:`riboconstruct.helper.enum` for details.)

**aptamer** = 0,
**hairpin** = 1,
**target_site** = 2,
**context_front** = 3,
**context_back** = 4,
**seq_constraint** = 5,
**access_constraint** = 6
"""


EXPRESSION_PLATFORM_TYPES = set((Type.hairpin, Type.target_site))


class Element(object):
    """Base class of the different riboswitch element types."""

    __all__ = ["type", "ident"]

    def __init__(self, state, pos, struct, seq):
        self.state = state
        self.pos = pos
        self.struct = struct
        self.seq = seq


class Aptamer(Element):
    """Aptamer element."""

    __all__ = ["type", "ident", "__repr__"]

    ident = "A"
    type = Type.aptamer

    def __init__(self, state, pos, struct, seq):
        if len(struct) != len(seq) != pos[1] - pos[0]:
            raise ValueError("The lengths of the aptamer given by its "
                             "position, structure and sequence do not match.")
        super(Aptamer, self).__init__(state, pos, struct, seq)

    def __repr__(self):
        return '(%s_%s %s "%s" "%s")' % (self.ident, State.get_str(self.state),
                                         self.pos, str(self.struct),
                                         str(self.seq))


class Hairpin(Element):
    """Hairpin element."""

    __all__ = ["type", "ident", "__repr__", "shift",
               "insert_bp_after", "insert_bp_before",
               "remove_bp",
               "insert_b", "remove_b"]

    ident = "H"
    type = Type.hairpin

    def __init__(self, state, pos, struct):
        if len(struct) != pos[1] - pos[0]:
            raise ValueError("The lengths of the hairpin given by its "
                             "position and structure do not match.")
        super(Hairpin, self).__init__(state, pos, struct, None)

    def __repr__(self):
        return '(%s_%s %s "%s")' % (self.ident, State.get_str(self.state),
                                    self.pos, str(self.struct))

    def shift(self, num_pos=1):
        """
        Returns a new hairpin shifted upstream by *num_pos* positions,
        if ``num_pos > 0``. If ``num_pos < 0``, the hairpin is shifted
        downstream.
        """
        pos = self.pos[0] + num_pos, self.pos[1] + num_pos
        return Hairpin(self.state, pos, self.struct)

    def insert_bp_after(self, bp_pos_id=0):
        """
        Returns a new hairpin with a new basepair inserted after
        *bp_pos_id*.
        """
        bp_pos_i, bp_pos_j = self.struct.bp_positions[bp_pos_id]
        new_bp_pos_i, new_bp_pos_j = bp_pos_i + 1, bp_pos_j + 1
        # copy
        struct = self.struct[:]
        struct.insert(new_bp_pos_i, '(')
        struct.insert(new_bp_pos_j, ')')
        pos = self.pos[0] - 1, self.pos[1] + 1
        return Hairpin(self.state, pos, rna.Structure(struct))

    def insert_bp_before(self, bp_pos_id=0):
        """
        Returns a new hairpin with a new basepair inserted before
        *bp_pos_id*.
        """
        bp_pos_i, bp_pos_j = self.struct.bp_positions[bp_pos_id]
        new_bp_pos_i, new_bp_pos_j = bp_pos_i - 1, bp_pos_j + 2
        if bp_pos_i == 0:
            new_bp_pos_i = 0
        else:
            new_bp_pos_i = bp_pos_i - 1
        # copy
        struct = self.struct[:]
        struct.insert(new_bp_pos_i, '(')
        struct.insert(new_bp_pos_j, ')')
        pos = self.pos[0] - 1, self.pos[1] + 1
        return Hairpin(self.state, pos, rna.Structure(struct))

    def remove_bp(self, bp_pos_id=0):
        """
        Returns a new hairpin with the basepair removed at *bp_pos_id*.
        """
        bp_pos_i, bp_pos_j = self.struct.bp_positions[bp_pos_id]
        # copy
        struct = self.struct[:]
        del struct[bp_pos_i]
        del struct[bp_pos_j - 1]
        pos = self.pos[0] + 1, self.pos[1] - 1
        return Hairpin(self.state, pos, rna.Structure(struct))

    def insert_b(self, pos):
        """
        Returns a new hairpin with a unpaired base inserted at *pos*.
        """
        # copy
        struct = self.struct[:]
        struct.insert(pos, '.')
        pos = self.pos[0], self.pos[1] + 1
        return Hairpin(self.state, pos, rna.Structure(struct))

    def remove_b(self, pos):
        """
        Returns a new hairpin with a unpaired base removed at *pos*.
        """
        # copy
        struct = self.struct[:]
        del struct[pos]
        pos = self.pos[0], self.pos[1] - 1
        return Hairpin(self.state, pos, rna.Structure(struct))


class TargetSite(Element):
    """Target site element."""

    __all__ = ["type", "ident", "__repr__", "shift"]

    ident = "TS"
    type = Type.target_site

    def __init__(self, pos, seq):
        if len(seq) != pos[1] - pos[0]:
            raise ValueError("The lengths of the target site given by its "
                             "position and sequence do not match.")
        super(TargetSite, self).__init__(None, pos, None, seq)

    def __repr__(self):
        return '(%s %s "%s")' % (self.ident, self.pos, str(self.seq))

    def shift_down(self, num_pos=1):
        pos = self.pos[0] - num_pos, self.pos[1] - num_pos
        return TargetSite(pos, self.seq)

    def shift_up(self, num_pos=1):
        pos = self.pos[0] + num_pos, self.pos[1] + num_pos
        return TargetSite(pos, self.seq)

    def shift(self, num_pos=1):
        """
        Returns a new target site element shifted upstream by *num_pos*
        positions, if ``num_pos > 0``. If ``num_pos < 0``, the target
        site is shifted downstream.
        """
        pos = self.pos[0] + num_pos, self.pos[1] + num_pos
        return TargetSite(pos, self.seq)


class Context(Element):
    """Base class for the context elements."""

    __all__ = ["type", "ident", "__repr__", "shift"]

    def __init__(self, pos, seq):
        if len(seq) != pos[1] - pos[0]:
            raise ValueError("The lengths of the context given by its "
                             "position and sequence do not match.")
        super(Context, self).__init__(None, pos, None, seq)

    def __repr__(self):
        return '(%s %s "%s")' % (self.ident, self.pos, str(self.seq))


class ContextFront(Context):
    """Context element located at the 5' end of the riboswitch."""

    __all__ = ["type", "ident", "__repr__", "shift"]

    ident = "Cf"
    type = Type.context_front

    def __init__(self, pos, seq):
        super(ContextFront, self).__init__(pos, seq)

    def shift_down(self, num_pos=1):
        pos = self.pos[0] - num_pos, self.pos[1] - num_pos
        return ContextFront(pos, self.seq)

    def shift_up(self, num_pos=1):
        pos = self.pos[0] + num_pos, self.pos[1] + num_pos
        return ContextFront(pos, self.seq)

    def shift(self, num_pos=1):
        """
        Returns a context front element shifted upstream by *num_pos*
        positions, if ``num_pos > 0``. If ``num_pos < 0``, the context
        front is shifted downstream.
        """
        pos = self.pos[0] + num_pos, self.pos[1] + num_pos
        return ContextFront(pos, self.seq)


class ContextBack(Context):
    """Context element located at the 3' end of the riboswitch."""

    __all__ = ["type", "ident", "__repr__", "shift"]

    ident = "Cb"
    type = Type.context_back

    def __init__(self, pos, seq):
        super(ContextBack, self).__init__(pos, seq)

    def shift(self, num_pos=1):
        """
        Returns a context back element shifted upstream by *num_pos*
        positions, if ``num_pos > 0``. If ``num_pos < 0``, the context
        back is shifted downstream.
        """
        pos = self.pos[0] + num_pos, self.pos[1] + num_pos
        return ContextBack(pos, self.seq)


class SequenceConstraint(Element):
    """Sequence constraint element."""

    ident = "SC"
    type = Type.seq_constraint

    def __init__(self, pos, seq):
        if len(seq) != pos[1] - pos[0]:
            raise ValueError("The lengths of the sequential constraint given "
                             "by its position and sequence do not match.")
        super(SequenceConstraint, self).__init__(None, pos, None, seq)

    def __repr__(self):
        return '(%s %s "%s")' % (self.ident, self.pos, str(self.seq))


class AccessConstraint(Element):
    """Accessibility constraint element."""

    ident = "AC"
    type = Type.access_constraint

    def __init__(self, pos):
        super(AccessConstraint, self).__init__(None, pos, None, None)

    def __repr__(self):
        return '(%s %s)' % (self.ident, self.pos)
