import types

from .. import rna
from ..helper import enum


State = enum("unbound", "bound", count=2)
"""
Struct-like representation of the riboswitch states.

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


class Type:
    count = 7

    (aptamer,
     hairpin,
     restriction_site,
     context_front,
     context_back,
     seq,
     access_constraint) = xrange(count)


FUNCTIONAL_SITE_TYPES = set((Type.hairpin, Type.restriction_site))


class RiboswitchElement(object):
    def __init__(self, m_type, state, pos, struct, seq=None):
        self.type = m_type
        self.state = state
        self.pos = pos
        self.struct = struct
        self.seq = seq


class Aptamer(RiboswitchElement):
    def __init__(self, state, pos, struct, seq):
        super(Aptamer, self).__init__(
            Type.aptamer, state, pos, struct, seq)

    def __repr__(self):
        return "(A_%s (%i,%i) '%s' '%s')" % (
            "b" if self.state == State.bound else "ub",
            self.pos[0], self.pos[1], str(self.struct), str(self.seq))


class Hairpin(RiboswitchElement):
    def __init__(self, state, pos, struct):
        super(Hairpin, self).__init__(
            Type.hairpin, state, pos, rna.Structure(struct))

    def __repr__(self):
        return "(H_%s (%i,%i) '%s')" % (
            "b" if self.state == State.bound else "ub",
            self.pos[0], self.pos[1], str(self.struct))

    def shift_down(self, num_pos=1):
        pos = self.pos[0] - num_pos, self.pos[1] - num_pos
        return Hairpin(self.state, pos, self.struct)

    def shift_up(self, num_pos=1):
        pos = self.pos[0] + num_pos, self.pos[1] + num_pos
        return Hairpin(self.state, pos, self.struct)

    def insert_bp_after(self, bp_pos_id):
        bp_pos_i, bp_pos_j = self.struct.bp_positions[bp_pos_id]
        new_bp_pos_i, new_bp_pos_j = bp_pos_i + 1, bp_pos_j + 1
        # copy
        struct = self.struct[:]
        struct.insert(new_bp_pos_i, '(')
        struct.insert(new_bp_pos_j, ')')
        pos = self.pos[0] - 1, self.pos[1] + 1
        return Hairpin(self.state, pos, struct)

    def insert_bp_before(self, bp_pos_id):
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
        return Hairpin(self.state, pos, struct)

    def remove_bp(self, bp_pos_id):
        bp_pos_i, bp_pos_j = self.struct.bp_positions[bp_pos_id]
        # copy
        struct = self.struct[:]
        del struct[bp_pos_i]
        del struct[bp_pos_j - 1]
        pos = self.pos[0] + 1, self.pos[1] - 1
        return Hairpin(self.state, pos, struct)

    def insert_b(self, pos):
        # copy
        struct = self.struct[:]
        struct.insert(pos, '.')
        pos = self.pos[0], self.pos[1] + 1
        return Hairpin(self.state, pos, struct)

    def remove_b(self, pos):
        # copy
        struct = self.struct[:]
        del struct[pos]
        pos = self.pos[0], self.pos[1] - 1
        return Hairpin(self.state, pos, struct)


class AccessConstraint(RiboswitchElement):
    ident = "Acc"

    def __init__(self, pos, seq):
        super(AccessConstraint, self).__init__(
            Type.access_constraint, None, pos, None, seq)

    def __repr__(self):
        return (
            "(%s (%i,%i) '%s')" % (
                self.ident, self.pos[0], self.pos[1], str(self.seq)))


class RestrictionSite(RiboswitchElement):
    def __init__(self, pos, seq):
        super(RestrictionSite, self).__init__(
            Type.restriction_site, None, pos, None, seq)

    def __repr__(self):
        return "(RS (%i,%i) '%s')" % (self.pos[0], self.pos[1], str(self.seq))

    def shift_down(self, num_pos=1):
        pos = self.pos[0] - num_pos, self.pos[1] - num_pos
        return RestrictionSite(pos, self.seq)

    def shift_up(self, num_pos=1):
        pos = self.pos[0] + num_pos, self.pos[1] + num_pos
        return RestrictionSite(pos, self.seq)


class Context(RiboswitchElement):
    def __init__(self, m_type, pos, struct, seq):
        super(Context, self).__init__(m_type, None, pos, struct, seq)


class ContextFront(Context):
    def __init__(self, pos, struct, seq):
        super(ContextFront, self).__init__(
            Type.context_front, pos, struct, seq)

    def __repr__(self):
        return "(Cf (%i,%i) '%s' '%s')" % (
            self.pos[0], self.pos[1], str(self.struct), str(self.seq))

    def shift_down(self, num_pos=1):
        pos = self.pos[0] - num_pos, self.pos[1] - num_pos
        return ContextFront(pos, self.struct, self.seq)

    def shift_up(self, num_pos=1):
        pos = self.pos[0] + num_pos, self.pos[1] + num_pos
        return ContextFront(pos, self.struct, self.seq)


class ContextBack(Context):
    def __init__(self, pos, struct, seq):
        super(ContextBack, self).__init__(
            Type.context_back, pos, struct, seq)

    def __repr__(self):
        return "(Cb (%i,%i) '%s' '%s')" % (
            self.pos[0], self.pos[1], str(self.struct), str(self.seq))


class SequenceConstraint(RiboswitchElement):
    def __init__(self, pos, seq):
        super(SequenceConstraint, self).__init__(
            Type.seq, None, pos, None, seq)

    def __repr__(self):
        return "(S (%i,%i) '%s')" % (self.pos[0], self.pos[1], str(self.seq))
