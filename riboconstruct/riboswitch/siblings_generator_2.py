from . import element as rs_e

from .. import riboswitch as rs
from .. import rna


HAIRPIN_LOOP_MAX_SIZE = 7
HAIRPIN_LOOP_MIN_SIZE = 3
HAIRPIN_STEM_MIN_SIZE = 5
HAIRPIN_STEM_MAX_SIZE = 30

UNBOUND_H_APT_DIFF = 2


class SimpleSiblingsGenerator(object):
    """
    Idea:
    The restriction site 'atg a gtatgt' marks the downstream end of the
    expression platform. The context front is directly connected, i.e. the
    context front can never leap into the expression platform
    """

    def __init__(self, parent):
        self.parent = parent

        rs_elements = self.parent._elements
        self.hairpins = rs_elements[rs_e.Type.hairpin]
        if self.hairpins[0].state != 0:
            self.hairpins = [self.hairpins[1], self.hairpins[0]]
        self.aptamers = rs_elements[rs_e.Type.aptamer]
        if self.aptamers[0].state != 0:
            self.aptamers = [self.aptamers[1], self.aptamers[0]]
        self.restriction_site = rs_elements[rs_e.Type.restriction_site][0]
        try:
            self.context_front = rs_elements[rs_e.Type.context_front][0]
        except IndexError:
            self.context_front = None

    def iter_siblings(self):
        # alter the hairpins
        for h, h_new in self.iter_altered_hairpins():
            sibling = self.parent.copy()
            sibling.replace(h, h_new)
            yield sibling

        # alter the expression platform (this implies a shift of the restriction
        # site and the context (if available))
        if self.valid_increase_ep():
            sibling = self.parent.copy()
            new_rs = self.restriction_site.shift_down()
            sibling.replace(self.restriction_site, new_rs)
            if self.context_front:
                new_cf = self.context_front.shift_down()
                sibling.replace(self.context_front, new_cf)
            yield sibling
        if self.valid_decrease_ep():
            sibling = self.parent.copy()
            new_rs = self.restriction_site.shift_up()
            sibling.replace(self.restriction_site, new_rs)
            if self.context_front:
                new_cf = self.context_front.shift_up()
                sibling.replace(self.context_front, new_cf)
            yield sibling

    def iter_altered_hairpins(self):
        # alter hairpins
        for hairpin in self.hairpins:
            h_shift_up = self.valid_hairpin_shift_up(hairpin)
            if h_shift_up:
                yield hairpin, hairpin.shift_up()
            h_shift_down = self.valid_hairpin_shift_down(hairpin)
            if h_shift_down:
                yield hairpin, hairpin.shift_down()
            if (h_shift_up and
                h_shift_down and
                len(hairpin.struct.bp_positions) < HAIRPIN_STEM_MAX_SIZE):
                yield hairpin, hairpin.insert_bp_before(0)
            if self.valid_hairpin_remove_bp(hairpin):
                yield hairpin, hairpin.remove_bp(0)
            # alter hairpin loop
            bp_pos_i, bp_pos_j = hairpin.struct.bp_positions[0]
            loop_size = bp_pos_j - bp_pos_i - 1
            if (loop_size < HAIRPIN_LOOP_MAX_SIZE and
                self.valid_increase_loop(hairpin)):
                yield hairpin, hairpin.insert_b(bp_pos_i + 1)
            if (loop_size > HAIRPIN_LOOP_MIN_SIZE and
                self.valid_decrease_loop(hairpin)):
                yield hairpin, hairpin.remove_b(bp_pos_i + 1)

    def valid_increase_ep(self):
        pos = self.restriction_site.pos
        new_pos = pos[0] - 1, pos[1] - 1
        hairpin_0, hairpin_1 = self.hairpins
        # after shifting, the restriction site should not prevent the hairpins
        # from forming a valid bp
        if (not self.valid_struct_seq_matching(hairpin_0.pos, new_pos,
                                               hairpin_0) or
            not self.valid_struct_seq_matching(hairpin_1.pos, new_pos,
                                               hairpin_1)):
            return False
        # after shifting, the bound hairpin should still be 'in touch'
        if not new_pos[1] > self.hairpins[rs_e.State.bound].pos[0]:
            return False
        # the size of the expression platform is bound to a maximum value
        pos = self.parent.pos_functional_site
        return pos[1] - pos[0] < rs.FUNCTIONAL_SITE_MAX_LENGTH

    def valid_decrease_ep(self):
        pos = self.restriction_site.pos
        new_pos = pos[0] + 1, pos[1] + 1
        hairpin_0, hairpin_1 = self.hairpins
        # after shifting, the restriction site should not prevent the
        # hairpins from forming a valid bp
        if (not self.valid_struct_seq_matching(hairpin_0.pos, new_pos,
                                               hairpin_0) or
            not self.valid_struct_seq_matching(hairpin_1.pos, new_pos,
                                               hairpin_1)):
            return False
        # after shifting, the bound hairpin should still be within the
        # expression platform
        if self.hairpins[rs_e.State.bound].pos[0] < new_pos[0]:
            return False
        # after shifting, the unbound hairpin should still be out of range
        return pos[1] < self.hairpins[rs_e.State.unbound].pos[0]

    def valid_hairpin_shift_down(self, hairpin):
        pos = hairpin.pos
        new_pos = pos[0] - 1, pos[1] - 1
        rs_pos = self.restriction_site.pos
        # if the restriction site is long enough: check that it doesn't prevent
        # the hairpin from forming a valid bp
        if (rs_pos[1] - rs_pos[0] >= rna.HAIRPIN_MIN_SIZE + 2 and
            not self.valid_struct_seq_matching(new_pos, rs_pos, hairpin)):
            return False
        # check that after shifting..
        if hairpin.state == rs_e.State.unbound:
            # ..the bound hairpin and the bound aptamer are still in range and
            # the restriction site is still uncovered
            aptamer_b = self.aptamers[rs_e.State.bound]
            if (new_pos[0] <= self.hairpins[rs_e.State.bound].pos[0] or
                pos[0] <= self.restriction_site.pos[1] or
                aptamer_b.pos[0] - new_pos[1] > UNBOUND_H_APT_DIFF):
                return False
        else:  # hairpin.state == bound
            # ..the unbound hairpin is still in range and the restriction site
            # is still covered
            if (pos[1] <= self.hairpins[rs_e.State.unbound].pos[0] or
                new_pos[1] <= self.restriction_site.pos[0]):
                return False
        # the shifted hairpin should still be within the expression platform
        return new_pos[0] >= rs_pos[0]

    def valid_hairpin_shift_up(self, hairpin):
        pos = hairpin.pos
        new_pos = pos[0] + 1, pos[1] + 1
        rs_pos = self.restriction_site.pos
        # if the restriction site is long enough: check that it doesn't prevent
        # the hairpin from forming a valid bp
        if (rs_pos[1] - rs_pos[0] >= rna.HAIRPIN_MIN_SIZE + 2 and
            not self.valid_struct_seq_matching(new_pos, rs_pos, hairpin)):
            return False
        # check that after shifting..
        if hairpin.state == rs_e.State.bound:
            # ..the unbound hairpin is still in range and the restriction site
            # is still covered with at least one base
            if (new_pos[0] >= self.hairpins[rs_e.State.unbound].pos[0] or
                new_pos[0] >= self.restriction_site.pos[1]):
                return False
        else:  # hairpin.state == unbound
            # ..the bound hairpin is still in range
            if new_pos[0] > self.hairpins[rs_e.State.bound].pos[1]:
                return False
        # check intersections with the aptamer in the same state
        return pos[1] < self.aptamers[hairpin.state].pos[0]

    def valid_increase_loop(self, hairpin):
        pos = hairpin.pos
        new_pos = pos[0], pos[1] + 1
        rs_pos = self.restriction_site.pos
        # if the restriction site is long enough: check that it doesn't prevent
        # the hairpin from forming a valid bp
        if (rs_pos[1] - rs_pos[0] >= rna.HAIRPIN_MIN_SIZE + 2 and
            not self.valid_struct_seq_matching(new_pos, rs_pos, hairpin)):
            return False
        # check intersections with the aptamer in the same state
        return pos[1] < self.aptamers[hairpin.state].pos[0]

    def valid_decrease_loop(self, hairpin):
        pos = hairpin.pos
        new_pos = pos[0], pos[1] - 1
        rs_pos = self.restriction_site.pos
        # if the restriction site is long enough: check that it doesn't prevent
        # the hairpin from forming a valid bp
        if (rs_pos[1] - rs_pos[0] >= rna.HAIRPIN_MIN_SIZE + 2 and
            not self.valid_struct_seq_matching(new_pos, rs_pos, hairpin)):
            return False
        # check that after shifting..
        if hairpin.state == rs_e.State.unbound:
            # ..the bound hairpin and the bound aptamer are still in range and
            # the restriction site is still uncovered
            aptamer_b = self.aptamers[rs_e.State.bound]
            if (aptamer_b.pos[0] - new_pos[1] > UNBOUND_H_APT_DIFF):
                return False
        else:  # hairpin.state == bound
            # ..the unbound hairpin is still in range and the restriction site
            # is still covered
            if (pos[1] <= self.hairpins[rs_e.State.unbound].pos[0] or
                new_pos[1] <= self.restriction_site.pos[0]):
                return False
        return True

    def valid_hairpin_remove_bp(self, hairpin):
        # the hairpin stem size should not drop below its min size
        num_basepairs = len(hairpin.struct.bp_positions)
        new_pos = hairpin.pos[0] + 1, hairpin.pos[1] - 1
        # check that after removing a bp..
        if hairpin.state == rs_e.State.bound:
            # ..the restriction site is still covered and the unbound hairpin
            # is still in range
            hairpin_ub = self.hairpins[rs_e.State.unbound]
            return (num_basepairs > HAIRPIN_STEM_MIN_SIZE and
                    new_pos[0] < self.restriction_site.pos[1] and
                    new_pos[0] < hairpin_ub.pos[0] and
                    new_pos[1] >= hairpin_ub.pos[0])
        else:  # hairpin.state == unbound
            # ..the bound aptamer and the bound hairpin are still in range
            aptamer_b = self.aptamers[rs_e.State.bound]
            return (num_basepairs > HAIRPIN_STEM_MIN_SIZE and
                    self.hairpins[rs_e.State.bound].pos[1] >= new_pos[0] and
                    aptamer_b.pos[0] - hairpin.pos[1] <= UNBOUND_H_APT_DIFF)

    def valid_hairpin_insert_bp(self, hairpin):
        # a bp can be insert if the hairpin can be shifted in both directions
        # (to insert a b on each site) and the hairpin stem max size is not
        # exceeded
        return (len(hairpin.struct.bp_positions) < HAIRPIN_STEM_MAX_SIZE and
                self.valid_hairpin_shift_up(hairpin) and
                self.valid_hairpin_shift_down(hairpin))

    # used to check if either a new h_pos or rs_pos prevents the hairpin from
    # forming a valid bp
    def valid_struct_seq_matching(self, h_pos, rs_pos, hairpin):
        if h_pos[1] <= rs_pos[0] or rs_pos[1] <= h_pos[0]:
            return True
        for pos in xrange(rs_pos[0], sum(rs_pos) / 2 + 1):
            if pos < h_pos[0] or pos >= h_pos[1]:
                continue
            bp_pos = hairpin.struct.basepairs[pos - h_pos[0]]
            if bp_pos is None:
                continue
            bp_pos += h_pos[0]
            if bp_pos < rs_pos[0] or bp_pos > rs_pos[1]:
                continue
            bp_id_i = self.restriction_site.seq[pos - rs_pos[0]]
            bp_id_j = self.restriction_site.seq[bp_pos - rs_pos[0] - 1]
            if rna.valid_bp_b_ids(bp_id_i, bp_id_j):
                continue
            return False
        return True
