from . import riboswitch as rs
from . import element as rs_e

from .. import rna


HAIRPIN_LOOP_MAX_SIZE = 7
HAIRPIN_LOOP_MIN_SIZE = 3
HAIRPIN_STEM_MIN_SIZE = 5
HAIRPIN_STEM_MAX_SIZE = 30

UNBOUND_H_APT_DIFF = 2


class SimpleSiblingsGenerator(object):
    def __init__(self, parent):
        self.parent = parent

        rs_es = self.parent.elements
        self.hairpins = rs_es[rs_e.Type.hairpin]
        if self.hairpins[0].state != 0:
            self.hairpins = [self.hairpins[1], self.hairpins[0]]
        self.aptamers = rs_es[rs_e.Type.aptamer]
        if self.aptamers[0].state != 0:
            self.aptamers = [self.aptamers[1], self.aptamers[0]]
        self.restriction_site = rs_es[rs_e.Type.restriction_site][0]
        try:
            self.context_front = rs_es[rs_e.Type.context_front][0]
        except IndexError:
            self.context_front = None

    def iter_siblings(self):
        # alter the riboswitch elements
        for rs_elem, rs_elem_new in self.iter_altered_riboswitch_elements():
            sibling = self.parent.copy()
            sibling.replace(rs_elem, rs_elem_new)
            yield sibling

        # alter the functional site size
        if not len(self.parent.elements[rs_e.Type.context_front]):
            # increase
            if self.valid_fs_shift_down():
                sibling = self.parent.copy()
                sibling.pos[0] -= 1
                sibling.pos_functional_site[0] = sibling.pos[0]
                sibling.pos_riboswitch[0] = sibling.pos[0]
                yield sibling

            # decrease
            if self.valid_fs_shift_up():
                sibling = self.parent.copy()
                sibling.pos[0] += 1
                sibling.pos_functional_site[0] = sibling.pos[0]
                sibling.pos_riboswitch[0] = sibling.pos[0]
                yield sibling

    def iter_altered_riboswitch_elements(self):
        # alter hairpins
        for hairpin in self.hairpins:
            h_shift_up = self.valid_hairpin_shift_up(hairpin)
            if h_shift_up:
                yield hairpin, hairpin.shift_up()
            h_shift_down = self.valid_hairpin_shift_down(hairpin)
            if h_shift_down:
                yield hairpin, hairpin.shift_down()
            # self.valid_hairpin_insert_bp()
            if (h_shift_up and
                h_shift_down and
                len(hairpin.struct.bp_positions) < HAIRPIN_STEM_MAX_SIZE):
                yield hairpin, hairpin.insert_bp_before(0)
            if self.valid_hairpin_remove_bp(hairpin):
                yield hairpin, hairpin.remove_bp(0)
            # # alter hairpin loop
            # bp_pos_i, bp_pos_j = hairpin.struct.bp_positions[0]
            # loop_size = bp_pos_j - bp_pos_i - 1
            # if loop_size < HAIRPIN_LOOP_MAX_SIZE:
            #     yield hairpin, hairpin.insert_b(bp_pos_i + 1)
            # if loop_size > HAIRPIN_LOOP_MIN_SIZE:
            #     yield hairpin, hairpin.remove_b(bp_pos_i + 1)

        # alter restriction site
        if self.valid_rs_shift_up():
            yield self.restriction_site, self.restriction_site.shift_up()
        if self.valid_rs_shift_down():
            yield self.restriction_site, self.restriction_site.shift_down()

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
        # check intersections with the context front (if available) or the size
        # of the functional site
        if self.context_front is None:
            return self.parent.pos_functional_site[0] < pos[0]
        else:
            return self.context_front.pos[1] < pos[0]

    def valid_hairpin_insert_bp(self, hairpin):
        # a bp can be insert if the hairpin can be shifted in both directions
        # (to insert a b on each site) and the hairpin stem max size is not
        # exceeded
        return (len(hairpin.struct.bp_positions) < HAIRPIN_STEM_MAX_SIZE and
                self.valid_hairpin_shift_up(hairpin) and
                self.valid_hairpin_shift_down(hairpin))

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

    def valid_rs_shift_up(self):
        # after shifting the restriction site should still not interfere with
        # the unbound hairpin and still be covered by the bound hairpin
        pos = self.restriction_site.pos
        if (pos[1] < self.hairpins[rs_e.State.unbound].pos[0] and
            pos[1] < self.hairpins[rs_e.State.bound].pos[1]):
            new_pos = pos[0] + 1, pos[1] + 1
            hairpin_0, hairpin_1 = self.hairpins
            # after shifting the restriction site should not prevent the
            # hairpins from forming a valid bp
            return (self.valid_struct_seq_matching(hairpin_0.pos, new_pos,
                                                   hairpin_0) and
                    self.valid_struct_seq_matching(hairpin_1.pos, new_pos,
                                                   hairpin_1))
        return False

    def valid_rs_shift_down(self):
        pos = self.restriction_site.pos
        new_pos = pos[0] - 1, pos[1] - 1
        hairpin_0, hairpin_1 = self.hairpins
        # after shifting the restriction site should not prevent the hairpins
        # from forming a valid bp
        if (not self.valid_struct_seq_matching(hairpin_0.pos, new_pos,
                                               hairpin_0) or
            not self.valid_struct_seq_matching(hairpin_1.pos, new_pos,
                                               hairpin_1)):
            return False
        # after shifting the restriction site should still be covered by the
        # bound hairpin and either not interfere with the context front or
        # increase the functional site to an invalid size
        if self.context_front is None:
            return (
                pos[0] > self.parent.pos_functional_site[0] and
                new_pos[1] > self.hairpins[rs_e.State.bound].pos[0])
        else:
            return (
                pos[0] > self.context_front.pos[1] and
                new_pos[1] > self.hairpins[rs_e.State.bound].pos[0])

    def valid_cf_shift_up(self):
        pos_end = self.context_front.pos[1]
        if pos_end >= self.restriction_site.pos[0]:
            return False
        for hairpin in self.hairpins:
            if pos_end >= hairpin.pos[0]:
                return False
        # TODO: implement check for seq constraint
        return True

    def valid_cf_shift_down(self):
        pos = self.parent.pos_functional_site
        return pos[1] - pos[0] < rs.FUNCTIONAL_SITE_MAX_LENGTH

    def valid_fs_shift_up(self):
        pos_start = self.parent.pos_functional_site[0]
        if pos_start >= self.restriction_site.pos[0]:
            return False
        for hairpin in self.hairpins:
            if pos_start >= hairpin.pos[0]:
                return False
        # TODO: implement check for seq constraint
        return True

    def valid_fs_shift_down(self):
        pos = self.parent.pos_functional_site
        return pos[1] - pos[0] < rs.FUNCTIONAL_SITE_MAX_LENGTH

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
