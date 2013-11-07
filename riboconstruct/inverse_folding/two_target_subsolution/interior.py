from . import filling as fl
from . import subsolution


class TwoTargetInteriorSubsolution(subsolution.TwoTargetSubsolution):
    def __init__(self, struct_id, substruct, task, predecessor):
        super(TwoTargetInteriorSubsolution, self).__init__(
            struct_id, substruct, task, predecessor)

        self.new_filling = fl.TwoTargetInteriorFilling

        bp_pos_i, bp_pos_j = self.bp_b_positions
        num_free_bases_left, num_free_bases_right = self.num_free_bases
        # dependencies
        if self.is_dependent:
            self.filling_b_index = []
            if self.depends_on(bp_pos_i):
                self.filling_b_index.append(0)
            if self.depends_on(bp_pos_j):
                self.filling_b_index.append(1)
            if self.depends_on(bp_pos_i + 1):
                self.filling_b_index.append(2)
            if (num_free_bases_left > 1 and
                self.depends_on(bp_pos_i + num_free_bases_left)):
                self.filling_b_index.append(3)
            if self.depends_on(bp_pos_j - 1):
                self.filling_b_index.append(len(self.filling_b_index))
            if (num_free_bases_right > 1 and
                self.depends_on(bp_pos_j - num_free_bases_right)):
                self.filling_b_index.append(len(self.filling_b_index))
        # b_pos to filling base index + b_pos to bp_pos_id
        self.b_pos_mapping[bp_pos_i] = (0, self.bp_pos_id)
        self.b_pos_mapping[bp_pos_j] = (1, self.bp_pos_id)
        i = 0
        self.b_pos_mapping[bp_pos_i + 1] = (2, self.bp_pos_id)
        if num_free_bases_left > 1:
            pos = bp_pos_i + num_free_bases_left
            self.b_pos_mapping[pos] = (3, self.bp_pos_id)
            i += 1
        pos = bp_pos_j - num_free_bases_right
        self.b_pos_mapping[pos] = (3 + i, self.bp_pos_id)
        if num_free_bases_right > 1:
            pos = bp_pos_j - 1
            self.b_pos_mapping[pos] = (4 + i, self.bp_pos_id)

    def append_min_filling(self, filling):
        if self.filling_is_constrained(filling):
            return
        if self.is_dependent:
            index_predecessor = filling.predecessor.index
            index = tuple(filling.base_ids[i] for i in self.filling_b_index)
            filling.index = sum(filter(None, (index_predecessor, index)), ())
        else:
            filling.index = filling.predecessor.index

        super(TwoTargetInteriorSubsolution, self).append_min_filling(filling)

    def append(self, filling):
        if self.filling_is_constrained(filling):
            return

        super(TwoTargetInteriorSubsolution, self).append(filling)

    def filling_is_constrained(self, filling):
        bp_pos_i, bp_pos_j = self.bp_b_positions
        if self.opp_struct_b_constrained(filling, bp_pos_i,
                                         filling.base_ids[0]):
            return True
        if self.opp_struct_b_constrained(filling, bp_pos_j,
                                         filling.base_ids[1]):
            return True
        num_free_bases_left, num_free_bases_right = self.num_free_bases
        # size_i >= 1 always given
        if self.opp_struct_b_constrained(filling, bp_pos_i + 1,
                                         filling.base_ids[2]):
            return True
        if (num_free_bases_left > 1 and
            self.opp_struct_b_constrained(filling,
                                          bp_pos_i + num_free_bases_left,
                                          filling.base_ids[3])):
            return True
        # size_j >= 1 always given
        if self.opp_struct_b_constrained(filling, bp_pos_j - 1,
                                         filling.base_ids[-1]):
            return True
        if (num_free_bases_right > 1 and
            self.opp_struct_b_constrained(filling,
                                          bp_pos_j - num_free_bases_right,
                                          filling.base_ids[-2])):
            return True
        return False
