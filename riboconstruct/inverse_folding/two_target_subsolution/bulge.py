from . import filling as fl
from . import subsolution


class TwoTargetBulgeSubsolution(subsolution.TwoTargetSubsolution):
    def __init__(self, struct_id, substruct, task, predecessor):
        super(TwoTargetBulgeSubsolution, self).__init__(
            struct_id, substruct, task, predecessor)

        self.new_filling = fl.TwoTargetBulgeFilling

        bp_pos_i, bp_pos_j = self.bp_b_positions
        # dependencies
        if self.is_dependent:
            self.filling_b_index = []
            if self.depends_on(bp_pos_i):
                self.filling_b_index.append(0)
            if self.depends_on(bp_pos_j):
                self.filling_b_index.append(1)
        # b_pos to filling base index + b_pos to bp_pos_id
        self.b_pos_mapping[bp_pos_i] = (0, self.bp_pos_id)
        self.b_pos_mapping[bp_pos_j] = (1, self.bp_pos_id)

    def append_min_filling(self, filling):
        if self.filling_is_constrained(filling):
            return
        if self.is_dependent:
            index_predecessor = filling.predecessor.index
            index = tuple(filling.base_ids[i] for i in self.filling_b_index)
            filling.index = sum(filter(None, (index_predecessor, index)), ())
        else:
            filling.index = filling.predecessor.index

        super(TwoTargetBulgeSubsolution, self).append_min_filling(filling)

    def append(self, filling):
        if self.filling_is_constrained(filling):
            return

        super(TwoTargetBulgeSubsolution, self).append(filling)

    def filling_is_constrained(self, filling):
        bp_pos_i, bp_pos_j = self.bp_b_positions
        if self.opp_struct_b_constrained(filling, bp_pos_i,
                                         filling.base_ids[0]):
            return True
        if self.opp_struct_b_constrained(filling, bp_pos_j,
                                         filling.base_ids[1]):
            return True
        return False
