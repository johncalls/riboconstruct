from . import filling as fl
from . import subsolution


class TwoTargetHairpinSubsolution(subsolution.TwoTargetSubsolution):
    def __init__(self, struct_id, substruct, task):
        super(TwoTargetHairpinSubsolution, self).__init__(
            struct_id, substruct, task)

        self.new_filling = fl.TwoTargetHairpinFilling

        bp_pos_i, bp_pos_j = self.bp_b_positions
        # dependencies
        if self.is_dependent:
            self.filling_b_index = []
            if self.depends_on(bp_pos_i):
                self.filling_b_index.append(0)
            if self.num_free_bases == 3:
                if self.depends_on(bp_pos_j):
                    self.filling_b_index.append(4)
            #     if self.depends_on(bp_pos_i + 1):
            #         self.filling_b_index.append((1, bp_pos_i + 1))
            #     if self.depends_on(bp_pos_i + 2):
            #         self.filling_b_index.append((2, bp_pos_i + 2))
            #     if self.depends_on(bp_pos_i + 3):
            #         self.filling_b_index.append((3, bp_pos_i + 3))
            if self.num_free_bases == 4:
                if self.depends_on(bp_pos_j):
                    self.filling_b_index.append(5)
                if self.depends_on(bp_pos_i + 1):
                    self.filling_b_index.append(1)
                if self.depends_on(bp_pos_i + 2):
                    self.filling_b_index.append(2)
                if self.depends_on(bp_pos_i + 3):
                    self.filling_b_index.append(3)
                if self.depends_on(bp_pos_i + 4):
                    self.filling_b_index.append(4)
            elif self.num_free_bases > 4:
                if self.depends_on(bp_pos_j):
                    self.filling_b_index.append(3)
                if self.depends_on(bp_pos_i + 1):
                    self.filling_b_index.append(1)
                if self.depends_on(bp_pos_j - 1):
                    self.filling_b_index.append(2)
        # b_pos to filling base index + b_pos to bp_pos_id
        if self.num_free_bases == 3:
            self.b_pos_mapping[bp_pos_i] = (0, self.bp_pos_id)
            self.b_pos_mapping[bp_pos_j] = (4, self.bp_pos_id)
        elif self.num_free_bases == 4:
            for i, b_pos in enumerate(xrange(bp_pos_i, bp_pos_j + 1)):
                self.b_pos_mapping[b_pos] = (i, self.bp_pos_id)
        else:
            self.b_pos_mapping[bp_pos_i] = (0, self.bp_pos_id)
            self.b_pos_mapping[bp_pos_i + 1] = (1, self.bp_pos_id)
            self.b_pos_mapping[bp_pos_j - 1] = (2, self.bp_pos_id)
            self.b_pos_mapping[bp_pos_j] = (3, self.bp_pos_id)

    def append_min_filling(self, filling):
        if self.filling_is_constrained(filling):
            return
        if self.is_dependent:
            filling.index = (
                tuple(filling.base_ids[i] for i in self.filling_b_index))

        super(TwoTargetHairpinSubsolution, self).append_min_filling(
            filling)

    def append(self, filling):
        if self.filling_is_constrained(filling):
            return

        super(TwoTargetHairpinSubsolution, self).append(filling)

    def filling_is_constrained(self, filling):
        bp_pos_i, bp_pos_j = self.bp_b_positions
        if self.opp_struct_b_constrained(filling, bp_pos_i,
                                         filling.base_ids[0]):
            return True
        if self.opp_struct_b_constrained(filling, bp_pos_j,
                                         filling.base_ids[-1]):
            return True
        if (self.num_free_bases == 4 and
            (self.opp_struct_b_constrained(filling, bp_pos_i + 1,
                                           filling.base_ids[1]) or
             self.opp_struct_b_constrained(filling, bp_pos_i + 2,
                                           filling.base_ids[2]) or
             self.opp_struct_b_constrained(filling, bp_pos_i + 3,
                                           filling.base_ids[3]) or
             self.opp_struct_b_constrained(filling, bp_pos_i + 4,
                                           filling.base_ids[4]))):
            return True
        elif (self.num_free_bases > 4 and
              (self.opp_struct_b_constrained(filling, bp_pos_i + 1,
                                             filling.base_ids[1]) or
               self.opp_struct_b_constrained(filling, bp_pos_j - 1,
                                             filling.base_ids[2]))):
            return True
        return False
