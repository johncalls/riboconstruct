from . import filling as fl
from . import subsolution


class TwoTargetMultiloopSubsolution(subsolution.TwoTargetSubsolution):
    def __init__(self, struct_id, substruct, task, predecessor):
        super(TwoTargetMultiloopSubsolution, self).__init__(
            struct_id, substruct, task, predecessor)

        self.predecessor = (
            sorted(predecessor,
                   key=lambda p: self.predecessor_pos_ids.index(p.bp_pos_id)))

        self.new_filling = fl.TwoTargetMultiloopFilling
        self.new_multiloop_stem_filling = fl.TwoTargetMultiloopStemFilling

        bp_pos_i, bp_pos_j = self.bp_b_positions
        num_free_bases_left = self.num_free_bases[-1]
        # dependencies
        if self.is_dependent:
            self.filling_b_index = []
            if self.depends_on(bp_pos_i):
                self.filling_b_index.append(0)
            if self.depends_on(bp_pos_j):
                self.filling_b_index.append(1)
            if num_free_bases_left:
                if self.depends_on(bp_pos_i + 1):
                    self.filling_b_index.append(2)
                if (num_free_bases_left > 1 and
                    self.depends_on(bp_pos_i + num_free_bases_left)):
                    self.filling_b_index.append(3)
        # b_pos to filling base index + b_pos to bp_pos_id
        self.b_pos_mapping[bp_pos_i] = (0, self.bp_pos_id)
        self.b_pos_mapping[bp_pos_j] = (1, self.bp_pos_id)
        if num_free_bases_left:
            self.b_pos_mapping[bp_pos_i + 1] = (2, self.bp_pos_id)
            if num_free_bases_left > 1:
                pos = bp_pos_i + num_free_bases_left
                self.b_pos_mapping[pos] = (3, self.bp_pos_id)

    def get_multiloop_stem(self, stem_id, predecessor=None):
        return TwoTargetMultiloopStem(self, stem_id, predecessor)

    def append_min_filling(self, filling):
        if self.filling_is_constrained(filling):
            return
        if self.is_dependent:
            index_predecessor = filling.predecessor.index
            index = tuple(filling.base_ids[i] for i in self.filling_b_index)
            filling.index = sum(filter(None, (index_predecessor, index)), ())
        else:
            filling.index = filling.predecessor.index

        super(TwoTargetMultiloopSubsolution, self).append_min_filling(
            filling)

    def append(self, filling):
        if self.filling_is_constrained(filling):
            return

        super(TwoTargetMultiloopSubsolution, self).append(filling)

    def filling_is_constrained(self, filling):
        bp_pos_i, bp_pos_j = self.bp_b_positions
        if self.opp_struct_b_constrained(filling, bp_pos_i,
                                         filling.base_ids[0]):
            return True
        if self.opp_struct_b_constrained(filling, bp_pos_j,
                                         filling.base_ids[1]):
            return True
        num_free_bases_left = self.num_free_bases[-1]
        if num_free_bases_left:
            if self.opp_struct_b_constrained(filling, bp_pos_i + 1,
                                             filling.base_ids[2]):
                return True
            if (num_free_bases_left > 1 and
                self.opp_struct_b_constrained(
                    filling, bp_pos_i + num_free_bases_left,
                    filling.base_ids[3])):
                return True
        return False


class TwoTargetMultiloopStem(subsolution.TwoTargetStem):
    def __init__(self, subsolution, stem_id, predecessor=None):
        super(TwoTargetMultiloopStem, self).__init__(
            subsolution, stem_id, predecessor)

        self.new_filling = fl.TwoTargetMultiloopStemFilling

        bp_pos_i, bp_pos_j = self.stem_bp_b_positions
        # dependencies
        if self.is_dependent:
            self.filling_b_index = []
            if self.num_free_bases_right:
                if self.depends_on(bp_pos_j + 1):
                    self.filling_b_index.append(2)
                if (self.stem_id > 0 and
                    self.num_free_bases_right > 1 and
                    self.depends_on(bp_pos_j + self.num_free_bases_right)):
                    self.filling_b_index.append(3)
        # b_pos to filling base index + b_pos to bp_pos_id
        # NOTE: bp_pos_i and bp_pos_j are overwritten with index for ml stem
        # fillings
        self.b_pos_mapping[bp_pos_i] = (0, self.stem_bp_pos_id)
        self.b_pos_mapping[bp_pos_j] = (1, self.stem_bp_pos_id)
        if self.num_free_bases_right:
            self.b_pos_mapping[bp_pos_j + 1] = (2, self.stem_bp_pos_id)
            if self.num_free_bases_right > 1:
                pos = bp_pos_j + self.num_free_bases_right
                self.b_pos_mapping[pos] = (3, self.stem_bp_pos_id)

    def set_index(self, filling):
        if self.is_dependent:
            index = tuple(filling.base_ids[i] for i in self.filling_b_index)
            if self.predecessor:
                index_predecessor = filling.predecessor.index
                index_stem = filling.stem_filling.index
                filling.index = (
                    sum(filter(None, (index_predecessor, index_stem, index)),
                        ()))
            else:
                index_stem = filling.stem_filling.index
                filling.index = sum(filter(None, (index_stem, index)), ())
        else:
            if self.predecessor:
                index_predecessor = filling.predecessor.index
                index_stem = filling.stem_filling.index
                filling.index = (
                    sum(filter(None, (index_predecessor, index_stem)), ()))
            else:
                filling.index = filling.stem_filling.index

    def append_min_filling(self, filling):
        if self.filling_is_constrained(filling):
            return
        self.set_index(filling)

        super(TwoTargetMultiloopStem, self).append_min_filling(filling)

    def append(self, filling):
        if self.filling_is_constrained(filling):
            return
        if not filling.index:
            filling.index = self.set_index(filling)

        super(TwoTargetMultiloopStem, self).append(filling)

    def filling_is_constrained(self, filling):
        bp_pos_j = self.stem_bp_b_positions[1]
        if self.num_free_bases_right:
            if self.opp_struct_b_constrained(filling, bp_pos_j + 1,
                                             filling.base_ids[2]):
                return True
            if (self.num_free_bases_right > 1 and
                self.opp_struct_b_constrained(
                    filling, bp_pos_j + self.num_free_bases_right,
                    filling.base_ids[3])):
                return True
        return False
