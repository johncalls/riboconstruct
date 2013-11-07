from . import filling as fl
from . import subsolution


class TwoTargetFinalSubsolution(subsolution.TwoTargetSubsolution):
    def __init__(self, struct_id, substruct, task, predecessor):
        super(TwoTargetFinalSubsolution, self).__init__(
            struct_id, substruct, task, predecessor)

        self.predecessor = (
            sorted(predecessor,
                   key=lambda p: self.predecessor_pos_ids.index(p.bp_pos_id)))

        self.new_filling = fl.TwoTargetFinalFilling
        self.new_exterior_stem_filling = fl.TwoTargetExteriorFilling

        bp_pos_i, bp_pos_j = self.stem_bp_b_positions
        num_free_bases_left = self.num_free_bases[-1]
        # b_pos to filling base index + b_pos to bp_pos_id
        self.b_pos_mapping[bp_pos_i] = (0, self.bp_pos_id)
        self.b_pos_mapping[bp_pos_j] = (1, self.bp_pos_id)
        if num_free_bases_left:
            self.b_pos_mapping[bp_pos_i - 1] = (2, self.bp_pos_id)

    @property
    def stem_bp_b_positions(self):
        return self.predecessor[-1].bp_b_positions

    def get_exterior_stem(self, stem_id, predecessor=None):
        return TwoTargetExteriorStem(self, stem_id, predecessor)

    def append_min_filling(self, filling):
        if self.filling_is_constrained(filling):
            return
        if self.is_dependent and self.num_free_bases[-1]:
            index_predecessor = filling.predecessor.index
            index = (filling.base_ids[2],)
            filling.index = sum(filter(None, (index_predecessor, index)), ())
        else:
            filling.index = filling.predecessor.index

        super(TwoTargetFinalSubsolution, self).append_min_filling(filling)

    def append(self, filling):
        if self.filling_is_constrained(filling):
            return

        super(TwoTargetFinalSubsolution, self).append(filling)

    def filling_is_constrained(self, filling):
        bp_pos_i = self.bp_b_positions[0]
        num_free_bases_left = self.num_free_bases[-1]
        if (num_free_bases_left and
            self.opp_struct_b_constrained(filling, bp_pos_i - 1,
                                          filling.base_ids[2])):
            return True
        return False


class TwoTargetExteriorStem(subsolution.TwoTargetStem):
    def __init__(self, subsolution, stem_id, predecessor=None):
        super(TwoTargetExteriorStem, self).__init__(
            subsolution, stem_id, predecessor)

        self.new_filling = fl.TwoTargetExteriorFilling

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
            if self.stem_id > 0 and self.num_free_bases_right > 1:
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

        super(TwoTargetExteriorStem, self).append_min_filling(filling)

    def append(self, filling):
        if self.filling_is_constrained(filling):
            return
        if not filling.index:
            self.set_index(filling)

        super(TwoTargetExteriorStem, self).append(filling)

    def filling_is_constrained(self, filling):
        bp_pos_j = self.stem_bp_b_positions[1]
        if self.num_free_bases_right:
            if self.opp_struct_b_constrained(filling, bp_pos_j + 1,
                                             filling.base_ids[2]):
                return True
            if (self.stem_id > 0 and
                self.num_free_bases_right > 1 and
                self.opp_struct_b_constrained(
                    filling, bp_pos_j + self.num_free_bases_right,
                    filling.base_ids[3])):
                return True
        return False
