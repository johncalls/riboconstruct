import itertools

from .fillings import filling

from .. import rna


class Subsolution(object):
    def __init__(self, substruct, seq_constraint, predecessor=None):
        self.substruct = substruct
        self.predecessor = predecessor
        self.seq_constraint = seq_constraint

        self.collection = Collection()

    @property
    def num_free_bases(self):
        return self.substruct.num_free_bases

    @property
    def bp_b_positions(self):
        return self.substruct.bp_b_positions

    @property
    def is_stem_end(self):
        return self.substruct.is_stem_end

    @property
    def bp_pos_id(self):
        return self.substruct.bp_pos_id

    @property
    def predecessor_pos_ids(self):
        return self.substruct.predecessor_pos_ids

    def b_constrained(self, b_pos, b_id):
        return self.seq_constraint[b_pos][b_id]

    def bp_constrained(self, bp_pos_i, bp_pos_j, bp_id_i, bp_id_j):
        return (self.seq_constraint[bp_pos_i][bp_id_i] or
                self.seq_constraint[bp_pos_j][bp_id_j])

    def iter_fillings(self):
        return self.collection.iter_fillings()

    def iter_bp_id_fillings(self, bp_id):
        return self.collection.iter_bp_id_fillings(bp_id)

    def get_fillings(self):
        return self.collection.get_fillings()

    def append_min_filling(self, filling):
        self.collection.append_min_filling(filling)

    def append(self, filling):
        self.collection.append(filling)

    def delete_predecessor(self):
        # remove predecessor in order to get rid of unnecessary fillings;
        # references to the needed fillings are kept in the current
        # subsolution's fillings
        del self.predecessor
        self.predecessor = None


class Collection(object):
    def __init__(self):
        self.collection = [[] for _ in xrange(rna.BasepairId.count)]

    def __getitem__(self, bp_id):
        return self.collection[bp_id]

    def iter_fillings(self):
        return itertools.chain.from_iterable(
            self.collection[bp_id] for bp_id in xrange(rna.BasepairId.count))

    def iter_bp_id_fillings(self, bp_id):
        return iter(self.collection[bp_id])

    def get_fillings(self):
        return sum(self.collection, [])

    def append_min_filling(self, filling):
        bp_id = filling.bp_id

        try:
            collected_filling = self.collection[bp_id][0]
        except IndexError:
            self.collection[bp_id] = [filling]
        else:
            if filling.energy < collected_filling.energy:
                self.collection[bp_id] = [filling]
            # elif helper.float_equal(filling.energy, collected_filling.energy):
            #     self.append(filling)

    def append(self, filling):
        bp_id = filling.bp_id

        try:
            self.collection[bp_id].append(filling)
        except KeyError:
            self.collection[bp_id] = [filling]


class Stem(object):
    def __init__(self, subsolution, stem_id, predecessor=None):
        self.subsolution = subsolution
        self.stem_id = stem_id
        self.predecessor = predecessor

        self.collection = StemCollection()

    @property
    def bp_pos_id(self):
        return self.subsolution.bp_pos_id

    @property
    def stem_subsolution(self):
        return self.subsolution.predecessor[self.stem_id]

    @property
    def stem_bp_pos_id(self):
        return self.stem_subsolution.bp_pos_id

    @property
    def stem_bp_b_positions(self):
        return self.stem_subsolution.bp_b_positions

    @property
    def num_free_bases_right(self):
        return self.subsolution.num_free_bases[self.stem_id]

    @property
    def seq_constraint(self):
        return self.subsolution.seq_constraint

    def b_constrained(self, b_pos, b_id):
        return self.subsolution.b_constrained(b_pos, b_id)

    def bp_constrained(self, bp_pos_i, bp_pos_j, bp_id_i, bp_id_j):
        return (
            self.subsolution.bp_constrained(
                bp_pos_i, bp_pos_j, bp_id_i, bp_id_j))

    def iter_fillings(self):
        return self.collection.iter_fillings()

    def get_fillings(self):
        return self.collection.get_fillings()

    def iter_bp_id_fillings(self, fixed_bp_id):
        return self.collection.iter_bp_id_fillings(fixed_bp_id)

    def append_min_filling(self, filling):
        self.collection.append_min_filling(filling)

    def append(self, filling):
        self.collection.append(filling)

    def delete_predecessor(self):
        # remove predecessor in order to get rid of unnecessary fillings;
        # references to the needed fillings are kept in the current
        # subsolution's fillings
        del self.predecessor
        self.predecessor = None


class StemCollection(object):
    def __init__(self):
        self.collection = [[[] for _ in xrange(rna.BasepairId.count)]
                           for _ in xrange(rna.BasepairId.count)]

    def iter_bp_id_fillings(self, fixed_bp_id):
        return (
            itertools.chain.from_iterable(
                self.collection[fixed_bp_id][bp_id]
                for bp_id in xrange(rna.BasepairId.count)))

    def append_min_filling(self, filling):
        fixed_bp_id = filling.fixed_bp_id
        bp_id = filling.bp_id

        try:
            collected_filling = self.collection[fixed_bp_id][bp_id][0]
        except IndexError:
            self.collection[fixed_bp_id][bp_id] = [filling]
        else:
            if (filling.energy < collected_filling.energy):
                self.collection[fixed_bp_id][bp_id] = [filling]
            # elif helper.float_equal(filling.energy, collected_filling.energy):
            #     self.append(filling)

    def append(self, filling):
        fixed_bp_id = filling.fixed_bp_id
        bp_id = filling.bp_id

        try:
            self.collection[fixed_bp_id][bp_id].append(filling)
        except KeyError:
            self.collection[fixed_bp_id][bp_id] = [filling]


class HairpinSubsolution(Subsolution):
    def __init__(self, substruct, seq_constraint):
        super(HairpinSubsolution, self).__init__(substruct, seq_constraint)

        self.new_filling = filling.HairpinFilling


class BulgeSubsolution(Subsolution):
    def __init__(self, substruct, seq_constraint, predecessor):
        super(BulgeSubsolution, self).__init__(
            substruct, seq_constraint, predecessor)

        self.new_filling = filling.BulgeFilling


class StackSubsolution(Subsolution):
    def __init__(self, substruct, seq_constraint, predecessor):
        super(StackSubsolution, self).__init__(
            substruct, seq_constraint, predecessor)

        self.new_filling = filling.StackFilling


class InteriorSubsolution(Subsolution):
    def __init__(self, substruct, seq_constraint, predecessor):
        super(InteriorSubsolution, self).__init__(
            substruct, seq_constraint, predecessor)

        self.new_filling = filling.InteriorFilling


class MultiloopSubsolution(Subsolution):
    def __init__(self, substruct, seq_constraint, predecessor):
        super(MultiloopSubsolution, self).__init__(
            substruct, seq_constraint, predecessor)

        self.new_filling = filling.MultiloopFilling
        self.new_multiloop_stem_filling = filling.MultiloopStemFilling
        self.predecessor = (
            sorted(predecessor,
                   key=lambda p: self.predecessor_pos_ids.index(p.bp_pos_id)))

    def get_multiloop_stem(self, stem_id, predecessor=None):
        return MultiloopStem(self, stem_id, predecessor)


class MultiloopStem(Stem):
    def __init__(self, subsolution, stem_id, predecessor=None):
        super(MultiloopStem, self).__init__(subsolution, stem_id, predecessor)

        self.new_filling = filling.MultiloopStemFilling


class FinalSubsolution(Subsolution):
    def __init__(self, substruct, seq_constraint, predecessor):
        super(FinalSubsolution, self).__init__(
            substruct, seq_constraint, predecessor)

        self.new_filling = filling.FinalFilling
        self.new_exterior_stem_filling = filling.ExteriorFilling
        self.predecessor = (
            sorted(predecessor,
                   key=lambda p: self.predecessor_pos_ids.index(p.bp_pos_id)))

    def get_exterior_stem(self, stem_id, predecessor=None):
        return ExteriorStem(self, stem_id, predecessor)


class ExteriorStem(Stem):
    def __init__(self, subsolution, stem_id, predecessor):
        super(ExteriorStem, self).__init__(subsolution, stem_id, predecessor)

        self.new_filling = filling.ExteriorFilling
