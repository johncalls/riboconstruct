import itertools

from .. import rna


class Structure(rna.Structure):
    def __init__(self, struct):
        # omit unnecessary checks
        if isinstance(struct, rna.Structure):
            self.struct = struct.struct
            self._basepairs = struct.basepairs
            self.__bp_positions = struct.bp_positions
        else:
            super(Structure, self).__init__(struct)

        self.predecessor_pos_ids = tuple(self.get_predecessor_pos_ids())

    @property
    def pos_2_bp_pos_id(self):
        try:
            return self.__pos_2_bp_pos_id
        except AttributeError:
            self.__pos_2_bp_pos_id = tuple(self.get_pos_2_bp_pos_id())
            return self.__pos_2_bp_pos_id

    @property
    def successor_pos_id(self):
        try:
            return self.__successor_pos_id
        except AttributeError:
            self.__successor_pos_id = tuple(self.get_successor_pos_id())
            return self.__successor_pos_id

    def get_predecessor_pos_ids(self):
        predecessor_pos_ids = []
        possible_ml_stem_pos = []
        exterior_loop_ids = []
        # first bp closes a hairpin: no predecessor
        predecessor_pos_ids.append((None, ))
        try:
            bp_pos_prev = self.bp_positions[0]
        except IndexError:
            # in the case that there is not one bp
            return tuple(predecessor_pos_ids)
        for bp_pos_id, bp_pos in enumerate(self.bp_positions[1:], 1):
            # check if the bp has a predecessor
            if bp_pos[1] > bp_pos_prev[1]:
                # check if the current bp closes the current ml
                if (possible_ml_stem_pos and
                    (bp_pos[1] >
                        self.bp_positions[possible_ml_stem_pos[-1]][1])):
                    # check if each bp pos id closes in fact a ml stem
                    ml_stem_ids = []
                    for stem_bp_pos_id in possible_ml_stem_pos:
                        # check if the bp pos id closes an exterior loop
                        if bp_pos[1] < self.bp_positions[stem_bp_pos_id][0]:
                            exterior_loop_ids.append(stem_bp_pos_id)
                        # otherwise it is a ml stem
                        else:
                            ml_stem_ids.append(stem_bp_pos_id)
                    # predecessor is part of the current ml
                    ml_stem_ids.append(bp_pos_id - 1)
                    predecessor_pos_ids.append(tuple(ml_stem_ids))
                    possible_ml_stem_pos = []
                else:
                    # add predecessor
                    predecessor_pos_ids.append((bp_pos_id - 1, ))
            # bp has no predecessor (see definition of the bp order)
            else:
                predecessor_pos_ids.append((None, ))
                # the previous bp could be part of the current ml
                possible_ml_stem_pos.append(bp_pos_id - 1)
            bp_pos_prev = bp_pos
        # check if there is more than one ending stem; this information is
        # stored under predecessor_pos_ids[len(self.bp_positions)]
        exterior_loop_ids.extend(possible_ml_stem_pos)
        if len(exterior_loop_ids):
            # also add the last closing stem as exterior loop
            last_bp_pos_id = len(self.bp_positions) - 1
            if last_bp_pos_id not in exterior_loop_ids:
                exterior_loop_ids.append(last_bp_pos_id)
            predecessor_pos_ids.append(tuple(exterior_loop_ids))
        else:
            predecessor_pos_ids.append((len(self.bp_positions) - 1, ))
        return predecessor_pos_ids

    def get_pos_2_bp_pos_id(self):
        pos_2_bp_pos_id = {}
        for bp_pos_id, (bp_pos_i, bp_pos_j) in enumerate(self.bp_positions):
            pos_2_bp_pos_id[bp_pos_i] = bp_pos_id
            pos_2_bp_pos_id[bp_pos_j] = bp_pos_id
        return [pos_2_bp_pos_id[i] if i in pos_2_bp_pos_id else None
                for i in xrange(len(self))]

    def get_successor_pos_id(self):
        successor_pos_id = {}
        for bp_pos_id in xrange(1, len(self.bp_positions) + 1):
            for pred_bp_pos_id in self.predecessor_pos_ids[bp_pos_id]:
                successor_pos_id[pred_bp_pos_id] = bp_pos_id
        successor_pos_id[len(self.bp_positions)] = None
        return [successor_pos_id[i]
                for i in xrange(len(self.bp_positions) + 1)]

    def iter_substructs(self):
        for bp_pos_id in xrange(len(self.bp_positions) + 1):
            if self.predecessor_pos_ids[bp_pos_id][0] is None:
                yield HairpinSubstructure(bp_pos_id, self)
            elif bp_pos_id == len(self.bp_positions):
                yield FinalSubstructure(bp_pos_id, self)
            elif len(self.predecessor_pos_ids[bp_pos_id]) > 1:
                yield MultiloopSubstructure(bp_pos_id, self)
            else:
                pred_bp_pos = (
                    self.bp_positions[self.predecessor_pos_ids[bp_pos_id][0]])
                bp_pos = self.bp_positions[bp_pos_id]
                if pred_bp_pos[0] - bp_pos[0] == 1:
                    if bp_pos[1] - pred_bp_pos[1] == 1:
                        yield StackSubstructure(bp_pos_id, self)
                    else:
                        yield BulgeSubstructure(bp_pos_id, self)
                elif bp_pos[1] - pred_bp_pos[1] == 1:
                    yield BulgeSubstructure(bp_pos_id, self)
                else:
                    yield InteriorSubstructure(bp_pos_id, self)

    def print_substruct(self, substruct):
        i, j = substruct.bp_b_positions
        return ''.join(itertools.islice(self.struct, i, j + 1))

    def is_stem_end(self, bp_pos_id):
        # stem end..
        # -if it is the last bp pos id
        # -if the successor is 'left' of the current bp or if the successor is
        #  closing a ml
        if bp_pos_id == len(self.bp_positions) - 1:
            return True
        succ_bp_pos_id = bp_pos_id + 1
        try:
            succ_bp_pos = self.bp_positions[succ_bp_pos_id]
        except IndexError:
            return False
        else:
            bp_pos_i, _ = self.bp_positions[bp_pos_id]
            return (len(self.predecessor_pos_ids[succ_bp_pos_id]) > 1 or
                    succ_bp_pos[1] < bp_pos_i)


class Substructure(object):
    def __init__(self, bp_pos_id, struct):
        self.bp_pos_id = bp_pos_id
        self.struct_type = None

        self.predecessor_pos_ids = struct.predecessor_pos_ids[self.bp_pos_id]
        try:
            self.bp_b_positions = struct.bp_positions[self.bp_pos_id]
        except IndexError:
            pass
        self.is_stem_end = struct.is_stem_end(bp_pos_id)

    def __len__(self):
        i, j = self.bp_b_positions
        return j - i + 1


class HairpinSubstructure(Substructure):
    def __init__(self, bp_pos_id, struct):
        super(HairpinSubstructure, self).__init__(bp_pos_id, struct)
        self.struct_type = rna.StructType.HAIRPIN
        self.num_free_bases = self.calc_num_free_bases()

    def calc_num_free_bases(self):
        return self.bp_b_positions[1] - self.bp_b_positions[0] - 1


class StackSubstructure(Substructure):
    def __init__(self, bp_pos_id, struct):
        super(StackSubstructure, self).__init__(bp_pos_id, struct)
        self.struct_type = rna.StructType.STACKING
        self.num_free_bases = 0


class BulgeSubstructure(Substructure):
    def __init__(self, bp_pos_id, struct):
        super(BulgeSubstructure, self).__init__(bp_pos_id, struct)
        self.struct_type = rna.StructType.BULGE
        self.num_free_bases = self.calc_num_free_bases(struct)

    def calc_num_free_bases(self, struct):
        pred_pos = struct.bp_positions[self.predecessor_pos_ids[0]]
        return (pred_pos[0] - self.bp_b_positions[0] - 1,
                self.bp_b_positions[1] - pred_pos[1] - 1)


class InteriorSubstructure(Substructure):
    def __init__(self, bp_pos_id, struct):
        super(InteriorSubstructure, self).__init__(bp_pos_id, struct)
        self.struct_type = rna.StructType.INTERIOR
        self.num_free_bases = self.calc_num_free_bases(struct)

    def calc_num_free_bases(self, struct):
        pred_pos = struct.bp_positions[self.predecessor_pos_ids[0]]
        return (pred_pos[0] - self.bp_b_positions[0] - 1,
                self.bp_b_positions[1] - pred_pos[1] - 1)


class MultiloopSubstructure(Substructure):
    def __init__(self, bp_pos_id, struct):
        super(MultiloopSubstructure, self).__init__(bp_pos_id, struct)
        self.struct_type = rna.StructType.MULTILOOP
        self.num_free_bases = self.calc_num_free_bases(struct)

    def calc_num_free_bases(self, struct):
        num_free_bases = []
        prev_stem_pos = struct.bp_positions[
            self.predecessor_pos_ids[0]]
        num_free_bases.append(self.bp_b_positions[1] - prev_stem_pos[1] - 1)
        for stem_pos_id in self.predecessor_pos_ids[1:]:
            stem_pos = struct.bp_positions[stem_pos_id]
            num_free_bases.append(prev_stem_pos[0] - stem_pos[1] - 1)
            prev_stem_pos = stem_pos
        num_free_bases.append(prev_stem_pos[0] - self.bp_b_positions[0] - 1)
        return tuple(num_free_bases)


class FinalSubstructure(Substructure):
    def __init__(self, bp_pos_id, struct):
        super(FinalSubstructure, self).__init__(bp_pos_id, struct)
        self.struct_type = rna.StructType.EXTERIOR
        self.bp_b_positions = 0, len(struct) - 1
        self.num_free_bases = self.calc_num_free_bases(struct)

    def calc_num_free_bases(self, struct):
        num_free_bases = []
        prev_stem_pos = struct.bp_positions[
            self.predecessor_pos_ids[0]]
        num_free_bases.append(len(struct) - prev_stem_pos[1] - 1)
        # if there is only one stem the for-loop will not be executed
        for stem_pos_id in self.predecessor_pos_ids[1:]:
            stem_pos = struct.bp_positions[stem_pos_id]
            num_free_bases.append(prev_stem_pos[0] - stem_pos[1] - 1)
            prev_stem_pos = stem_pos
        num_free_bases.append(prev_stem_pos[0])
        return tuple(num_free_bases)
