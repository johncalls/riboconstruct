from ... import rna


def get_min_list(fillings, new_energy):
    append = False
    try:
        if new_energy < fillings[0].energy:
            # new_energy is smaller than any energy in the list: return new one
            del fillings
            fillings = []
            append = True
        # elif helper.float_equal(new_energy, fillings[0].energy):
        #     # new_energy is equal the other ones: simply return it
        #     append = True
    except IndexError:
        # list is empty, simply append
        append = True
    return append, fillings


class Filling(object):
    def __init__(self, bp_pos_id, bp_id, energy, predecessor=None):
        self.bp_pos_id = bp_pos_id
        self.bp_id = bp_id
        self.energy = energy
        self.predecessor = predecessor

        self.base_ids = []

    @property
    def bp_b_ids(self):
        return rna.BASEPAIRS[self.bp_id]


class HairpinFilling(Filling):
    def __init__(self, bp_pos_id, bp_id, b_ids, size, energy):
        super(HairpinFilling, self).__init__(bp_pos_id, bp_id, energy)

        bp_id_i, bp_id_j = self.bp_b_ids
        self.base_ids.append(bp_id_i)
        self.base_ids.extend(b_ids)
        self.base_ids.append(bp_id_j)
        self.size = size

    def __repr__(self):
        if self.size > 4:
            return ''.join((rna.BASES[self.base_ids[0]],
                            rna.BASES[self.base_ids[1]],
                            rna.BASES[rna.BaseId.UNSPEC] * (self.size - 2),
                            rna.BASES[self.base_ids[2]],
                            rna.BASES[self.base_ids[3]]))
        else:  # self.size == 3 or self.size == 4
            return ''.join((rna.BASES[self.base_ids[0]],
                            rna.BASES[rna.BaseId.UNSPEC] * self.size,
                            rna.BASES[self.base_ids[-1]]))

    def get_seq_list(self, seq, struct):
        i, j = struct.bp_positions[self.bp_pos_id]
        seq[i] = rna.BASES[self.base_ids[0]]
        seq[j] = rna.BASES[self.base_ids[-1]]
        if self.size > 4:
            seq[i + 1] = rna.BASES[self.base_ids[1]]
            seq[j - 1] = rna.BASES[self.base_ids[2]]
        elif self.size == 4:
            seq[i + 1] = rna.BASES[self.base_ids[1]]
            seq[i + 2] = rna.BASES[self.base_ids[2]]
            seq[i + 3] = rna.BASES[self.base_ids[3]]
            seq[i + 4] = rna.BASES[self.base_ids[4]]
        else:  # self.size == 3
            pass
            # unspec = rna.BASES[rna.BaseId.UNSPEC]
            # b_id = self.base_ids[1]
            # seq[i + 1] = rna.BASES[b_id] if b_id != -1 else unspec
            # b_id = self.base_ids[2]
            # seq[i + 2] = rna.BASES[b_id] if b_id != -1 else unspec
            # b_id = self.base_ids[3]
            # seq[i + 3] = rna.BASES[b_id] if b_id != -1 else unspec


class BulgeFilling(Filling):
    def __init__(self, bp_pos_id, bp_id, size, energy, pred):
        super(BulgeFilling, self).__init__(bp_pos_id, bp_id, energy, pred)

        bp_id_i, bp_id_j = self.bp_b_ids
        self.base_ids.append(bp_id_i)
        self.base_ids.append(bp_id_j)
        self.size = size

    def __repr__(self):
        # left bulge..
        if self.size[0]:
            bulge_seq = rna.BASES[rna.BaseId.UNSPEC] * self.size[0]
            return ''.join((rna.BASES[self.base_ids[0]],
                            bulge_seq,
                            repr(self.predecessor),
                            rna.BASES[self.base_ids[1]]))
        # right bulge..
        else:
            bulge_seq = rna.BASES[rna.BaseId.UNSPEC] * self.size[1]
            return ''.join((rna.BASES[self.base_ids[0]],
                            repr(self.predecessor),
                            bulge_seq,
                            rna.BASES[self.base_ids[1]]))

    def get_seq_list(self, seq, struct):
        i, j = struct.bp_positions[self.bp_pos_id]
        seq[i] = rna.BASES[self.base_ids[0]]
        seq[j] = rna.BASES[self.base_ids[1]]
        self.predecessor.get_seq_list(seq, struct)


class StackFilling(Filling):
    def __init__(self, bp_pos_id, bp_id, energy, pred):
        super(StackFilling, self).__init__(bp_pos_id, bp_id, energy, pred)

        bp_id_i, bp_id_j = self.bp_b_ids
        self.base_ids.append(bp_id_i)
        self.base_ids.append(bp_id_j)

    def __repr__(self):
        return ''.join((rna.BASES[self.base_ids[0]],
                        repr(self.predecessor),
                        rna.BASES[self.base_ids[1]]))

    def get_seq_list(self, seq, struct):
        i, j = struct.bp_positions[self.bp_pos_id]
        seq[i] = rna.BASES[self.base_ids[0]]
        seq[j] = rna.BASES[self.base_ids[1]]
        self.predecessor.get_seq_list(seq, struct)


class InteriorFilling(Filling):
    def __init__(self, bp_pos_id, bp_id, b_ids, size, energy, pred):
        super(InteriorFilling, self).__init__(bp_pos_id, bp_id, energy, pred)

        bp_id_i, bp_id_j = self.bp_b_ids
        self.base_ids.append(bp_id_i)
        self.base_ids.append(bp_id_j)
        self.base_ids.extend(b_ids)
        self.size = size

    def __repr__(self):
        size_left, size_right = self.size
        if size_left > 1:
            seq_left = ''.join([rna.BASES[self.base_ids[2]],
                                rna.BASES[rna.BaseId.UNSPEC] * (size_left - 2),
                                rna.BASES[self.base_ids[3]]])
        else:
            seq_left = rna.BASES[self.base_ids[2]]
        if size_right > 1:
            seq_right = ''.join((
                rna.BASES[self.base_ids[-2]],
                rna.BASES[rna.BaseId.UNSPEC] * (size_right - 2),
                rna.BASES[self.base_ids[-1]]))
        else:
            seq_right = rna.BASES[self.base_ids[-1]]
        return ''.join((rna.BASES[self.base_ids[0]],
                        seq_left, repr(self.predecessor), seq_right,
                        rna.BASES[self.base_ids[1]]))

    def get_seq_list(self, seq, struct):
        i, j = struct.bp_positions[self.bp_pos_id]
        seq[i] = rna.BASES[self.base_ids[0]]
        seq[j] = rna.BASES[self.base_ids[1]]
        size_left, size_right = self.size
        seq[i + 1] = rna.BASES[self.base_ids[2]]
        if size_left > 1:
            seq[i + size_left] = rna.BASES[self.base_ids[3]]
        seq[j - 1] = rna.BASES[self.base_ids[-1]]
        if size_right > 1:
            seq[j - size_right] = rna.BASES[self.base_ids[-2]]
        self.predecessor.get_seq_list(seq, struct)


class MultiloopFilling(Filling):
    """
    Represents the filling of a multiloop, i.e. recursively the filling of its
    free bases and its stems as well as the filling of the ml closing bp and
    the adjacent left free bases.

    ML:

                           (i',..,j')->..->(..)   (ml stem fillings)
                       i'-1                   .
                           .                .
        (free_b_ids_left)   .             .
        (num_free_b_left)    .          .
                              i+1     .
                                 (i j)       (bp_id)
                           ..,i-1     j+1,..

    """
    def __init__(self, bp_pos_id, bp_id, free_b_ids_left, num_free_b_left,
                 energy, pred):
        super(MultiloopFilling, self).__init__(bp_pos_id, bp_id, energy, pred)

        bp_id_i, bp_id_j = self.bp_b_ids
        self.base_ids.append(bp_id_i)
        self.base_ids.append(bp_id_j)
        self.base_ids.extend(free_b_ids_left)
        self.num_free_b_left = num_free_b_left

    def __repr__(self):
        if self.num_free_b_left:
            if self.num_free_b_left > 1:
                free_bases = ''.join((
                    rna.BASES[self.base_ids[2]],
                    rna.BASES[rna.BaseId.UNSPEC] * (self.num_free_b_left - 2),
                    rna.BASES[self.base_ids[3]]))
            else:
                free_bases = rna.BASES[self.base_ids[2]]
        else:
            free_bases = ''

        return ''.join((rna.BASES[self.base_ids[0]],
                        free_bases,
                        repr(self.predecessor),
                        rna.BASES[self.base_ids[1]]))

    def get_seq_list(self, seq, struct):
        i, j = struct.bp_positions[self.bp_pos_id]
        seq[i] = rna.BASES[self.base_ids[0]]
        seq[j] = rna.BASES[self.base_ids[1]]
        if self.num_free_b_left:
            seq[i + 1] = rna.BASES[self.base_ids[2]]
            if self.num_free_b_left > 1:
                seq[i + self.num_free_b_left] = rna.BASES[self.base_ids[3]]
        self.predecessor.get_seq_list(seq, struct)


class MultiloopStemFilling(Filling):
    """
    Represents the filling of a stem in a ml, i.e. recursively the filling of
    the stem itself and the filling of the stem on the right as well as the
    filling of the stem's adjacent free bases on the right side.

    ML stem:

                       (stem_filling)      (pred)
                         .       .       .         .
                          .     .         .       .
        (stem_bp_id)       (i j)           (i' j')
                     ..,i-1     j+1,..,i'-1       j'+1

                             (free_b_ids_right)
                             (num_free_b_right)

    (fixed_bp_id) ^= bp id assigned to the ml closing bp
    """
    def __init__(self, fixed_bp_id, bp_pos_id, free_b_ids_right,
                 num_free_b_right, energy, stem_filling, pred=None):
        super(MultiloopStemFilling, self).__init__(
            bp_pos_id, stem_filling.bp_id, energy, pred)

        self.fixed_bp_id = fixed_bp_id
        bp_id_i, bp_id_j = stem_filling.bp_b_ids
        self.base_ids.append(bp_id_i)
        self.base_ids.append(bp_id_j)
        self.base_ids.extend(free_b_ids_right)
        self.num_free_b_right = num_free_b_right
        self.stem_filling = stem_filling

    def __repr__(self):
        if self.predecessor:
            pred = repr(self.predecessor)
        else:
            pred = ''

        if self.num_free_b_right:
            if len(self.base_ids) > 3:
                free_bases = ''.join((
                    rna.BASES[self.base_ids[2]],
                    rna.BASES[rna.BaseId.UNSPEC] * (self.num_free_b_right - 2),
                    rna.BASES[self.base_ids[3]]))
            else:
                free_bases = rna.BASES[self.base_ids[2]]
        else:
            free_bases = ''

        return ''.join([repr(self.stem_filling), free_bases, pred])

    def get_seq_list(self, seq, struct):
        i, j = struct.bp_positions[self.stem_filling.bp_pos_id]
        if self.num_free_b_right:
            seq[j + 1] = rna.BASES[self.base_ids[2]]
            if self.num_free_b_right > 1:
                seq[j + self.num_free_b_right] = rna.BASES[self.base_ids[3]]
        self.stem_filling.get_seq_list(seq, struct)
        if self.predecessor:
            self.predecessor.get_seq_list(seq, struct)


class ExteriorFilling(Filling):
    """
    Exterior loop:

                               (stem_filling)      (pred)
                                 .       .       .         .
                                  .     .         .       .
                   (bp_id)         (i j)           (i' j')
        (free_b_ids_right)   ..,i-1     j+1,..,i'-1       j'+1

                                     (num_free_b_right)

    (fixed_bp_id) ^= bp id assigned to the first exterior stem
    """
    def __init__(self, fixed_bp_id, bp_pos_id, free_b_ids_right,
                 num_free_b_right, energy, stem_filling, pred=None):
        super(ExteriorFilling, self).__init__(
            bp_pos_id, stem_filling.bp_id, energy, pred)

        self.fixed_bp_id = fixed_bp_id
        bp_id_i, bp_id_j = stem_filling.bp_b_ids
        self.base_ids.append(bp_id_i)
        self.base_ids.append(bp_id_j)
        self.base_ids.extend(free_b_ids_right)
        self.num_free_b_right = num_free_b_right
        self.stem_filling = stem_filling

    def __repr__(self):
        if self.predecessor:
            pred = repr(self.predecessor)
        else:
            pred = ''

        if self.num_free_b_right:
            if len(self.base_ids) > 3:
                free_bases = ''.join((
                    rna.BASES[self.base_ids[2]],
                    rna.BASES[rna.BaseId.UNSPEC] * (self.num_free_b_right - 2),
                    rna.BASES[self.base_ids[3]]))
            else:
                free_bases = ''.join((
                    rna.BASES[self.base_ids[2]],
                    rna.BASES[rna.BaseId.UNSPEC] * (
                        self.num_free_b_right - 1)))
        else:
            free_bases = ''

        return ''.join([repr(self.stem_filling), free_bases, pred])

    def get_seq_list(self, seq, struct):
        i, j = struct.bp_positions[self.stem_filling.bp_pos_id]
        if self.num_free_b_right:
            seq[j + 1] = rna.BASES[self.base_ids[2]]
            if len(self.base_ids) > 3:
                seq[j + self.num_free_b_right] = rna.BASES[self.base_ids[3]]
        self.stem_filling.get_seq_list(seq, struct)
        if self.predecessor:
            self.predecessor.get_seq_list(seq, struct)


class FinalFilling(Filling):
    """
    Final filling:

                             (stem_filling)      (pred)
                               .       .       .         .
                                .     .         .       .
                                 (i j)           (i' j')
         (free_b_id_left)  ..,i-1     j+1,..,i'-1       j'+1

                (num_free_b_left)

    """
    def __init__(self, bp_pos_id, bp_id, free_b_id_left, num_free_b_left,
                 energy, pred):
        super(FinalFilling, self).__init__(bp_pos_id, bp_id, energy, pred)

        bp_id_i, bp_id_j = self.bp_b_ids
        self.base_ids.append(bp_id_i)
        self.base_ids.append(bp_id_j)
        if num_free_b_left:
            self.base_ids.append(free_b_id_left)
        self.num_free_b_left = num_free_b_left

    def __repr__(self):
        if self.num_free_b_left:
            return ''.join((
                rna.BASES[rna.BaseId.UNSPEC] * (self.num_free_b_left - 1),
                rna.BASES[self.base_ids[2]],
                repr(self.predecessor)))
        else:
            return repr(self.predecessor)

    def get_seq_list(self, seq, struct):
        i, j = struct.bp_positions[-1]
        if self.num_free_b_left:
            seq[i - 1] = rna.BASES[self.base_ids[2]]
        self.predecessor.get_seq_list(seq, struct)
