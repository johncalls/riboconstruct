import math

from . import energies
from . import generator

from ... import rna


def interior_energy(size_left, size_right, bp_id, pred_bp_id,
                    bp_pos_i, bp_pos_j, seq):
    if size_left == size_right == 1:
        b1 = seq[bp_pos_i + 1]
        b2 = seq[bp_pos_j - 1]
        return energies.INTERIOR_LOOP_1_1_ENERGY[bp_id][b1][pred_bp_id][b2]
    elif size_left == 1 and size_right == 2:
        b1 = seq[bp_pos_i + 1]
        b2 = seq[bp_pos_j - 1]
        b3 = seq[bp_pos_j - 2]
        return energies.INTERIOR_LOOP_1_2_ENERGY[bp_id][b2][b1][pred_bp_id][b3]
    elif size_left == 2 and size_right == 1:
        b1 = seq[bp_pos_i + 1]
        b2 = seq[bp_pos_i + 2]
        b3 = seq[bp_pos_j - 1]
        return energies.INTERIOR_LOOP_1_2_ENERGY[pred_bp_id][b1][b3][bp_id][b2]
    elif size_left == size_right == 2:
        b1 = seq[bp_pos_i + 1]
        b2 = seq[bp_pos_i + 2]
        b4 = seq[bp_pos_j - 2]
        b3 = seq[bp_pos_j - 1]
        return (
            energies.INTERIOR_LOOP_2_2_ENERGY[
                bp_id][pred_bp_id][b1][b3][b2][b4])
    else:
        bp_id_i, bp_id_j = rna.BASEPAIRS[bp_id]
        pred_bp_id_i, pred_bp_id_j = rna.BASEPAIRS[pred_bp_id]

        if size_left == 1:
            i1 = seq[bp_pos_i + 1]
            i2 = seq[bp_pos_j - 1]
            i3 = seq[bp_pos_j - size_right]
            return (
                energies.MISMATCH_ENERGIES_INTERIOR[
                    pred_bp_id_j][i2][pred_bp_id_i][i1] +
                energies.MISMATCH_ENERGIES_INTERIOR[bp_id_i][i1][bp_id_j][i3])
        elif size_right == 1:
            i1 = seq[bp_pos_i + 1]
            i2 = seq[bp_pos_i + size_left]
            i3 = seq[bp_pos_j - 1]
            return (
                energies.MISMATCH_ENERGIES_INTERIOR[
                    pred_bp_id_j][i3][pred_bp_id_i][i2] +
                energies.MISMATCH_ENERGIES_INTERIOR[bp_id_i][i1][bp_id_j][i3])
        else:
            i1 = seq[bp_pos_i + 1]
            i2 = seq[bp_pos_i + size_left]
            i3 = seq[bp_pos_j - size_right]
            i4 = seq[bp_pos_j - 1]
            return (
                energies.MISMATCH_ENERGIES_INTERIOR[bp_id_i][i1][bp_id_j][i4] +
                energies.MISMATCH_ENERGIES_INTERIOR[
                    pred_bp_id_j][i3][pred_bp_id_i][i2])


class InteriorFillingsGenerator(generator.FillingsGenerator):
    def __init__(self, subsolution):
        super(InteriorFillingsGenerator, self).__init__(subsolution)

    def generate_fillings(self, bp_id):
        bp_pos_i, bp_pos_j = self.bp_b_positions
        bp_id_i, bp_id_j = rna.BASEPAIRS[bp_id]

        if self.bp_constrained(bp_pos_i, bp_pos_j, bp_id_i, bp_id_j):
            return

        size_left, size_right = self.num_free_bases
        size = size_left + size_right

        if (self.is_stem_end and
            (bp_id == rna.BasepairId.AU or
             bp_id == rna.BasepairId.GU or
             bp_id == rna.BasepairId.UA or
             bp_id == rna.BasepairId.UG)):
            energy_start = rna.TERMINAL_AU
        else:
            energy_start = 0.0

        for pred_filling in self.predecessor.iter_fillings():
            pred_bp_id_i, pred_bp_id_j = rna.BASEPAIRS[pred_filling.bp_id]

            energy = pred_filling.energy + energy_start

            # special cases
            if size_left == size_right == 1:
                for b1 in xrange(rna.BaseId.count):
                    if self.b_constrained(bp_pos_i + 1, b1):
                        continue

                    for b2 in xrange(rna.BaseId.count):
                        if self.b_constrained(bp_pos_j - 1, b2):
                            continue

                        energy_tmp = (
                            energy +
                            energies.INTERIOR_LOOP_1_1_ENERGY[bp_id][b1][
                            pred_filling.bp_id][b2])

                        self.append_min_filling(
                            self.new_filling(
                                self.bp_pos_id, bp_id, (b1, b2), (1, 1),
                                energy_tmp, pred_filling))

            elif size_left == 1 and size_right == 2:
                for b2 in xrange(rna.BaseId.count):
                    if self.b_constrained(bp_pos_j - 2, b2):
                        continue

                    for b1 in xrange(rna.BaseId.count):
                        if self.b_constrained(bp_pos_i + 1, b1):
                            continue

                        for b3 in xrange(rna.BaseId.count):
                            if self.b_constrained(bp_pos_j - 1, b3):
                                continue

                            energy_tmp = (
                                energy +
                                energies.INTERIOR_LOOP_1_2_ENERGY[bp_id][b2][
                                b1][pred_filling.bp_id][b3])

                            self.append_min_filling(
                                self.new_filling(
                                    self.bp_pos_id, bp_id, (b1, b3, b2),
                                    (1, 2), energy_tmp, pred_filling))

            elif size_left == 2 and size_right == 1:
                for b1 in xrange(rna.BaseId.count):
                    if self.b_constrained(bp_pos_i + 1, b1):
                        continue

                    for b3 in xrange(rna.BaseId.count):
                        if self.b_constrained(bp_pos_j - 1, b3):
                            continue

                        for b2 in xrange(rna.BaseId.count):
                            if self.b_constrained(bp_pos_i + 2, b2):
                                continue

                            energy_tmp = (
                                energy +
                                energies.INTERIOR_LOOP_1_2_ENERGY[
                                pred_filling.bp_id][b1][b3][bp_id][b2])

                            self.append_min_filling(
                                self.new_filling(
                                    self.bp_pos_id, bp_id, (b1, b2, b3),
                                    (2, 1), energy_tmp, pred_filling))

            elif size_left == size_right == 2:
                for b1 in xrange(rna.BaseId.count):
                    if self.b_constrained(bp_pos_i + 1, b1):
                        continue

                    for b4 in xrange(rna.BaseId.count):
                        if self.b_constrained(bp_pos_j - 2, b4):
                            continue

                        for b2 in xrange(rna.BaseId.count):
                            if self.b_constrained(bp_pos_i + 2, b2):
                                continue

                            for b3 in xrange(rna.BaseId.count):
                                if self.b_constrained(bp_pos_j - 1, b3):
                                    continue

                                energy_tmp = (
                                    energy +
                                    energies.INTERIOR_LOOP_2_2_ENERGY[bp_id][
                                    pred_filling.bp_id][b1][b4][b2][b3])

                                self.append_min_filling(
                                    self.new_filling(
                                        self.bp_pos_id, bp_id,
                                        (b1, b2, b3, b4), (2, 2),
                                        energy_tmp, pred_filling))

            # standard case
            else:
                # asymmetry penalty
                if size_left != size_right:
                    asym = 0.5 * abs(size_left - size_right)
                    if asym > 3.0:
                        asym = 3.0
                    energy += asym

                # loop destabilizing energies
                if size <= 30:
                    energy = (
                        energy +
                        energies.LOOP_DESTABILIZING_ENERGIES[size - 1][
                        rna.StructType.INTERIOR])
                else:
                    energy = (
                        energy +
                        energies.LOOP_DESTABILIZING_ENERGIES[30 - 1][
                            rna.StructType.INTERIOR] +
                        1.75 * rna.RT * math.log(size / 30.0))

                # there are dependencies if the IL has size 1 at one side
                # of the loop: this base is involved in the terminal
                # mismatches at both closings
                #     (pred_bp_id_i, pred_bp_id_j)
                #                                  i2
                #  i1                              i3
                #          (bp_id_i, bp_id_j)
                if size_left == 1:
                    for i1 in xrange(rna.BaseId.count):
                        if self.b_constrained(bp_pos_i + 1, i1):
                            continue

                        for i2 in xrange(rna.BaseId.count):
                            if self.b_constrained(pred_bp_id_j + 1, i2):
                                continue

                            for i3 in xrange(rna.BaseId.count):
                                if self.b_constrained(bp_pos_j - 1, i3):
                                    continue

                                energy_tmp = (
                                    energy +
                                    energies.MISMATCH_ENERGIES_INTERIOR[
                                    pred_bp_id_j][i2][pred_bp_id_i][i1] +
                                    energies.MISMATCH_ENERGIES_INTERIOR[
                                    bp_id_i][i1][bp_id_j][i3])

                                self.append_min_filling(
                                    self.new_filling(
                                        self.bp_pos_id, bp_id, (i1, i2, i3),
                                        (1, size_right), energy_tmp,
                                        pred_filling))

                #     (pred_bp_id_i, pred_bp_id_j)
                #  i2                              i3
                #  i1
                #          (bp_id_i, bp_id_j)
                elif size_right == 1:
                    for i1 in xrange(rna.BaseId.count):
                        if self.b_constrained(bp_pos_i + 1, i1):
                            continue

                        for i3 in xrange(rna.BaseId.count):
                            if self.b_constrained(bp_pos_j - 1, i3):
                                continue

                            for i2 in xrange(rna.BaseId.count):
                                if self.b_constrained(pred_bp_id_i - 1, i2):
                                    continue

                                energy_tmp = (
                                    energy +
                                    energies.MISMATCH_ENERGIES_INTERIOR[
                                    pred_bp_id_j][i3][pred_bp_id_i][i2] +
                                    energies.MISMATCH_ENERGIES_INTERIOR[
                                    bp_id_i][i1][bp_id_j][i3])

                                self.append_min_filling(
                                    self.new_filling(
                                        self.bp_pos_id, bp_id, (i1, i2, i3),
                                        (size_left, 1), energy_tmp,
                                        pred_filling))

                #     (pred_bp_id_i, pred_bp_id_j)
                #  i2                              i3
                #  i1                              i4
                #          (bp_id_i, bp_id_j)
                else:
                    for i1 in xrange(rna.BaseId.count):
                        if self.b_constrained(bp_pos_i + 1, i1):
                            continue

                        for i4 in xrange(rna.BaseId.count):
                            if self.b_constrained(bp_pos_j + 1, i4):
                                continue

                            energy_tmp = (
                                energy +
                                energies.MISMATCH_ENERGIES_INTERIOR[bp_id_i][
                                i1][bp_id_j][i4])

                            self.append_min_filling(
                                self.new_filling(
                                    self.bp_pos_id, bp_id,
                                    (i1, rna.BaseId.UNSPEC, rna.BaseId.UNSPEC,
                                     i4),
                                    (size_left, size_right), energy_tmp,
                                    pred_filling))

                    # TODO: fix!?
                    tmp_fillings = list(self.iter_bp_id_fillings(bp_id))
                    self.collection.collection[bp_id] = {}

                    # NOTE: global energy added in the previous for-loop
                    # and stored in each tmp_fillings entry
                    for tmp_filling in tmp_fillings:
                        for i2 in xrange(rna.BaseId.count):
                            if self.b_constrained(pred_bp_id_i, i2):
                                continue

                            for i3 in xrange(rna.BaseId.count):
                                if self.b_constrained(pred_bp_id_j + 1, i3):
                                    continue

                                energy_tmp = (
                                    energies.MISMATCH_ENERGIES_INTERIOR[
                                    pred_bp_id_j][i3][pred_bp_id_i][i2] +
                                    tmp_filling.energy)

                                _, _, i1, _, _, i4 = tmp_filling.base_ids

                                self.append_min_filling(
                                    self.new_filling(
                                        self.bp_pos_id, bp_id,
                                        (i1, i2, i3, i4),
                                        (size_left, size_right),
                                        energy_tmp, pred_filling))
