import itertools
import math

from . import energies
from . import generator

from ... import rna


def hairpin_energy(bp_pos_id, seq, struct):

    bp_pos_i, bp_pos_j = struct.bp_positions[bp_pos_id]
    bp_id_i, bp_id_j = seq[bp_pos_i], seq[bp_pos_j]
    # TODO: use inverse BasepairId?
    bp_id = getattr(rna.BasepairId, '%s%s' % (bp_id_i, bp_id_j))
    size = bp_id_j - bp_id_i - 1

    if size > 3:
        b_id_i = seq[bp_id_i + 1]
        b_id_j = seq[bp_id_j - 1]

        energy = (
            energies.MISMATCH_ENERGIES_HAIRPIN[bp_id_i][b_id_i][bp_id_j][
            b_id_j])

        if size == 4:
            embedded_seq = (
                ''.join(itertools.islice(seq, bp_pos_i + 1, bp_pos_j)))
            energy += energies.tetra_loop_energy(embedded_seq)

    elif size == 3:
        if (bp_id == rna.BasepairId.AU or
            bp_id == rna.BasepairId.GU or
            bp_id == rna.BasepairId.UA or
            bp_id == rna.BasepairId.UG):
            energy += rna.TERMINAL_AU

    else:
        energy = 0.0

        # TODO: special case 'CCC'

    if (struct.is_stem_end(bp_pos_id) and
        (bp_id == rna.BasepairId.AU or
         bp_id == rna.BasepairId.GU or
         bp_id == rna.BasepairId.UA or
         bp_id == rna.BasepairId.UG)):
        energy += rna.TERMINAL_AU

    return energy


class HairpinFillingsGenerator(generator.FillingsGenerator):
    def __init__(self, subsolution):
        super(HairpinFillingsGenerator, self).__init__(subsolution)

    def generate_fillings(self, bp_id):
        size = self.num_free_bases
        bp_pos_i, bp_pos_j = self.bp_b_positions
        bp_id_i, bp_id_j = rna.BASEPAIRS[bp_id]

        if self.bp_constrained(bp_pos_i, bp_pos_j, bp_id_i, bp_id_j):
            return

        if size <= 30:
            energy = energies.LOOP_DESTABILIZING_ENERGIES[size - 1][
                rna.StructType.HAIRPIN]
        else:
            energy = (
                energies.LOOP_DESTABILIZING_ENERGIES[30 - 1][self.struct_type])
            energy = energy + 1.75 * rna.RT * math.log(size / 30.0)
        if (self.is_stem_end and
            (bp_id == rna.BasepairId.AU or
             bp_id == rna.BasepairId.GU or
             bp_id == rna.BasepairId.UA or
             bp_id == rna.BasepairId.UG)):
            energy += rna.TERMINAL_AU

        if size > 3:
            if size != 4:
                for b_id_i in xrange(rna.BaseId.count):
                    if self.b_constrained(bp_pos_i + 1, b_id_i):
                        continue

                    for b_id_j in xrange(rna.BaseId.count):
                        if self.b_constrained(bp_pos_j - 1, b_id_j):
                            continue

                        energy_tmp = (
                            energy +
                            energies.MISMATCH_ENERGIES_HAIRPIN[bp_id_i][b_id_i][
                                bp_id_j][b_id_j])

                        self.append_min_filling(
                            self.new_filling(
                                self.bp_pos_id, bp_id, (b_id_i, b_id_j),
                                size, energy_tmp))

            else:
                b_pos_i1, b_pos_i2, b_pos_j2, b_pos_j1 = (
                    xrange(bp_pos_i + 1, bp_pos_j))
                for b_id_i1 in xrange(rna.BaseId.count):
                    if self.b_constrained(b_pos_i1, b_id_i1):
                        continue

                    for b_id_j1 in xrange(rna.BaseId.count):
                        if self.b_constrained(b_pos_j1, b_id_j1):
                            continue

                        mismatch_energy = (
                            energies.MISMATCH_ENERGIES_HAIRPIN[
                            bp_id_i][b_id_i1][bp_id_j][b_id_j1])
                        for b_id_i2 in xrange(rna.BaseId.count):
                            if self.b_constrained(b_pos_i2, b_id_i2):
                                continue

                            for b_id_j2 in xrange(rna.BaseId.count):
                                if self.b_constrained(b_pos_j2, b_id_j2):
                                    continue

                                embedded_seq = ''.join((rna.BASES[bp_id_i],
                                                        rna.BASES[b_id_i1],
                                                        rna.BASES[b_id_i2],
                                                        rna.BASES[b_id_j2],
                                                        rna.BASES[b_id_j1],
                                                        rna.BASES[bp_id_j]))
                                energy_tmp = (
                                    energy +
                                    energies.tetra_loop_energy(
                                        ''.join(embedded_seq)) +
                                    mismatch_energy)

                                self.append_min_filling(
                                    self.new_filling(
                                        self.bp_pos_id, bp_id,
                                        (b_id_i1, b_id_i2, b_id_j2, b_id_j1),
                                        size, energy_tmp))

        else:  # size == 3
            if (bp_id == rna.BasepairId.AU or
                bp_id == rna.BasepairId.GU or
                bp_id == rna.BasepairId.UA or
                bp_id == rna.BasepairId.UG):
                energy += rna.TERMINAL_AU

            # for b_id_1 in xrange(rna.BaseId.count):
            #     for b_id_2 in xrange(rna.BaseId.count):
            #         for b_id_3 in xrange(rna.BaseId.count):
            #             energy_tmp = energy

            #             if b_id_1 == b_id_2 == b_id_3 == rna.BaseId.C:
            #                 energy_tmp += rna.C_TRI_LOOP

            #             self.append_min_filling(
            #                 self.new_filling(
            #                     self.bp_pos_id, bp_id, (b_id_1, b_id_2, b_id_3),
            #                     size, energy_tmp))

            # NOTE: unless the three bases are not 'CCC' their energetical
            # contribution is equal: add '*' instead of actual base fillings
            self.append_min_filling(
                self.new_filling(
                    self.bp_pos_id, bp_id,
                    (rna.BaseId.UNSPEC, rna.BaseId.UNSPEC, rna.BaseId.UNSPEC),
                    size, energy))
