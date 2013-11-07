import math

from . import energies
from . import generator

from ... import rna


def bulge_energy(size, bp_id, pred_bp_id):
    bp_id_i, bp_id_j = rna.BASEPAIRS[bp_id]
    pred_bp_id_i, pred_bp_id_j = rna.BASEPAIRS[pred_bp_id]

    energy = 0.0

    # loop destabilizing energies
    if size <= 30:
        energy = energies.LOOP_DESTABILIZING_ENERGIES[size - 1][
            rna.StructType.BULGE]
    else:
        energy = (
            energies.LOOP_DESTABILIZING_ENERGIES[30 - 1][
                rna.StructType.BULGE])
        energy += 1.75 * rna.RT * math.log(size / 30.0)
    # bulge energy
    if size == 1:
        energy += (
            energies.STACKING_ENERGIES[bp_id_i][pred_bp_id_i][bp_id_j][
                pred_bp_id_j])
    else:
        # terminal energy
        if (bp_id == rna.BasepairId.AU or
            bp_id == rna.BasepairId.GU or
            bp_id == rna.BasepairId.UA or
            bp_id == rna.BasepairId.UG):
            energy += rna.TERMINAL_AU
        if (pred_bp_id == rna.BasepairId.AU or
            pred_bp_id == rna.BasepairId.GU or
            pred_bp_id == rna.BasepairId.UA or
            pred_bp_id == rna.BasepairId.UG):
            energy += rna.TERMINAL_AU

    return energy


class BulgeFillingsGenerator(generator.FillingsGenerator):
    def __init__(self, subsolution):
        super(BulgeFillingsGenerator, self).__init__(subsolution)

    def generate_fillings(self, bp_id):
        bp_pos_i, bp_pos_j = self.bp_b_positions
        bp_id_i, bp_id_j = rna.BASEPAIRS[bp_id]

        if self.bp_constrained(bp_pos_i, bp_pos_j, bp_id_i, bp_id_j):
            return

        size = self.num_free_bases
        # left bulge..
        if size[0]:
            size = size[0]
        # right bulge
        else:
            size = size[1]

        if size <= 30:
            energy = (
                energies.LOOP_DESTABILIZING_ENERGIES[size - 1][
                rna.StructType.BULGE])
        else:
            energy = (
                energies.LOOP_DESTABILIZING_ENERGIES[30 - 1][
                    rna.StructType.BULGE])
            energy += 1.75 * rna.RT * math.log(size / 30.0)
        if (bp_id == rna.BasepairId.AU or
            bp_id == rna.BasepairId.GU or
            bp_id == rna.BasepairId.UA or
            bp_id == rna.BasepairId.UG):
            if size == 1:
                energy += rna.TERMINAL_AU
            if self.is_stem_end:
                energy += rna.TERMINAL_AU

        for pred_filling in self.predecessor.iter_fillings():
            pred_bp_id_i, pred_bp_id_j = pred_filling.bp_b_ids

            energy_tmp = energy + pred_filling.energy
            if size == 1:
                energy_tmp += (
                    energies.STACKING_ENERGIES[
                        bp_id_i][pred_bp_id_i][bp_id_j][pred_bp_id_j])
            else:
                pred_bp_id = pred_filling.bp_id
                if (pred_bp_id == rna.BasepairId.AU or
                    pred_bp_id == rna.BasepairId.GU or
                    pred_bp_id == rna.BasepairId.UA or
                    pred_bp_id == rna.BasepairId.UG):
                    energy_tmp += rna.TERMINAL_AU

            self.append_min_filling(
                self.new_filling(
                    self.bp_pos_id, bp_id, self.num_free_bases, energy_tmp,
                    pred_filling))
