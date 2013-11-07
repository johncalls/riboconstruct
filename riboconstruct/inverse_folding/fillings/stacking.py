from . import energies
from . import generator

from ... import rna


def stacking_energy(bp_id, pred_bp_id):
    bp_id_i, bp_id_j = rna.BASEPAIRS[bp_id]
    pred_bp_id_i, pred_bp_id_j = rna.BASEPAIRS[pred_bp_id]

    return (
        energies.STACKING_ENERGIES[
            bp_id_i][pred_bp_id_i][bp_id_j][pred_bp_id_j])


class StackFillingsGenerator(generator.FillingsGenerator):
    def __init__(self, subsolution):
        super(StackFillingsGenerator, self).__init__(subsolution)

    def generate_fillings(self, bp_id):
        bp_pos_i, bp_pos_j = self.bp_b_positions
        bp_id_i, bp_id_j = rna.BASEPAIRS[bp_id]

        if self.bp_constrained(bp_pos_i, bp_pos_j, bp_id_i, bp_id_j):
            return

        if (self.is_stem_end and
            (bp_id == rna.BasepairId.AU or
             bp_id == rna.BasepairId.GU or
             bp_id == rna.BasepairId.UA or
             bp_id == rna.BasepairId.UG)):
            energy = rna.TERMINAL_AU
        else:
            energy = 0.0

        for pred_filling in self.predecessor.iter_fillings():
            energy_tmp = (
                energy +
                stacking_energy(bp_id, pred_filling.bp_id) +
                pred_filling.energy)

            self.append_min_filling(
                self.new_filling(
                    self.bp_pos_id, bp_id, energy_tmp, pred_filling))
