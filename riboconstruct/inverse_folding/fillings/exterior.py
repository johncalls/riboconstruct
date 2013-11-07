from . import energies
from . import generator

from ... import rna


def exterior_energy(seq, struct):
    stem_pos_ids = struct.predecessor_pos_ids[len(struct.bp_positions)]

    energy = 0.0

    # first stem
    stem_pos_id = stem_pos_ids[0]
    stem_bp_pos_i, stem_bp_pos_j = struct.bp_positions[stem_pos_id]
    stem_bp_id_i, stem_bp_id_j = seq[stem_bp_pos_i], seq[stem_bp_pos_j]

    size = len(struct) - stem_bp_pos_j - 1
    if size == 1:
        stem_b_pos_j = stem_bp_pos_j + 1
        stem_b_id_j = seq[stem_b_pos_j]

        energy_fb = (
            energies.SINGLE_BASE_STACKING_ENERGY[0][
                stem_bp_id_j][stem_bp_id_i][stem_b_id_j])
        energy = energy_fb

    # the following stems
    for stem_pos_id in stem_pos_ids[1:]:
        stem_bp_pos_i, stem_bp_pos_j = struct.bp_positions[stem_pos_id]
        stem_bp_id_i, stem_bp_id_j = seq[stem_bp_pos_i], seq[stem_bp_pos_j]

        (pred_stem_bp_pos_i,
         pred_stem_bp_pos_j) = struct.bp_positions[stem_pos_id - 1]
        pred_stem_bp_id_i = seq[pred_stem_bp_pos_i]
        pred_stem_bp_id_j = seq[pred_stem_bp_pos_j]

        size = pred_stem_bp_pos_i - stem_bp_pos_j
        if size == 1:
            stem_b_pos_j = stem_bp_pos_j + 1
            stem_b_id_j = seq[stem_b_pos_j]

            energy_pred_fb = (
                energies.SINGLE_BASE_STACKING_ENERGY[1][
                    pred_stem_bp_id_j][pred_stem_bp_id_i][stem_b_id_j])
            energy_fb = (
                energies.SINGLE_BASE_STACKING_ENERGY[
                    0][stem_bp_id_j][stem_bp_id_i][stem_b_id_j])
            energy += (
                energy_fb if energy_fb < energy_pred_fb else energy_pred_fb)

        elif size > 1:
            stem_b_pos_j = stem_bp_pos_j + 1
            stem_b_id_j = seq[stem_b_pos_j]
            pred_stem_b_pos_i = pred_stem_bp_pos_i - 1
            pred_stem_b_id_i = seq[pred_stem_b_pos_i]

            energy_pred_fb = (
                energies.SINGLE_BASE_STACKING_ENERGY[1][
                    pred_stem_bp_id_j][pred_stem_bp_id_i][pred_stem_b_id_i])
            energy_fb = (
                energies.SINGLE_BASE_STACKING_ENERGY[0][
                    stem_bp_id_j][stem_bp_id_i][stem_b_id_j])
            energy += energy_fb + energy_pred_fb

    # the last stem
    stem_pos_id = stem_pos_ids[-1]
    stem_bp_pos_i, stem_bp_pos_j = struct.bp_positions[stem_pos_id]
    stem_bp_id_i, stem_bp_id_j = seq[stem_bp_pos_i], seq[stem_bp_pos_j]

    size = stem_bp_pos_i
    if size == 1:
        stem_b_pos_i = stem_bp_pos_i - 1
        stem_b_id_i = seq[stem_b_pos_i]
        energy_fb = energies.SINGLE_BASE_STACKING_ENERGY[1][
            stem_bp_id_j][stem_bp_id_i][stem_b_id_i]
        energy += energy_fb

    return energy


class FinalFillingsGenerator(generator.FillingsGenerator):
    def __init__(self, subsolution):
        super(FinalFillingsGenerator, self).__init__(subsolution)

    def evaluate(self):
        self.generate_fillings()
        self.delete_predecessor()

    def generate_fillings(self):
        #######################################################################
        ### evaluate the first exterior loop with fixed bp assignment #########
        #######################################################################

        fillings_generator = (
            ExteriorFillingsGenerator(self.get_exterior_stem(0)))
        stem = self.predecessor[0]
        fixed_bp_pos_i, fixed_bp_pos_j = stem.bp_b_positions

        # at least one dangling free base at the end of the sequence
        if self.num_free_bases[0]:
            stem_b_pos_j = stem.bp_b_positions[1] + 1

            for fixed_bp_id, (fixed_bp_id_i, fixed_bp_id_j) in \
                enumerate(rna.BASEPAIRS):
                if self.bp_constrained(fixed_bp_pos_i, fixed_bp_pos_j,
                                       fixed_bp_id_i, fixed_bp_id_j):
                    continue

                for stem_b_id_j in xrange(rna.BaseId.count):
                    if self.b_constrained(stem_b_pos_j, stem_b_id_j):
                        continue

                    energy_fb = (
                        energies.SINGLE_BASE_STACKING_ENERGY[0][fixed_bp_id_j][
                        fixed_bp_id_i][stem_b_id_j])
                    energy_all = energy_fb

                    for stem_filling in stem.iter_bp_id_fillings(fixed_bp_id):
                        fillings_generator.append_min_filling(
                            self.new_exterior_stem_filling(
                                fixed_bp_id, self.bp_pos_id, (stem_b_id_j,),
                                self.num_free_bases[0], energy_all,
                                stem_filling))

        # no dangling free base at the end
        else:
            for fixed_bp_id, (fixed_bp_id_i, fixed_bp_id_j) in \
                enumerate(rna.BASEPAIRS):
                if self.bp_constrained(fixed_bp_pos_i, fixed_bp_pos_j,
                                       fixed_bp_id_i, fixed_bp_id_j):
                    continue

                energy_all = 0.0
                for stem_filling in stem.iter_bp_id_fillings(fixed_bp_id):
                    fillings_generator.append(
                        self.new_exterior_stem_filling(
                            fixed_bp_id, self.bp_pos_id, (None,), 0, energy_all,
                            stem_filling))

        fillings_generator.delete_predecessor()

        #######################################################################
        ### evaluate all but the first exterior loop (which has been ##########
        ### evaluated before)                                        ##########
        #######################################################################

        for i, stem_subsolution in enumerate(self.predecessor[1:], 1):
            new_fillings_generator = (
                ExteriorFillingsGenerator(
                    self.get_exterior_stem(
                        i, fillings_generator.exterior_stem)))

            new_fillings_generator.evaluate()
            new_fillings_generator.delete_predecessor()

            fillings_generator = new_fillings_generator

        #######################################################################
        ### evaluate the final best energies (for the fixed first bp ##########
        ### assignment)                                              ##########
        #######################################################################

        pred_fillings = fillings_generator.iter_bp_id_fillings

        # at least one dangling free base at the beginning of the sequence
        if self.num_free_bases[-1]:
            stem_b_pos_i = self.predecessor[-1].bp_b_positions[0] - 1

            for fixed_bp_id in xrange(rna.BasepairId.count):
                for pred_filling in pred_fillings(fixed_bp_id):
                    (stem_bp_id_i,
                     stem_bp_id_j) = rna.BASEPAIRS[pred_filling.bp_id]

                    for stem_b_id_i in xrange(rna.BaseId.count):
                        if self.b_constrained(stem_b_pos_i, stem_b_id_i):
                            continue

                        energy_pred = (
                            pred_filling.energy +
                            pred_filling.stem_filling.energy)
                        energy_fb = (
                            energies.SINGLE_BASE_STACKING_ENERGY[1][
                            stem_bp_id_j][stem_bp_id_i][stem_b_id_i])
                        energy_all = energy_pred + energy_fb

                        self.append_min_filling(
                            self.new_filling(
                                self.bp_pos_id, fixed_bp_id, stem_b_id_i,
                                self.num_free_bases[-1], energy_all,
                                pred_filling))

        # no dangling free base at the end
        else:
            for fixed_bp_id in xrange(rna.BasepairId.count):
                for pred_filling in pred_fillings(fixed_bp_id):
                    energy_all = (
                        pred_filling.energy + pred_filling.stem_filling.energy)

                    self.append_min_filling(
                        self.new_filling(
                            self.bp_pos_id, fixed_bp_id, None, 0, energy_all,
                            pred_filling))


class ExteriorFillingsGenerator(generator.FillingsGenerator):
    def __init__(self, stem):
        super(ExteriorFillingsGenerator, self).__init__(stem)

    @property
    def exterior_stem(self):
        ''' The exterior stem instance the generator is wrapping. '''
        return self.subsolution

    def evaluate(self):
        for fixed_bp_id in xrange(rna.BasepairId.count):
            for bp_id in xrange(rna.BasepairId.count):
                # bp constraint already handled during stem generation, i.e.
                # stem_fillings(bp_id) will be empty if bp_id is constrained
                self.generate_fillings(fixed_bp_id, bp_id)

    def generate_fillings(self, fixed_bp_id, bp_id):
        bp_id_i, bp_id_j = rna.BASEPAIRS[bp_id]
        bp_pos_i, bp_pos_j = self.stem_subsolution.bp_b_positions

        stem_fillings = self.stem_subsolution.iter_bp_id_fillings

        for pred_filling in self.predecessor.iter_bp_id_fillings(fixed_bp_id):
            pred_bp_id_i, pred_bp_id_j = rna.BASEPAIRS[pred_filling.bp_id]

            pred_energy = pred_filling.energy + pred_filling.stem_filling.energy

            # at least one free base between current and previous stem
            if self.num_free_bases_right:
                stem_b_pos_j = self.stem_bp_b_positions[1] + 1

                # exactly one free base between current and previous stem
                if self.num_free_bases_right == 1:
                    for stem_b_id_j in xrange(rna.BaseId.count):
                        if self.b_constrained(stem_b_pos_j, stem_b_id_j):
                            continue

                        energy_pred_fb = (
                            energies.SINGLE_BASE_STACKING_ENERGY[1][
                            pred_bp_id_j][pred_bp_id_i][stem_b_id_j])
                        energy_fb = (
                            energies.SINGLE_BASE_STACKING_ENERGY[0][bp_id_j][
                            bp_id_i][stem_b_id_j])
                        energy_fb = (energy_fb
                                     if energy_fb < energy_pred_fb
                                     else energy_pred_fb)
                        energy_all = pred_energy + energy_fb

                        for stem_filling in stem_fillings(bp_id):
                            self.append_min_filling(
                                self.new_filling(
                                    fixed_bp_id, self.bp_pos_id, (stem_b_id_j,),
                                    1, energy_all, stem_filling, pred_filling))

                # more than one free base between current and previous stem
                else:
                    pred_stem_b_pos_i = (
                        self.predecessor.stem_bp_b_positions[0] - 1)

                    for pred_stem_b_id_i in xrange(rna.BaseId.count):
                        if self.b_constrained(
                            pred_stem_b_pos_i, pred_stem_b_id_i):
                            continue

                        energy_pred_fb = (
                            energies.SINGLE_BASE_STACKING_ENERGY[1][
                            pred_bp_id_j][pred_bp_id_i][pred_stem_b_id_i])
                        for stem_b_id_j in xrange(rna.BaseId.count):
                            if self.b_constrained(stem_b_pos_j, stem_b_id_j):
                                continue

                            energy_fb = (
                                energies.SINGLE_BASE_STACKING_ENERGY[0][
                                bp_id_j][bp_id_i][stem_b_id_j])
                            energy_all = (
                                pred_energy + energy_pred_fb + energy_fb)

                            for stem_filling in stem_fillings(bp_id):
                                self.append_min_filling(
                                    self.new_filling(
                                        fixed_bp_id, self.bp_pos_id,
                                        (stem_b_id_j, pred_stem_b_id_i),
                                        self.num_free_bases_right,
                                        energy_all, stem_filling,
                                        pred_filling))

            # no free base(s) between current and previous stem
            else:
                energy_all = pred_energy

                for stem_filling in stem_fillings(bp_id):
                    self.append_min_filling(
                        self.new_filling(
                            fixed_bp_id, self.bp_pos_id, (None,), 0, energy_all,
                            stem_filling, pred_filling))
