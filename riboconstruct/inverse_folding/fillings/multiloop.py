from . import energies
from . import generator

from ... import rna


def multiloop_energy(bp_pos_id, seq, struct):
    bp_pos_i, bp_pos_j = struct.bp_positions[bp_pos_id]
    stem_pos_ids = struct.predecessor_pos_ids[bp_pos_id]

    bp_id_i, bp_id_j = seq[bp_pos_i], seq[bp_pos_j]

    # first stem
    stem_pos_id = stem_pos_ids[0]
    stem_bp_pos_i, stem_bp_pos_j = struct.bp_positions[stem_pos_id]
    stem_bp_id_i, stem_bp_id_j = seq[stem_bp_pos_i], seq[stem_bp_pos_j]

    size = bp_id_j - stem_bp_pos_j
    if size == 1:
        stem_b_pos_j = stem_bp_pos_j + 1
        stem_b_id_j = seq[stem_b_pos_j]

        energy_pred_fb = (
            energies.SINGLE_BASE_STACKING_ENERGY[1][
                bp_id_i][bp_id_j][stem_b_id_j])
        energy_fb = (
            energies.SINGLE_BASE_STACKING_ENERGY[0][
                stem_bp_id_j][stem_bp_id_i][stem_b_id_j])
        energy = energy_pred_fb if energy_pred_fb < energy_fb else energy_fb
    elif size > 1:
        b_pos_j = bp_pos_j - 1
        b_id_j = seq[b_pos_j]

        energy_pred_fb = (
            energies.SINGLE_BASE_STACKING_ENERGY[1][bp_id_i][bp_id_j][b_id_j])
        energy_fb = (
                energies.SINGLE_BASE_STACKING_ENERGY[0][
                    stem_bp_id_j][stem_bp_id_i][stem_b_id_j])
        energy = energy_pred_fb + energy_fb

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
                energies.SINGLE_BASE_STACKING_ENERGY[1][pred_stem_bp_id_j][
                pred_stem_bp_id_i][pred_stem_b_id_i])
            energy_fb = (
                energies.SINGLE_BASE_STACKING_ENERGY[0][stem_bp_id_j][
                stem_bp_id_i][stem_b_id_j])
            energy += energy_fb + energy_pred_fb

    # the ml closing bp
    (pred_stem_bp_pos_i,
     pred_stem_bp_pos_j) = struct.bp_positions[stem_pos_ids[-1]]
    pred_stem_bp_id_i = seq[pred_stem_bp_pos_i]
    pred_stem_bp_id_j = seq[pred_stem_bp_pos_j]
    b_pos_i = bp_pos_i + 1
    b_id_i = seq[b_pos_i]

    size = pred_stem_bp_pos_i - bp_pos_i
    if size == 1:
        b_pos_i = bp_pos_i + 1
        b_id_i = seq[b_pos_i]

        energy_pred_fb = (
            energies.SINGLE_BASE_STACKING_ENERGY[1][
                pred_stem_bp_id_j][pred_stem_bp_id_i][b_id_i])
        energy_fb = (
            energies.SINGLE_BASE_STACKING_ENERGY[0][bp_id_i][bp_id_j][b_id_i])
        energy += energy_fb if energy_fb < energy_pred_fb else energy_pred_fb
    elif size > 1:
        b_pos_i = bp_pos_i + 1
        b_id_i = seq[b_pos_i]
        pred_stem_b_pos_i = pred_stem_bp_pos_i - 1
        pred_stem_b_id_i = seq[pred_stem_b_pos_i]

        energy_pred_fb = (
            energies.SINGLE_BASE_STACKING_ENERGY[1][
                pred_stem_bp_id_j][pred_stem_bp_id_i][pred_stem_b_id_i])
        energy_fb = (
            energies.SINGLE_BASE_STACKING_ENERGY[0][bp_id_i][bp_id_j][b_id_i])
        energy += energy_fb + energy_pred_fb

    return energy


class MultiloopFillingsGenerator(generator.FillingsGenerator):
    def __init__(self, subsolution):
        super(MultiloopFillingsGenerator, self).__init__(subsolution)

    def evaluate(self):
        self.generate_fillings()
        self.delete_predecessor()

    def generate_fillings(self):
        #######################################################################
        ### evaluate the first ml stem with fixed bp ml closing bp assignment #
        #######################################################################

        fillings_generator = (
            MultiloopStemFillingsGenerator(self.get_multiloop_stem(0)))
        predecessor = self.predecessor[0]
        bp_pos_i, bp_pos_j = self.bp_b_positions

        # at least one dangling free base at the end of the sequence
        if self.num_free_bases[0]:
            stem_b_pos_j = predecessor.bp_b_positions[1] + 1

            if self.num_free_bases[0] == 1:
                for fixed_bp_id, (fixed_bp_id_i, fixed_bp_id_j) in \
                    enumerate(rna.BASEPAIRS):
                    if self.bp_constrained(bp_pos_i, bp_pos_j,
                                           fixed_bp_id_i, fixed_bp_id_j):
                        continue

                    for stem_bp_id, (stem_bp_id_i, stem_bp_id_j) in \
                        enumerate(rna.BASEPAIRS):

                        for stem_b_id_j in xrange(rna.BaseId.count):
                            if self.b_constrained(stem_b_pos_j, stem_b_id_j):
                                continue

                            energy_pred_fb = (
                                energies.SINGLE_BASE_STACKING_ENERGY[1][
                                fixed_bp_id_i][fixed_bp_id_j][stem_b_id_j])
                            energy_fb = (
                                energies.SINGLE_BASE_STACKING_ENERGY[0][
                                stem_bp_id_j][stem_bp_id_i][stem_b_id_j])
                            energy_fb = (energy_fb
                                         if energy_fb < energy_pred_fb
                                         else energy_pred_fb)
                            energy_all = energy_fb

                            for stem_filling in \
                                predecessor.iter_bp_id_fillings(stem_bp_id):
                                fillings_generator.append_min_filling(
                                    self.new_multiloop_stem_filling(
                                        fixed_bp_id, self.bp_pos_id,
                                        (stem_b_id_j,), 1, energy_all,
                                        stem_filling))
            # num_free_bases > 1
            else:
                b_pos_j = bp_pos_j - 1

                for fixed_bp_id, (fixed_bp_id_i, fixed_bp_id_j) in \
                    enumerate(rna.BASEPAIRS):
                    if self.bp_constrained(bp_pos_i, bp_pos_j,
                                           fixed_bp_id_i, fixed_bp_id_j):
                        continue

                    for stem_bp_id, (stem_bp_id_i, stem_bp_id_j) in\
                        enumerate(rna.BASEPAIRS):

                        for b_id_j in xrange(rna.BaseId.count):
                            if self.b_constrained(b_pos_j, b_id_j):
                                continue

                            energy_pred_fb = (
                                energies.SINGLE_BASE_STACKING_ENERGY[1][
                                fixed_bp_id_i][fixed_bp_id_j][b_id_j])
                            for stem_b_id_j in xrange(rna.BaseId.count):
                                if self.b_constrained(stem_b_pos_j,
                                                      stem_b_id_j):
                                    continue

                                energy_fb = (
                                    energies.SINGLE_BASE_STACKING_ENERGY[0][
                                    stem_bp_id_j][stem_bp_id_i][stem_b_id_j])
                                energy_all = energy_pred_fb + energy_fb

                                for stem_filling in \
                                    predecessor.iter_bp_id_fillings(stem_bp_id):
                                    fillings_generator.append_min_filling(
                                        self.new_multiloop_stem_filling(
                                            fixed_bp_id, self.bp_pos_id,
                                            (stem_b_id_j, b_id_j),
                                            self.num_free_bases[0], energy_all,
                                            stem_filling))

        # no dangling free base at the end
        else:
            for fixed_bp_id, (fixed_bp_id_i, fixed_bp_id_j) in \
                enumerate(rna.BASEPAIRS):
                if self.bp_constrained(bp_pos_i, bp_pos_j,
                                       fixed_bp_id_i, fixed_bp_id_j):
                    continue

                energy_all = 0.0
                for stem_bp_id in xrange(rna.BasepairId.count):
                    for stem_filling in \
                        predecessor.iter_bp_id_fillings(stem_bp_id):
                        fillings_generator.append(
                            self.new_multiloop_stem_filling(
                                fixed_bp_id, self.bp_pos_id, (None,), 0,
                                energy_all, stem_filling))

        fillings_generator.delete_predecessor()

        #######################################################################
        ### evaluate all but the first ml stem (which has been evaluated ######
        ### before)                                                      ######
        #######################################################################

        for i, stem_subsolution in enumerate(self.predecessor[1:], 1):
            new_fillings_generator = (
                MultiloopStemFillingsGenerator(
                    self.get_multiloop_stem(
                        i, fillings_generator.multiloop_stem)))

            for fixed_bp_id in xrange(rna.BasepairId.count):
                new_fillings_generator.evaluate(fixed_bp_id)
            new_fillings_generator.delete_predecessor()

            fillings_generator = new_fillings_generator

        #######################################################################
        ### evaluate the final ml energy (for the fixed ml closing bp) ########
        #######################################################################

        offset = 3.4
        free_base_penalty = 0.0
        helix_penalty = 0.4
        energy_ml = (offset +
                     free_base_penalty * sum(self.num_free_bases) +
                     helix_penalty * (len(self.predecessor) + 1))

        pred_fillings = fillings_generator.iter_bp_id_fillings

        # at least one dangling free base left of the ml closing bp (right side
        # has been evaluated at the beginning)
        if self.num_free_bases[-1]:
            b_pos_i = bp_pos_i + 1

            # exactly one free base between ml closing bp and previous ml stem
            if self.num_free_bases[-1] == 1:
                for fixed_bp_id, (fixed_bp_id_i, fixed_bp_id_j) in \
                    enumerate(rna.BASEPAIRS):
                    if (fixed_bp_id == rna.BasepairId.AU or
                        fixed_bp_id == rna.BasepairId.GU or
                        fixed_bp_id == rna.BasepairId.UA or
                        fixed_bp_id == rna.BasepairId.UG):
                        energy_terminal = rna.TERMINAL_AU
                        if self.is_stem_end:
                            energy_terminal += rna.TERMINAL_AU
                    else:
                        energy_terminal = 0.0

                    for pred_filling in pred_fillings(fixed_bp_id):
                        (pred_bp_id_i,
                         pred_bp_id_j) = rna.BASEPAIRS[pred_filling.bp_id]

                        energy_pred = (energy_ml + energy_terminal +
                                       pred_filling.energy +
                                       pred_filling.stem_filling.energy)

                        for b_id_i in xrange(rna.BaseId.count):
                            if self.b_constrained(b_pos_i, b_id_i):
                                continue

                            energy_pred_fb = (
                                energies.SINGLE_BASE_STACKING_ENERGY[1][
                                pred_bp_id_j][pred_bp_id_i][b_id_i])
                            # free base adjacent to ml closing bp
                            energy_fb = (
                                energies.SINGLE_BASE_STACKING_ENERGY[0][
                                fixed_bp_id_i][fixed_bp_id_j][b_id_i])
                            energy_fb = (energy_fb
                                         if energy_fb < energy_pred_fb
                                         else energy_pred_fb)
                            energy_all = energy_pred + energy_fb

                            self.append_min_filling(
                                self.new_filling(
                                    self.bp_pos_id, fixed_bp_id, (b_id_i,), 1,
                                    energy_all, pred_filling))

            # more than one free base between ml closing bp and previous ml
            # stem
            else:
                pred_stem_b_pos_i = self.predecessor[-1].bp_b_positions[0] - 1

                for fixed_bp_id, (fixed_bp_id_i, fixed_bp_id_j) in \
                    enumerate(rna.BASEPAIRS):
                    if (fixed_bp_id == rna.BasepairId.AU or
                        fixed_bp_id == rna.BasepairId.GU or
                        fixed_bp_id == rna.BasepairId.UA or
                        fixed_bp_id == rna.BasepairId.UG):
                        energy_terminal = rna.TERMINAL_AU
                        if self.is_stem_end:
                            energy_terminal += rna.TERMINAL_AU
                    else:
                        energy_terminal = 0.0

                    for pred_filling in pred_fillings(fixed_bp_id):
                        (pred_bp_id_i,
                         pred_bp_id_j) = rna.BASEPAIRS[pred_filling.bp_id]

                        energy_pred = (energy_ml + energy_terminal +
                                       pred_filling.energy +
                                       pred_filling.stem_filling.energy)

                        for pred_stem_b_id_i in xrange(rna.BaseId.count):
                            if self.b_constrained(pred_stem_b_pos_i,
                                                  pred_stem_b_id_i):
                                continue

                            energy_pred_fb = (
                                energies.SINGLE_BASE_STACKING_ENERGY[1][
                                pred_bp_id_j][pred_bp_id_i][pred_stem_b_id_i])
                            for b_id_i in xrange(rna.BaseId.count):
                                if self.b_constrained(b_pos_i, b_id_i):
                                    continue

                                # free base adjacent to ml closing bp
                                energy_fb = (
                                    energies.SINGLE_BASE_STACKING_ENERGY[0][
                                    fixed_bp_id_i][fixed_bp_id_j][b_id_i])
                                energy_all = (
                                    energy_pred + energy_pred_fb + energy_fb)

                                self.append_min_filling(
                                    self.new_filling(
                                        self.bp_pos_id, fixed_bp_id,
                                        (b_id_i, pred_stem_b_id_i),
                                        self.num_free_bases[-1], energy_all,
                                        pred_filling))

        # no dangling free base left of the ml closing bp
        else:
            for fixed_bp_id, (fixed_bp_id_i, fixed_bp_id_j) in \
                enumerate(rna.BASEPAIRS):
                if (fixed_bp_id == rna.BasepairId.AU or
                    fixed_bp_id == rna.BasepairId.GU or
                    fixed_bp_id == rna.BasepairId.UA or
                    fixed_bp_id == rna.BasepairId.UG):
                    energy_terminal = rna.TERMINAL_AU
                    if self.is_stem_end:
                        energy_terminal += rna.TERMINAL_AU
                else:
                    energy_terminal = 0.0

                for pred_filling in pred_fillings(fixed_bp_id):
                    energy_all = (energy_ml + energy_terminal +
                                  pred_filling.energy +
                                  pred_filling.stem_filling.energy)

                    self.append_min_filling(
                        self.new_filling(
                            self.bp_pos_id, fixed_bp_id, (None,), 0, energy_all,
                            pred_filling))


class MultiloopStemFillingsGenerator(generator.FillingsGenerator):
    def __init__(self, stem):
        super(MultiloopStemFillingsGenerator, self).__init__(stem)

    @property
    def multiloop_stem(self):
        ''' The multiloop stem instance the generator is wrapping. '''
        return self.subsolution

    def evaluate(self, fixed_bp_id):
        for bp_id in xrange(rna.BasepairId.count):
            # bp constraint already handled during stem generation, i.e.
            # stem_fillings(bp_id) will be empty if bp_id is constrained
            self.generate_fillings(fixed_bp_id, bp_id)

    def generate_fillings(self, fixed_bp_id, bp_id):
        bp_id_i, bp_id_j = rna.BASEPAIRS[bp_id]
        stem_fillings = self.stem_subsolution.iter_bp_id_fillings

        for pred_filling in self.predecessor.iter_bp_id_fillings(fixed_bp_id):
            pred_bp_id_i, pred_bp_id_j = rna.BASEPAIRS[pred_filling.bp_id]

            energy_pred = pred_filling.energy + pred_filling.stem_filling.energy

            # at least one free base between current and previous ml stem
            if self.num_free_bases_right:
                stem_b_pos_j = self.stem_bp_b_positions[1] + 1

                # exactly one free base between current and previous ml
                # stem
                if self.num_free_bases_right == 1:
                    for b_id_j in xrange(rna.BaseId.count):
                        if self.b_constrained(stem_b_pos_j, b_id_j):
                            continue

                        energy_pred_fb = (
                            energies.SINGLE_BASE_STACKING_ENERGY[1][
                            pred_bp_id_j][pred_bp_id_i][b_id_j])
                        energy_fb = (
                            energies.SINGLE_BASE_STACKING_ENERGY[0][bp_id_j][
                            bp_id_i][b_id_j])
                        energy_fb = (energy_fb
                                     if energy_fb < energy_pred_fb
                                     else energy_pred_fb)
                        energy_all = energy_pred + energy_fb

                        for stem_filling in stem_fillings(bp_id):
                            self.append_min_filling(
                                self.new_filling(
                                    fixed_bp_id, self.bp_pos_id,
                                    (b_id_j,), 1, energy_all,
                                    stem_filling, pred_filling))

                # more than one free base between current and previous ml
                # stem
                else:
                    pred_stem_b_pos_i = (
                        self.predecessor.stem_bp_b_positions[0] - 1)

                    for pred_stem_b_id_i in xrange(rna.BaseId.count):
                        if self.b_constrained(pred_stem_b_pos_i,
                                              pred_stem_b_id_i):
                            continue

                        energy_pred_fb = (
                            energies.SINGLE_BASE_STACKING_ENERGY[1][
                            pred_bp_id_j][pred_bp_id_i][pred_stem_b_id_i])
                        for b_id_j in xrange(rna.BaseId.count):
                            if self.b_constrained(stem_b_pos_j, b_id_j):
                                continue

                            energy_fb = (
                                energies.SINGLE_BASE_STACKING_ENERGY[0][
                                bp_id_j][bp_id_i][b_id_j])
                            energy_all = (
                                energy_pred + energy_pred_fb + energy_fb)

                            for stem_filling in stem_fillings(bp_id):
                                self.append_min_filling(
                                    self.new_filling(
                                        fixed_bp_id, self.bp_pos_id,
                                        (b_id_j, pred_stem_b_id_i),
                                        self.num_free_bases_right,
                                        energy_all, stem_filling,
                                        pred_filling))

            # no free base(s) between current and previous stem
            else:
                energy_all = energy_pred
                for stem_filling in stem_fillings(bp_id):
                    self.append_min_filling(
                        self.new_filling(
                            fixed_bp_id, self.bp_pos_id, (None,), 0, energy_all,
                            stem_filling, pred_filling))
