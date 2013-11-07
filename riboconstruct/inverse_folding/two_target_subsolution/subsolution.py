import itertools

from . import filling as fl

from .. import subsolution

from ... import rna


class TwoTargetSubsolution(subsolution.Subsolution):
    def __init__(self, struct_id, substruct, task, predecessor=None):
        self.struct_id = struct_id  # TODO: get rid of when no longer needed
        self.substruct = substruct
        self.predecessor = predecessor

        # TODO: get rid of self.struct
        self.struct = task.target_structs[struct_id]
        self.seq_constraint = task.seq_constraint
        self.dependencies = task.dependencies(struct_id)
        self.opposite_subsolutions = task.opposite_subsolutions(struct_id)
        self.b_pos_mapping = task.b_pos_mapping(struct_id)

        self.collection = TwoTargetCollection()

        bp_pos_i, bp_pos_j = self.bp_b_positions
        self.dep_b_positions = (
            self.dependencies.intersection(xrange(bp_pos_i, bp_pos_j + 1)))

        self.is_dependent = len(self.dep_b_positions)

    def num_fillings(self):
        return len(self.collection)

    def b_constrained(self, b_pos, b_id):
        return self.seq_constraint[b_pos][b_id]

    def bp_constrained(self, bp_pos_i, bp_pos_j, bp_id_i, bp_id_j):
        return (self.seq_constraint[bp_pos_i][bp_id_i] or
                self.seq_constraint[bp_pos_j][bp_id_j])

    def opp_struct_b_constrained(self, filling, b_pos, b_id):
        return False
        for opposite_subsolution, common_dep_positions in self.iter_common_dep_positions():
            b_pos_2 = opposite_subsolution.struct.basepairs[b_pos]
            if b_pos_2 is None or b_pos_2 not in common_dep_positions:
                continue
            pos_i, pos_j = self.bp_b_positions
            if b_pos_2 < pos_i or pos_j < b_pos_2:
                continue
            b_id_2 = self.get_base(filling, b_pos_2)
            return not rna.valid_bp_b_ids(b_id, b_id_2)
        return False

    def iter_common_dep_positions(self):
        for opposite_subsolution in self.opposite_subsolutions:
            common_dep_positions = (
                self.dep_b_positions.intersection(
                    opposite_subsolution.dep_b_positions))
            if len(common_dep_positions):
                yield opposite_subsolution, common_dep_positions
                self.clear_equality_keys()
                opposite_subsolution.clear_equality_keys()

    def depends_on(self, b_pos):
        return b_pos in self.dependencies

    def get_equality_key(self, filling, b_positions):
        try:
            return self.__equality_keys[filling]
        except AttributeError:
            self.__equality_keys = {}
        except KeyError:
            pass
        key = self.calc_equality_key(filling, b_positions)
        self.__equality_keys[filling] = key
        return key

    def clear_equality_keys(self):
        try:
            del self.__equality_keys
        except AttributeError:
            pass

    def get_base(self, filling, b_pos):
        filling_b_id, bp_pos_id = self.b_pos_mapping[b_pos]
        if filling.bp_pos_id == bp_pos_id:
            return filling.base_ids[filling_b_id]
        if (type(filling) == fl.TwoTargetExteriorFilling or
            type(filling) == fl.TwoTargetMultiloopStemFilling):
            # b_pos somewhere around the stem
            if filling.stem_filling.bp_pos_id == bp_pos_id:
                return filling.base_ids[filling_b_id]
            elif filling.predecessor:
                next_stem_bp_pos_id = filling.predecessor.stem_filling.bp_pos_id
                # b_pos somewhere within the stem
                if bp_pos_id > next_stem_bp_pos_id:
                    return self.get_base(filling.stem_filling, b_pos)
                # b_pos somewhere right of the stem
                return self.get_base(filling.predecessor, b_pos)
            # b_pos somewhere within the stem
            return self.get_base(filling.stem_filling, b_pos)
        return self.get_base(filling.predecessor, b_pos)

    def calc_equality_key(self, filling, b_positions):
        min_b_pos = min(b_positions)
        return sum([self.get_base(filling, b_pos) * 4 ** (b_pos - min_b_pos)
                    for b_pos in b_positions])

    def delete_predecessor(self):
        if self.is_dependent:
            for opposite_subsolution, common_dep_positions in self.iter_common_dep_positions():
                if not len(common_dep_positions):
                    continue
                common_dep_positions = tuple(sorted(common_dep_positions))

                for filling in self.get_fillings():
                    if not filling.unopt:
                        continue

                    key = self.get_equality_key(filling, common_dep_positions)

                    delete = False
                    no_matching_filling = True
                    for opp_filling in opposite_subsolution.iter_fillings():
                        opp_key = (
                            opposite_subsolution.get_equality_key(
                                opp_filling, common_dep_positions))

                        if key == opp_key:
                            no_matching_filling = False
                            if opp_filling.unopt and filling.unopt:
                                delete = True
                            # if at least one opposite filling is not unopt than
                            # do not delete
                            else:
                                delete = False
                                break
                    if no_matching_filling or delete:
                        self.delete_filling(filling)

                for opp_filling in opposite_subsolution.get_fillings():
                    if not opp_filling.unopt:
                        continue

                    opp_key = (
                        opposite_subsolution.get_equality_key(
                            opp_filling, common_dep_positions))

                    delete = False
                    no_matching_filling = True
                    for filling in self.iter_fillings():
                        key = (
                            self.get_equality_key(
                                filling, common_dep_positions))

                        if opp_key == key:
                            no_matching_filling = False
                            if opp_filling.unopt and filling.unopt:
                                delete = True
                            # if at least one opposite filling is not unopt than
                            # do not delete
                            else:
                                delete = False
                                break

                    if no_matching_filling or delete:
                        opposite_subsolution.delete_filling(opp_filling)

                self.clear_equality_keys()

        del self.collection.min_fillings
        super(TwoTargetSubsolution, self).delete_predecessor()

    def delete_filling(self, filling):
        self.collection[filling.bp_id][filling.index].remove(filling)


class TwoTargetCollection(subsolution.Collection):
    def __init__(self):
        self.collection = [{} for _ in xrange(rna.BasepairId.count)]
        self.min_fillings = [[] for _ in xrange(rna.BasepairId.count)]

    def __len__(self):
        return sum(
            (len(self.collection[bp_id])
             for bp_id in xrange(rna.BasepairId.count)))

    def iter_fillings(self):
        return itertools.chain.from_iterable(
            itertools.chain.from_iterable(self.collection[bp_id].itervalues())
            for bp_id in xrange(rna.BasepairId.count))

    def iter_bp_id_fillings(self, bp_id):
        return itertools.chain.from_iterable(
            self.collection[bp_id].itervalues())

    def get_fillings(self):
        return list(self.iter_fillings())

    def append_min_filling(self, filling):
        bp_id = filling.bp_id
        index = filling.index

        try:
            collected_filling = self.collection[bp_id][index][0]
        except KeyError:
            self.collection[bp_id][index] = [filling]
        else:
            if filling.energy < collected_filling.energy:
                # NOTE: working? -- should work since matching fillings in the
                # opposite subsolution are sorted out when no match was found
                # for unopt_filling in self.collection[bp_id][index]:
                #     unopt_filling.unopt = True
                self.collection[bp_id][index] = [filling]
            # elif helper.float_equal(filling.energy, collected_filling.energy):
            #     self.append(filling)
            else:
                return

        self.set_minimum(filling)

    def append(self, filling):
        bp_id = filling.bp_id
        index = filling.index

        try:
            self.collection[bp_id][index].append(filling)
        except KeyError:
            self.collection[bp_id][index] = [filling]

        self.set_minimum(filling)

    def set_minimum(self, filling):
        bp_id = filling.bp_id

        try:
            min_filling = self.min_fillings[bp_id][0]
        except IndexError:
            filling.unopt = False
            self.min_fillings[bp_id] = [filling]
        else:
            if filling.energy < min_filling.energy:
                for unopt_filling in self.min_fillings[bp_id]:
                    unopt_filling.unopt = True
                filling.unopt = False
                self.min_fillings[bp_id] = [filling]
            # elif helper.float_equal(filling.energy, min_filling.energy):
            #     filling.unopt = False
            #     self.min_fillings[fixed_bp_id][bp_id].append(filling)


class TwoTargetStem(subsolution.Stem):
    def __init__(self, subsolution, stem_id, predecessor=None):
        self.subsolution = subsolution
        self.stem_id = stem_id
        self.predecessor = predecessor

        self.collection = TwoTargetStemCollection()

        bp_pos_i, bp_pos_j = self.stem_bp_b_positions
        if self.predecessor:
            self.dep_b_positions = (
                self.subsolution.dependencies.intersection(
                    xrange(bp_pos_i, bp_pos_j + self.num_free_bases_right + 1)).
                union(self.predecessor.dep_b_positions))
        else:
            self.dep_b_positions = (
                self.subsolution.dependencies.intersection(
                    xrange(bp_pos_i, bp_pos_j + self.num_free_bases_right + 1)))

        self.is_dependent = len(self.dep_b_positions)

    @property
    def b_pos_mapping(self):
        return self.subsolution.b_pos_mapping

    def opp_struct_b_constrained(self, filling, b_pos, b_id):
        return self.subsolution.opp_struct_b_constrained(filling, b_pos, b_id)

    def iter_common_dep_positions(self):
        for opposite_subsolution in self.subsolution.opposite_subsolutions:
            common_dep_positions = (
                self.dep_b_positions.intersection(
                    opposite_subsolution.dep_b_positions))
            if len(common_dep_positions):
                yield opposite_subsolution, tuple(common_dep_positions)
                self.subsolution.clear_equality_keys()
                opposite_subsolution.clear_equality_keys()

    def depends_on(self, bp_pos):
        return self.subsolution.depends_on(bp_pos)

    def delete_predecessor(self):
        super(TwoTargetStem, self).delete_predecessor()
        return

        if self.is_dependent:
            for opposite_subsolution, common_dep_positions in self.iter_common_dep_positions():
                print '\t%i. stem: (s_id=%i, bp_pos=%s), opp ss: (s_id=%i, bp_pos=%s)' % (self.stem_id + 1, self.subsolution.struct_id, str(self.stem_bp_b_positions), opposite_subsolution.struct_id, str(opposite_subsolution.bp_b_positions))
                print '\tcommon positions:', common_dep_positions
                assert(len(list(self.iter_fillings())) == len(self.get_fillings()))
                print "\tclean up:"
                print "\t\tbefore", len(list(self.iter_fillings())), len(list(opposite_subsolution.iter_fillings()))

                if not len(common_dep_positions):
                    continue

                for filling in self.get_fillings():
                    if not filling.unopt:
                        continue

                    key = (
                        self.subsolution.get_equality_key(
                            filling, common_dep_positions))

                    delete = True
                    no_matching_filling = True
                    for opp_filling in opposite_subsolution.iter_fillings():
                        opp_key = (
                            opposite_subsolution.get_equality_key(
                                opp_filling, common_dep_positions))

                        if key == opp_key:
                            no_matching_filling = False
                            if opp_filling.unopt and filling.unopt:
                                delete = True
                            # if at least one opposite filling is not unopt than
                            # do not delete
                            else:
                                delete = False
                                break
                    if no_matching_filling or delete:
                        self.delete_filling(filling)

                for opp_filling in opposite_subsolution.get_fillings():
                    if not opp_filling.unopt:
                        continue

                    opp_key = (
                        opposite_subsolution.get_equality_key(
                            opp_filling, common_dep_positions))

                    delete = True
                    no_matching_filling = True
                    for filling in self.iter_fillings():
                        key = (
                            self.subsolution.get_equality_key(
                                filling, common_dep_positions))

                        if opp_key == key:
                            no_matching_filling = False
                            if opp_filling.unopt and filling.unopt:
                                delete = True
                            # if at least one opposite filling is not unopt than
                            # do not delete
                            else:
                                delete = False
                                break
                    if no_matching_filling or delete:
                        opposite_subsolution.delete_filling(opp_filling)

                bp_pos_i, bp_pos_j = opposite_subsolution.bp_b_positions
                for bp_id, (bp_id_i, bp_id_j) in enumerate(rna.BASEPAIRS):
                    if self.bp_constrained(bp_pos_i, bp_pos_j,
                                           bp_id_i, bp_id_j):
                        continue
                    if not sum(not filling.unopt for filling in opposite_subsolution.iter_bp_id_fillings(bp_id)):
                        return

                self.subsolution.clear_equality_keys()

            bp_pos_i, bp_pos_j = self.stem_bp_b_positions
            for bp_id, (bp_id_i, bp_id_j) in enumerate(rna.BASEPAIRS):
                if self.bp_constrained(bp_pos_i, bp_pos_j,
                                       bp_id_i, bp_id_j):
                    continue
                if not sum(not filling.unopt for filling in self.iter_bp_id_fillings(bp_id)):
                    return
                assert(sum(not filling.unopt for filling in self.iter_bp_id_fillings(bp_id)))

        del self.collection.min_fillings

    def delete_filling(self, filling):
        self.collection[filling.fixed_bp_id][filling.bp_id][
            filling.index].remove(filling)


class TwoTargetStemCollection(object):
    def __init__(self):
        self.collection = [[{} for _ in xrange(rna.BasepairId.count)]
                           for _ in xrange(rna.BasepairId.count)]
        self.min_fillings = [[[] for _ in xrange(rna.BasepairId.count)]
                             for _ in xrange(rna.BasepairId.count)]

    def __getitem__(self, fixed_bp_id):
        return self.collection[fixed_bp_id]

    def iter_fillings(self):
        for bp_id in xrange(rna.BasepairId.count):
            for filling in self.iter_bp_id_fillings(bp_id):
                yield filling

    def iter_bp_id_fillings(self, fixed_bp_id):
        return (
            itertools.chain.from_iterable(
                itertools.chain.from_iterable(
                    self.collection[fixed_bp_id][bp_id].itervalues())
                for bp_id in xrange(rna.BasepairId.count)))

    def get_fillings(self):
        return list(self.iter_fillings())

    def append_min_filling(self, filling):
        fixed_bp_id = filling.fixed_bp_id
        bp_id = filling.bp_id
        index = filling.index

        try:
            collected_filling = self.collection[fixed_bp_id][bp_id][index][0]
        except KeyError:
            self.collection[fixed_bp_id][bp_id][index] = [filling]
        else:
            if filling.energy < collected_filling.energy:
                # NOTE: working? -- should work since matching fillings in the
                # opposite subsolution are sorted out when no match was found
                # for unopt_filling in self.collection[fixed_bp_id][bp_id][index]:
                #     unopt_filling.unopt = True
                self.collection[fixed_bp_id][bp_id][index] = [filling]
            # elif helper.float_equal(filling.energy, collected_filling.energy):
            #     self.append(filling)
            else:
                return

        self.set_minimum(filling)

    def append(self, filling):
        fixed_bp_id = filling.fixed_bp_id
        bp_id = filling.bp_id
        index = filling.index

        try:
            self.collection[fixed_bp_id][bp_id][index].append(filling)
        except KeyError:
            self.collection[fixed_bp_id][bp_id][index] = [filling]

        self.set_minimum(filling)

    def set_minimum(self, filling):
        fixed_bp_id = filling.fixed_bp_id
        bp_id = filling.bp_id

        try:
            min_filling = self.min_fillings[fixed_bp_id][bp_id][0]
        except IndexError:
            filling.unopt = False
            self.min_fillings[fixed_bp_id][bp_id] = [filling]
        else:
            if filling.energy < min_filling.energy:
                for unopt_filling in self.min_fillings[fixed_bp_id][bp_id]:
                    unopt_filling.unopt = True
                filling.unopt = False
                self.min_fillings[fixed_bp_id][bp_id] = [filling]
            # elif helper.float_equal(filling.energy, min_filling.energy):
            #     filling.unopt = False
            #     self.min_fillings[fixed_bp_id][bp_id].append(filling)
