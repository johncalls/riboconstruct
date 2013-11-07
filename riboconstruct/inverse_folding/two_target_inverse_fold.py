#!/usr/bin/env python

import itertools
import random


from . import dependency_graph
from . import inverse_fold
from . import structure
from . import two_target_local_refinement as loc_ref

from .fillings import bulge as fill_gen_bulge
from .fillings import exterior as fill_gen_exterior
from .fillings import hairpin as fill_gen_hairpin
from .fillings import interior as fill_gen_interior
from .fillings import multiloop as fill_gen_multiloop
from .fillings import stacking as fill_gen_stacking

from .two_target_subsolution import bulge
from .two_target_subsolution import exterior
from .two_target_subsolution import hairpin
from .two_target_subsolution import interior
from .two_target_subsolution import multiloop
from .two_target_subsolution import stacking

from .. import rna


MAX_RAND_B_ID = rna.BaseId.count - 1


def calculate_sequence_constraint(structs, iupac_seq):
    def matching_base_ids(bp_id_i):
        if bp_id_i == rna.BaseId.A:
            return (rna.BaseId.U,)
        elif bp_id_i == rna.BaseId.C:
            return (rna.BaseId.G,)
        elif bp_id_i == rna.BaseId.G:
            return (rna.BaseId.C, rna.BaseId.U)
        elif bp_id_i == rna.BaseId.U:
            return (rna.BaseId.A, rna.BaseId.G)
        else:
            return tuple()

    def rec_bases_mutation(pos_id_next, next_valid_base_ids):
        pos_id = pos_id_next
        valid_base_ids = next_valid_base_ids
        pos_id_next += 1
        iupac_id = getattr(rna.IUPAC_Id, iupac_seq[positions[pos_id]])
        if pos_id_next == len(positions):
            for b_id in valid_base_ids:
                if not rna.base_valid_IUPAC(iupac_id, b_id):
                    continue
                subseq[pos_id] = b_id
                yield subseq
            return
        for b_id in valid_base_ids:
            if not rna.base_valid_IUPAC(iupac_id, b_id):
                continue
            subseq[pos_id] = b_id
            for _ in rec_bases_mutation(pos_id_next,
                                        matching_base_ids(b_id)):
                yield subseq

    seq_constraint = [None] * len(iupac_seq)
    paired_positions = []

    dep_graph = dependency_graph.DependencyGraph()
    struct_0, struct_1 = structs
    for bp_pos in struct_0.bp_positions:
        dep_graph.add_edge(bp_pos)
    for bp_pos in struct_1.bp_positions:
        dep_graph.add_edge(bp_pos)

    for subgraph in dep_graph.iter_subgraphs():
        positions = list(subgraph.iter_nodes_sorted())
        paired_positions.extend(positions)
        for pos in positions:
            seq_constraint[pos] = [True] * rna.BaseId.count
        subseq = [None] * len(positions)
        iupac_id = getattr(rna.IUPAC_Id, iupac_seq[positions[0]])
        for b_id in xrange(rna.BaseId.count):
            if not rna.base_valid_IUPAC(iupac_id, b_id):
                continue
            subseq[0] = b_id
            for _ in rec_bases_mutation(1, matching_base_ids(b_id)):
                for pos, b_id in zip(positions, subseq):
                    seq_constraint[pos][b_id] = False

    for pos in set(xrange(len(iupac_seq))).difference(paired_positions):
        iupac_id = getattr(rna.IUPAC_Id, iupac_seq[pos])
        seq_constraint[pos] = [not rna.base_valid_IUPAC(iupac_id, b_id)
                               for b_id in xrange(rna.BaseId.count)]

    return seq_constraint


def calculate_struct_properties(struct):
    def add_positions_hairpin(substruct):
        if substruct.num_free_bases == 3:
            pass
        elif substruct.num_free_bases == 4:
            bp_pos_i, _ = substruct.bp_b_positions
            dependencies.add(bp_pos_i + 1)
            dependencies.add(bp_pos_i + 2)
            dependencies.add(bp_pos_i + 3)
            dependencies.add(bp_pos_i + 4)
        elif substruct.num_free_bases > 4:
            bp_pos_i, bp_pos_j = substruct.bp_b_positions
            dependencies.add(bp_pos_i + 1)
            dependencies.add(bp_pos_j - 1)

    def add_positions_bulge(substruct):
        pred_bp_pos_i, pred_bp_pos_j = (
            struct.bp_positions[substruct.predecessor_pos_ids[0]])
        dependencies.add(pred_bp_pos_i)
        dependencies.add(pred_bp_pos_j)

    def add_positions_stack(substruct):
        pred_bp_pos_i, pred_bp_pos_j = (
            struct.bp_positions[substruct.predecessor_pos_ids[0]])
        dependencies.add(pred_bp_pos_i)
        dependencies.add(pred_bp_pos_j)

    def add_positions_interior(substruct):
        bp_pos_i, bp_pos_j = substruct.bp_b_positions
        pred_bp_pos_i, pred_bp_pos_j = (
            struct.bp_positions[substruct.predecessor_pos_ids[0]])
        dependencies.add(pred_bp_pos_i)
        dependencies.add(pred_bp_pos_j)
        num_free_bases = substruct.num_free_bases
        dependencies.add(bp_pos_i + 1)
        if num_free_bases[0] >= 2:
            dependencies.add(pred_bp_pos_i - 1)
        dependencies.add(bp_pos_j - 1)
        if num_free_bases[1] >= 2:
            dependencies.add(pred_bp_pos_j - 1)

    def add_positions_multiloop(substruct):
        for stem_pos_id, num_free_bases in itertools.izip(
            substruct.predecessor_pos_ids, substruct.num_free_bases):
            stem_bp_pos_i, stem_bp_pos_j = struct.bp_positions[stem_pos_id]
            dependencies.add(stem_bp_pos_i)
            dependencies.add(stem_bp_pos_j)
            if num_free_bases > 0:
                dependencies.add(stem_bp_pos_j + 1)
                if num_free_bases > 1:
                    dependencies.add(stem_bp_pos_j + num_free_bases)
        bp_pos_i, _ = substruct.bp_b_positions
        num_free_bases = substruct.num_free_bases[-1]
        if num_free_bases > 0:
            dependencies.add(bp_pos_i + 1)
            if num_free_bases > 1:
                dependencies.add(bp_pos_i + num_free_bases)

    def add_positions_exterior(substruct):
        it = itertools.izip(substruct.predecessor_pos_ids,
                            substruct.num_free_bases)
        stem_pos_id, num_free_bases = it.next()
        stem_bp_pos_i, stem_bp_pos_j = struct.bp_positions[stem_pos_id]
        dependencies.add(stem_bp_pos_i)
        dependencies.add(stem_bp_pos_j)
        if num_free_bases > 0:
            dependencies.add(stem_bp_pos_j + 1)
        for stem_pos_id, num_free_bases in it:
            stem_bp_pos_i, stem_bp_pos_j = struct.bp_positions[stem_pos_id]
            dependencies.add(stem_bp_pos_i)
            dependencies.add(stem_bp_pos_j)
            dependencies.add(stem_bp_pos_j + 1)
            if num_free_bases > 1:
                dependencies.add(stem_bp_pos_j + num_free_bases)
        num_free_bases = substruct.num_free_bases[-1]
        if num_free_bases > 0:
            # NOTE: stem_bp_pos_i is defined in the for-loop which is
            # visited at least once
            dependencies.add(stem_bp_pos_i - 1)

    ########################################################################

    dependencies = set()

    for substruct in struct.iter_substructs():
        if substruct.struct_type == rna.StructType.HAIRPIN:
            add_positions_hairpin(substruct)
        elif substruct.struct_type == rna.StructType.STACKING:
            add_positions_stack(substruct)
        elif substruct.struct_type == rna.StructType.BULGE:
            add_positions_bulge(substruct)
        elif substruct.struct_type == rna.StructType.INTERIOR:
            add_positions_interior(substruct)
        elif substruct.struct_type == rna.StructType.MULTILOOP:
            add_positions_multiloop(substruct)
        elif substruct.struct_type == rna.StructType.EXTERIOR:
            add_positions_exterior(substruct)

    return dependencies


class Task(object):
    def __init__(self, target_structs, seq_constraint):
        self.target_structs = target_structs
        self.seq_constraint = seq_constraint

        self.dependent_positions = (
            calculate_struct_properties(self.target_structs[0]),
            calculate_struct_properties(self.target_structs[1]))
        self.__b_pos_mapping = ({}, {})

        self.subsolution_collector = [{}, {}]

    @property
    def subsolutions(self):
        return (self.subsolution_collector[0].values(),
                self.subsolution_collector[1].values())

    def seq_constraints(self, struct_id):
        return self.__seq_constraints[struct_id]

    def dependencies(self, struct_id):
        return self.dependent_positions[(struct_id + 1) & 1]

    def opposite_subsolutions(self, struct_id):
        return self.subsolution_collector[(struct_id + 1) & 1].values()

    def b_pos_mapping(self, struct_id):
        return self.__b_pos_mapping[struct_id]

    def iter_substructs(self):
        substruct_iter_0 = self.target_structs[0].iter_substructs()
        substruct_iter_1 = self.target_structs[1].iter_substructs()
        substruct_0 = substruct_iter_0.next()
        substruct_1 = substruct_iter_1.next()
        while 1:
            if substruct_0.struct_type == rna.StructType.EXTERIOR:
                while substruct_1.struct_type != rna.StructType.EXTERIOR:
                    yield 1, substruct_1
                    substruct_1 = substruct_iter_1.next()
                yield 0, substruct_0
                yield 1, substruct_1
                break
            if substruct_1.struct_type == rna.StructType.EXTERIOR:
                while substruct_0.struct_type != rna.StructType.EXTERIOR:
                    yield 0, substruct_0
                    substruct_0 = substruct_iter_0.next()
                yield 1, substruct_1
                yield 0, substruct_0
                break
            bp_pos_0 = substruct_0.bp_b_positions
            bp_pos_1 = substruct_1.bp_b_positions
            if bp_pos_0[0] > bp_pos_1[0]:
                yield 0, substruct_0
                substruct_0 = substruct_iter_0.next()
            elif bp_pos_0[0] < bp_pos_1[0]:
                yield 1, substruct_1
                substruct_1 = substruct_iter_1.next()
            # bp_pos_0[0] == bp_pos_1[0]
            else:
                if bp_pos_0[1] > bp_pos_1[1]:
                    yield 0, substruct_0
                    substruct_0 = substruct_iter_0.next()
                elif bp_pos_0[1] < bp_pos_1[1]:
                    yield 1, substruct_1
                    substruct_1 = substruct_iter_1.next()
                # bp_pos_0 == bp_pos_1
                else:
                    yield 0, substruct_0
                    substruct_0 = substruct_iter_0.next()
                    yield 1, substruct_1
                    substruct_1 = substruct_iter_1.next()
        # make sure both iterators have 'run dry'
        for substruct in substruct_iter_0:
            yield 0, substruct
        for substruct in substruct_iter_1:
            yield 1, substruct


def two_target_inverse_fold(target_structs, seq_constraint):
    def evaluate(struct_id, substruct):
        def get_fillings_generator(predecessor=None):
            if substruct.struct_type == rna.StructType.HAIRPIN:
                return fill_gen_hairpin.HairpinFillingsGenerator(
                    hairpin.TwoTargetHairpinSubsolution(
                        struct_id, substruct, task))
            elif substruct.struct_type == rna.StructType.BULGE:
                return fill_gen_bulge.BulgeFillingsGenerator(
                    bulge.TwoTargetBulgeSubsolution(
                        struct_id, substruct, task, predecessor))
            elif substruct.struct_type == rna.StructType.STACKING:
                return fill_gen_stacking.StackFillingsGenerator(
                    stacking.TwoTargetStackSubsolution(
                        struct_id, substruct, task, predecessor))
            elif substruct.struct_type == rna.StructType.INTERIOR:
                return fill_gen_interior.InteriorFillingsGenerator(
                    interior.TwoTargetInteriorSubsolution(
                        struct_id, substruct, task, predecessor))
            elif substruct.struct_type == rna.StructType.MULTILOOP:
                return fill_gen_multiloop.MultiloopFillingsGenerator(
                    multiloop.TwoTargetMultiloopSubsolution(
                        struct_id, substruct, task, predecessor))
            elif substruct.struct_type == rna.StructType.EXTERIOR:
                return fill_gen_exterior.FinalFillingsGenerator(
                    exterior.TwoTargetFinalSubsolution(
                        struct_id, substruct, task, predecessor))

        if substruct.struct_type == rna.StructType.MULTILOOP:
            pred = []
            for stem_bp_pos_id in substruct.predecessor_pos_ids:
                pred.append(
                    task.subsolution_collector[struct_id][stem_bp_pos_id])
                del task.subsolution_collector[struct_id][stem_bp_pos_id]
            fillings_generator = get_fillings_generator(pred)
            task.subsolution_collector[struct_id][substruct.bp_pos_id] = (
                fillings_generator.subsolution)
        elif substruct.struct_type == rna.StructType.EXTERIOR:
            pred = []
            for stem_bp_pos_id in substruct.predecessor_pos_ids:
                pred.append(
                    task.subsolution_collector[struct_id][stem_bp_pos_id])
                del task.subsolution_collector[struct_id][stem_bp_pos_id]
            fillings_generator = get_fillings_generator(pred)
            task.subsolution_collector[struct_id][substruct.bp_pos_id] = (
                fillings_generator.subsolution)
        elif substruct.struct_type == rna.StructType.HAIRPIN:
            fillings_generator = get_fillings_generator()
            task.subsolution_collector[struct_id][substruct.bp_pos_id] = (
                fillings_generator.subsolution)
        else:
            pred = (
                task.subsolution_collector[struct_id][substruct.bp_pos_id - 1])
            del task.subsolution_collector[struct_id][substruct.bp_pos_id - 1]
            fillings_generator = get_fillings_generator(pred)
            task.subsolution_collector[struct_id][substruct.bp_pos_id] = (
                fillings_generator.subsolution)

        fillings_generator.evaluate()

    ############################################################################

    task = Task(target_structs, seq_constraint)

    if task.target_structs[0] == task.target_structs[1]:
        subsol = inverse_fold.inverse_fold(task.target_structs[0])
        return (subsol, subsol)

    else:
        for bp_belonging, substruct in task.iter_substructs():
            evaluate(bp_belonging, substruct)

        subsol_0, subsol_1 = task.subsolutions
        # assert(len(subsol_0) == 1)
        # assert(len(subsol_1) == 1)
        subsol_0 = subsol_0[0]
        subsol_1 = subsol_1[0]

        return task.subsolutions[0][0], task.subsolutions[1][0]


def iter_combined_seq_lists(subsolutions, target_structs, seq_constraint):
        subsol_0, subsol_1 = subsolutions

        check_doubles = set()

        common_dep_positions = (
            subsol_0.dep_b_positions.intersection(subsol_1.dep_b_positions))
        seq_list = [rna.BASES[rna.BaseId.UNSPEC]] * len(seq_constraint)
        for filling in subsol_0.iter_fillings():
            filling.get_seq_list(seq_list, target_structs[0])
            key = subsol_0.get_equality_key(filling, common_dep_positions)
            found_equal = False
            for opp_filling in subsol_1.iter_fillings():
                opp_key = (
                    subsol_1.get_equality_key(
                        opp_filling, common_dep_positions))
                if common_dep_positions:
                    if key == opp_key:
                        found_equal = True
                        opp_filling.get_seq_list(seq_list, target_structs[1])
                        seq = ''.join(seq_list)
                        if seq not in check_doubles:
                            yield seq_list
                            check_doubles.add(seq)
                else:
                    opp_filling.get_seq_list(seq_list, target_structs[1])
                    yield seq_list
            if not found_equal:
                pass
            #     print filling, filling.bp_id, filling.unopt
            #     raise ValueError("No matching filling in opposite subsolution. "
            #                      "Something went somewhere terribly wrong.")


def fill_missing_bases(seq_lists, seq_constraint, unspec_b_positions):
    def rand_valid_base_id(pos):
        # randint includes 2nd argument into range
        b_id = random.randint(0, MAX_RAND_B_ID)
        # while not valid: try to find a valid base
        while seq_constraint[pos][b_id]:
            b_id = random.randint(0, MAX_RAND_B_ID)
        return b_id

    for seq_list in seq_lists:
        for pos in unspec_b_positions:
            seq_list[pos] = rna.BASES[rand_valid_base_id(pos)]

        yield rna.Sequence(seq_list)


def generate_sequences(target_structs, seq_constraint=None, local_refinement=True):
    target_structs = (
        structure.Structure(target_structs[0]),
        structure.Structure(target_structs[1]))
    if seq_constraint is None:
        seq_constraint = (
            calculate_sequence_constraint(
                target_structs,
                rna.IUPACSequence('N' * len(target_structs[0]))))
    else:
        seq_constraint = (
            calculate_sequence_constraint(target_structs,
                                          rna.IUPACSequence(seq_constraint)))

    subsolutions = two_target_inverse_fold(target_structs, seq_constraint)
    seq_lists = (
        iter_combined_seq_lists(subsolutions, target_structs, seq_constraint))
    unspec_b_positions = (
        subsolutions[0].dep_b_positions.union(
            subsolutions[1].dep_b_positions).symmetric_difference(
                xrange(len(seq_constraint))))
    sequences = fill_missing_bases(seq_lists, seq_constraint, unspec_b_positions)

    if local_refinement:
        for seq in sequences:
            # local refinement: yield seq, cost and steps
            yield loc_ref.local_search(seq, target_structs, seq_constraint)
    else:
        for seq in sequences:
            yield seq, None, None
