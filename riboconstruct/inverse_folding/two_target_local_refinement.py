#!/usr/bin/env python

import bisect
import math
import multiprocessing as mp
import time
import random

import RNA as vienna_rna

from . import dependency_graph
from .fillings import bulge
from .fillings import exterior
from .fillings import hairpin
from .fillings import interior
from .fillings import multiloop
from .fillings import stacking
from .fillings import energies
from .. import helper
from .. import rna

class EvalSeqContainer(object):
    def __init__(self):
        self.container = {}

    def add(self, seq):
        def bp_pr(i, j):
            if i > j:
                i, j = j, i
            return bppm[iindx[i + 1] - (j + 1)]

        def pos_bp_pr_iter():
            """Yield probability for base i being paired in each structure."""
            struct_0, struct_1 = target_structs
            for i in xrange(len(seq)):
                a = struct_0.basepairs[i]
                b = struct_1.basepairs[i]
                if a is not None:
                    if b is not None:
                        yield bp_pr(i, a), bp_pr(i, b)
                    else:
                        sum_ = sum((bp_pr(i, j) for j in xrange(size)))
                        yield bp_pr(i, a), 1 - sum_
                elif b is not None:
                    sum_ = sum((bp_pr(i, j) for j in xrange(size)))
                    yield 1 - sum_, bp_pr(i, b)
                else:
                    sum_ = 1 - sum((bp_pr(i, j) for j in xrange(size)))
                    yield sum_, sum_

        seq_str = str(seq)
        if seq_str in self.container:
            return
        size = len(seq)
        self.container[seq_str] = {}

        # with LOCK_VIENNA_RNA:
        energy_ensemble = vienna_rna.pf_fold_par(seq_str, None, None, 1, 0,
                                                 0)
        bppm = vienna_rna.doubleArray_frompointer(vienna_rna.export_bppm())
        i = vienna_rna.get_iindx(len(seq_str))
        iindx = vienna_rna.intArray_frompointer(i)
        energy_structs = (
            vienna_rna.energy_of_struct(seq_str, str(target_structs[0])),
            vienna_rna.energy_of_struct(seq_str, str(target_structs[1])))
        struct_mfe = " " * len(seq_str)
        energy_mfe = vienna_rna.fold_par(seq_str, struct_mfe, None, 0, 0)
        self.container[seq_str]["pos_bp_pr"] = list(pos_bp_pr_iter())
        self.container[seq_str]["energy_ensemble"] = energy_ensemble
        self.container[seq_str]["energy_mfe"] = energy_mfe
        self.container[seq_str]["struct_mfe"] = struct_mfe
        self.container[seq_str]["seq_pr"] = (
            math.exp((energy_ensemble - energy_structs[0]) / kT),
            math.exp((energy_ensemble - energy_structs[1]) / kT))

    def get(self, seq):
        seq_str = str(seq)
        try:
            return self.container[seq_str]
        except KeyError:
            self.add(seq_str)
            return self.container[seq_str]

    def reset(self):
        self.container = {}


### Optimization criterias

def optimization_criterias(seq):
    size = len(seq)
    eval_container = eval_seq_container.get(seq)

    def structs_neg_logs_iter():
        pos_bp_pr = eval_container["pos_bp_pr"]
        for i, (pr_a, pr_b) in enumerate(pos_bp_pr):
            try:
                pr_a = -math.log(pr_a) / 10.0
            except ValueError:
                pr_a = 1.0
            try:
                pr_b = -math.log(pr_b) / 10.0
            except ValueError:
                pr_b = 1.0
            yield pr_a, pr_b

    def min_pr_defect():
        pos_bp_pr = eval_container["pos_bp_pr"]
        return size - sum(pr_a * pr_b for pr_a, pr_b in pos_bp_pr)

    def min_log_pr_defect():
        return sum(min(log_a * log_b * 4.0, 1.0)
                   for log_a, log_b in structs_neg_logs_iter())

    def min_avg_pr_defect():
        pos_bp_pr = eval_container["pos_bp_pr"]
        return size - sum((pr_a + pr_b) / 2.0
                          for pr_a, pr_b in pos_bp_pr)

    def min_avg_log_pr_defect():
        return sum(min((log_a + log_b) / 2.0, 1.0)
                   for log_a, log_b in structs_neg_logs_iter())

    def min_worst_pr_defect():
        pos_bp_pr = eval_container["pos_bp_pr"]
        return size - sum(min(pr_a, pr_b) for pr_a, pr_b in pos_bp_pr)

    def min_worst_log_pr_defect():
        return sum(min(log_a, log_b, 1.0)
                   for log_a, log_b in structs_neg_logs_iter())

    def struct_prob_variance():
        seq_pr = eval_container["seq_pr"]
        x_0 = -math.log(seq_pr[0])
        x_1 = -math.log(seq_pr[1])
        x_avg = (x_0 + x_1) / 2.
        return x_avg + VAR_ALPHA * (x_0 ** 2 + x_1 ** 2 - 2 * x_avg ** 2)

    return (min_pr_defect, min_log_pr_defect,
            min_avg_pr_defect, min_avg_log_pr_defect,
            min_worst_pr_defect, min_worst_log_pr_defect,
            struct_prob_variance)


### Local variables and parameters

class SearchStrategy:
    num = 3

    adaptive_walk, full_local_search, stochastic_local_search = xrange(num)


class NeighbourSelection:
    num = 3

    random, energy_dependent, propability_dependent = xrange(num)


class OptimizationCriteriaId:
    num = 7

    (min_pr_defect, min_log_pr_defect,
     min_avg_pr_defect, min_avg_log_pr_defect,
     min_worst_pr_defect, min_worst_log_pr_defect,
     struct_prob_variance) = xrange(num)


SLS_ACCEPT_PROB = 0.025
# TODO: evaluate
MAX_SEARCH_TIME = 60
MAX_STEPS = 3600
GASCONST = 1.98717
VAR_ALPHA = 0.75 / 2.
TEMPERATURE = 37.
K0 = 273.15
kT = (TEMPERATURE + K0) * GASCONST / 1000.

SEARCH_STRATEGY = SearchStrategy.stochastic_local_search
NEIGHBOUR_SELECTION = NeighbourSelection.propability_dependent
OPTIMIZATION_CRITERIA = OptimizationCriteriaId.struct_prob_variance

BASES = [(b_id, b) for b_id, b in enumerate(rna.BASES[:-1])]
BASEPAIRS = [(bp_id, bp) for bp_id, bp in enumerate(rna.BASEPAIRS)]

start_seq = None
seq_constraint = None
target_structs = (None, None)
eval_seq_container = EvalSeqContainer()

LOCK_VIENNA_RNA = mp.Lock()


################################################################################
###                                   Main                                   ###
################################################################################

def local_search(start_seq_, target_structs_, seq_constraint_,
                 context_front=None, context_back=None):
    global start_seq
    global seq_constraint
    global target_structs

    rna.check_struct_seq_match(target_structs_[0], start_seq_)
    rna.check_struct_seq_match(target_structs_[1], start_seq_)

    start_seq = start_seq_
    target_structs = target_structs_
    seq_constraint = seq_constraint_

    # # TODO: has to be checked
    # preset_dangles = RNA.dangles
    # if preset_dangles != 0:
    #     RNA.dangles = 1

    if (SEARCH_STRATEGY == SearchStrategy.adaptive_walk or
        SEARCH_STRATEGY == SearchStrategy.stochastic_local_search):
        seq, cost, steps = local_search_sls_pf()
    elif SEARCH_STRATEGY == SearchStrategy.full_local_search:
        seq, cost, steps = local_search_fls_pf()
    else:
        raise ValueError("Specified search strategy not valid.")

    eval_seq_container.reset()
    vienna_rna.free_pf_arrays()
    vienna_rna.free_arrays()

    # RNA.dangles = preset_dangles
    return seq, cost, steps


################################################################################
###                    The different local search methods                    ###
################################################################################

def get_cost(seq):
    eval_seq_container.add(seq)
    criterias = optimization_criterias(seq)
    return criterias[OPTIMIZATION_CRITERIA]()


# TODO: introduce stop condition based on a cost threshold

# search_strategy={aw, sls}, fold_type=pf
def local_search_sls_pf():
    using_sls = SEARCH_STRATEGY == SearchStrategy.stochastic_local_search

    def evaluate(seq):
        cost = get_cost(seq)
        rnd = random.random()
        if cost < local_search_sls_pf.cost:
            local_search_sls_pf.seq = rna.Sequence(seq)
            local_search_sls_pf.cost = cost
            local_search_sls_pf.returned_with_break = True
            return True
        if using_sls and rnd < SLS_ACCEPT_PROB:
            local_search_sls_pf.seq = rna.Sequence(seq)
            local_search_sls_pf.cost = cost
            local_search_sls_pf.returned_with_break = True
            return True
        return False

    # dependent bases (bases appearing in a bp in at least one struct)
    subgraphs = generate_dependency_subgraphs(target_structs)
    dep_b_positions = list()
    for subgraph in subgraphs:
        dep_b_positions.extend(subgraph.iter_nodes())
    # independent bases (bases unpaired in both structs)
    indep_b_positions = set(xrange(len(start_seq))).difference(dep_b_positions)
    # pool of bases/subgraphs to be rec_pos_mutationd
    mutation_positions_list = list()
    mutation_positions_list.extend(subgraphs)
    mutation_positions_list.extend(indep_b_positions)

    best_seq = start_seq
    best_cost = get_cost(best_seq)

    current_seq = best_seq
    current_cost = best_cost

    steps = 0
    start_time = time.time()
    while (time.time() - start_time < MAX_SEARCH_TIME and
           steps < MAX_STEPS):
        if NEIGHBOUR_SELECTION == NeighbourSelection.random:
            # random mutation positions order
            random.shuffle(mutation_positions_list)
        elif NEIGHBOUR_SELECTION == NeighbourSelection.propability_dependent:
            mutation_positions_list = [item for weight, item in
                                       get_unfitness(mutation_positions_list,
                                                     current_seq)]
        elif (NEIGHBOUR_SELECTION == NeighbourSelection.energy_dependent):
            # mutation positions ordered by average energy difference
            mutation_positions_list = (
                sorted_by_avg_energy_diff(mutation_positions_list, current_seq))
        else:
            raise ValueError("Specified neighour selection strategy not valid.")

        cont = False

        # NOTE: accept the first mutation that is better than the old one
        # (for a small prob. also accept mutations worse than the old one)
        for mutation_positions in mutation_positions_list:
            local_search_sls_pf.seq = current_seq
            local_search_sls_pf.cost = current_cost
            local_search_sls_pf.returned_with_break = False

            mutated_seq = rna.Sequence(current_seq)

            # check if mutation_positions is a single base position..
            if type(mutation_positions) == int:
                b_pos = mutation_positions
                random.shuffle(BASES)
                for b_id, b in BASES:
                    if b == mutated_seq[b_pos]:
                        continue
                    if seq_constraint[b_pos][b_id]:
                        continue
                    # rec_pos_mutation
                    mutated_seq[b_pos] = b
                    # evaluate and break if a better sequence has been found
                    # (with a certain probability worse solutions are accepted)
                    if evaluate(mutated_seq):
                        break

            # ..or a list of dependent base positions part of basepairs in at
            # least one of the struct
            else:
                # mutate the positions
                for _ in mutate_seq(mutated_seq, mutation_positions, True):
                    # evaluate and break if a better sequence has been found
                    # (with a certain probability worse solutions are accepted)
                    if evaluate(mutated_seq):
                        break

            if local_search_sls_pf.returned_with_break:
                # check if the found seq is the currently best one
                if local_search_sls_pf.cost + helper.CMP_THR < best_cost:
                    best_seq = local_search_sls_pf.seq
                    best_cost = local_search_sls_pf.cost
                current_seq = local_search_sls_pf.seq
                current_cost = local_search_sls_pf.cost
                cont = True
                break

            steps += 1

        # if the sequence was not altered: break
        if not cont:
        # if id(best_seq) == id(current_seq):
            break
        # if the sequence was improved: set new best sequence
        if current_cost + helper.CMP_THR < best_cost:
            best_seq = current_seq
            best_cost = current_cost

    return best_seq, best_cost, steps


# search_strategy=fls, fold_type=pf
def local_search_fls_pf():
    def evaluate(seq):
        cost = get_cost(seq)
        if cost < local_search_fls_pf.cost:
            local_search_fls_pf.seq = rna.Sequence(seq)
            local_search_fls_pf.cost = cost

    # dependent bases (bases appearing in a bp in at least one struct)
    subgraphs = generate_dependency_subgraphs(target_structs)
    dep_b_positions = list()
    for subgraph in subgraphs:
        dep_b_positions.extend(subgraph.iter_nodes())
    # independent bases (bases unpaired in both structs)
    indep_b_positions = set(xrange(len(start_seq))).difference(dep_b_positions)
    # pool of bases/subgraphs to be rec_pos_mutationd
    mutation_positions_list = list()
    mutation_positions_list.extend(subgraphs)
    mutation_positions_list.extend(indep_b_positions)

    best_seq = start_seq
    best_cost = get_cost(best_seq)

    current_seq = best_seq
    current_cost = best_cost

    steps = 0
    start_time = time.time()
    while (time.time() - start_time < MAX_SEARCH_TIME and
           steps < MAX_STEPS):
        # randomize the base/subgraph pool
        random.shuffle(mutation_positions_list)

        for mutation_positions in mutation_positions_list:
            local_search_fls_pf.seq = best_seq
            local_search_fls_pf.cost = best_cost

            mutated_seq = rna.Sequence(best_seq)

            # check if mutation_positions is a single base position..
            if type(mutation_positions) == int:
                b_pos = mutation_positions
                random.shuffle(BASES)
                for b_id, b in BASES:
                    if b == mutated_seq[b_pos]:
                        continue
                    if seq_constraint[b_pos][b_id]:
                        continue
                    # mutate
                    mutated_seq[b_pos] = b
                    # evaluate
                    evaluate(mutated_seq)

            # ..or a list of dependent base positions part of basepairs in at
            # least one of the struct
            else:
                # mutate the positions
                for _ in mutate_seq(mutated_seq, mutation_positions, True):
                    # evaluate
                    evaluate(mutated_seq)

            if local_search_fls_pf.cost + helper.CMP_THR < current_cost:
                current_seq = local_search_fls_pf.seq
                current_cost = local_search_fls_pf.cost

            steps += 1

        # if the sequence was not improved: break
        if id(best_seq) == id(current_seq):
            break
        # otherwise: set new best sequence
        else:
            best_seq = current_seq
            best_cost = current_cost

    return best_seq, best_cost, steps


################################################################################
###   Different helper functions and classes used during the local search    ###
################################################################################

def generate_dependency_subgraphs(structs):
    dep_graph = dependency_graph.DependencyGraph()
    struct_0, struct_1 = structs
    for bp_pos in struct_0.bp_positions:
        dep_graph.add_edge(bp_pos)
    for bp_pos in struct_1.bp_positions:
        dep_graph.add_edge(bp_pos)
    return list(dep_graph.iter_subgraphs())


def mutate_seq(seq, dep_subgraph, rand=True):
    def b_base_matching(bp_b):
        if bp_b == 'A':
            return [(rna.BaseId.U, 'U')]
        elif bp_b == 'C':
            return [(rna.BaseId.G, 'G')]
        elif bp_b == 'G':
            return [(rna.BaseId.C, 'C'), (rna.BaseId.U, 'U')]
        elif bp_b == 'U':
            return [(rna.BaseId.A, 'A'), (rna.BaseId.C, 'G')]
        else:
            return list()

    def rec_pos_mutation(b_index_next, bases_next):
        b_index = b_index_next
        bases = bases_next
        b_index_next += 1
        b_pos = mutation_positions[b_index]
        if rand:
            random.shuffle(bases)
        # check stop criterion
        if b_index_next == len(mutation_positions):
            for b_id, b in bases:
                if seq_constraint[b_pos][b_id]:
                    continue
                if dep_subgraph.is_cyclic and seq[b_pos] != b:
                    continue
                seq[b_pos] = b
                yield seq
            return
        for b_id, b in bases:
            if seq_constraint[b_pos][b_id]:
                continue
            seq[b_pos] = b
            for _ in rec_pos_mutation(b_index_next, b_base_matching(b)):
                yield seq

    if dep_subgraph.is_cyclic:
        mutation_positions = tuple(dep_subgraph.get_rand_nodes_iter())
        b_pos = mutation_positions[0]
        if rand:
            random.shuffle(BASES)
        for b_id, b in BASES:
            if seq_constraint[b_pos][b_id]:
                continue
            seq[b_pos] = b
            for _ in rec_pos_mutation(1, b_base_matching(b)):
                yield seq
    else:
        mutation_pos_iterators = dep_subgraph.get_rand_nodes_iter()
        try:
            iter_0, iter_1 = mutation_pos_iterators
            mutation_positions_0 = tuple(iter_0)
            mutation_positions_1 = tuple(iter_1)
            single = False
        except ValueError:
            mutation_positions_0 = tuple(mutation_pos_iterators[0])
            single = True

        mutation_positions = mutation_positions_0
        b_pos = mutation_positions[0]

        if rand:
            random.shuffle(BASES)
        for b_id, b in BASES:
            if seq_constraint[b_pos][b_id]:
                continue
            seq[b_pos] = b
            for _ in rec_pos_mutation(1, b_base_matching(b)):
                if single:
                    yield seq
                else:
                    mutation_positions = mutation_positions_1
                    b = seq[mutation_positions[0]]
                    b_index_next = 1
                    for _ in rec_pos_mutation(b_index_next, b_base_matching(b)):
                        yield seq
                    mutation_positions = mutation_positions_0


def get_unfitness(mutation_positions_list, seq):
    unfitness = []
    size = len(seq)
    eval_seq_container.add(seq)
    eval_seq = eval_seq_container.get(seq)

    seq_unfitness = size - sum(bp_pr_0 * bp_pr_1
                               for bp_pr_0, bp_pr_1 in eval_seq["pos_bp_pr"])

    for mutation_positions in mutation_positions_list:
        if type(mutation_positions) == int:
            pos_bp_pr = eval_seq["pos_bp_pr"][mutation_positions]
            unfitness.append(
                ((1. - pos_bp_pr[0] * pos_bp_pr[1]) / seq_unfitness,
                 mutation_positions))
        else:
            strech_prob = 0.
            for pos in mutation_positions.iter_nodes():
                pos_bp_pr = eval_seq["pos_bp_pr"][pos]
                strech_prob += 1. - pos_bp_pr[0] * pos_bp_pr[1]
            unfitness.append((strech_prob / seq_unfitness, mutation_positions))

    return sorted(unfitness, key=lambda f: f[0], reverse=True)


def get_unfitness_prob_selector(mutation_positions_list, seq):
    def choice(rnd=random.random, bis=bisect.bisect):
        return items[bis(added_unfitness, rnd() * last_sum)][0]

    items = get_unfitness(mutation_positions_list, seq)
    added_unfitness = []
    last_sum = 0
    for weight, item in items:
        last_sum += weight
        added_unfitness.append(last_sum)

    return choice


# TODO: not really usefule --> remove
def sorted_by_avg_energy_diff(mutation_positions_list, seq):
    def get_b_part_energy(b_pos, seq, struct):
        # if the free base is not adjacent to a stem it has no influence on the
        # energy
        b_pos_pred = struct[b_pos - 1]
        b_pos_succ = struct[b_pos + 1]
        if b_pos == 0 and b_pos_succ == '.':
            return 0.0
        elif b_pos == len(struct) - 1 and b_pos_pred == '.':
            return 0.0
        elif b_pos_pred == '.' and b_pos_succ == '.':
            return 0.0

        # finding the group of free bases
        k = b_pos
        while k >= 0 and struct[k] == '.':
            k -= 1
        group_j = k + 1
        k = b_pos
        while k < len(struct) and struct[k] == '.':
            k += 1
        group_i = k - 1

        if group_j == 0 or group_i == len(struct) - 1:
            return exterior.exterior_energy(seq, struct)

        bp_pos_pred = group_j - 1
        bp_pos_id_pred = struct.pos_2_bp_pos_id[bp_pos_pred]
        bp_pos_succ = group_i + 1
        bp_pos_id_succ = struct.pos_2_bp_pos_id[bp_pos_succ]

        # b_pos in hairpin
        if bp_pos_id_succ == bp_pos_id_pred:
            return hairpin.hairpin_energy(bp_pos_id_succ, seq, struct)
        # b_pos 'on the right side' of a bulge/il/ml
        elif struct[bp_pos_pred] == struct[bp_pos_succ] == ')':
            return get_part_energy(bp_pos_id_succ, bp_pos_id_pred, seq)
        # b_pos 'on the left side' of a bulge/il/ml
        elif struct[bp_pos_pred] == struct[bp_pos_succ] == '(':
            return get_part_energy(bp_pos_id_pred, bp_pos_id_succ, seq)
        elif (struct.successor_pos_id[bp_pos_id_pred] ==
              struct.successor_pos_id[bp_pos_id_succ]):
            # ml
            if struct.successor_pos_id[bp_pos_id_pred] != len(struct.basepairs):
                return (
                    multiloop.multiloop_energy(
                        struct.successor_pos_id[bp_pos_id_pred], seq, struct))
            # el
            else:
                return exterior.exterior_energy(seq, struct)
        # TODO: remove after testing
        else:
            raise NotImplementedError('CHECK THIS OUT')

    def get_bp_part_energy(bp_pos_id, seq, struct):
        pred_bp_pos_ids = struct.predecessor_pos_ids[bp_pos_id]
        succ_bp_pos_id = struct.successor_pos_id[bp_pos_id]

        # hl
        if pred_bp_pos_ids[0] == None:
            energy = hairpin.hairpin_energy(bp_pos_id, seq, struct)
            if succ_bp_pos_id is not None:
                return energy + get_part_energy(succ_bp_pos_id, bp_pos_id, seq)
        # el
        elif succ_bp_pos_id is None:
            energy = exterior.exterior_energy(seq, struct)
        # ml
        elif len(pred_bp_pos_ids) > 1:
            return (multiloop.multiloop_energy(bp_pos_id, seq, struct) +
                    get_part_energy(succ_bp_pos_id, bp_pos_id, seq))
        # bulge / interior / stack
        else:
            return (get_part_energy(bp_pos_id, pred_bp_pos_ids[0], seq) +
                    get_part_energy(succ_bp_pos_id, bp_pos_id, seq))

    def get_part_energy(bp_pos_id, pred_bp_pos_id, seq, struct):
        try:
            bp_pos_i, bp_pos_j = struct.bp_positions[bp_pos_id]
        except IndexError:
            bp_pos_i, bp_pos_j = 0, len(struct) - 1
        pred_bp_pos_i, pred_bp_pos_j = struct.bp_positions[pred_bp_pos_id]

        bp_i = seq[bp_pos_i]
        bp_j = seq[bp_pos_j]
        pred_bp_i = seq[pred_bp_pos_i]
        pred_bp_j = seq[pred_bp_pos_j]

        bp_id = getattr(rna.BasepairId, '%s%s' % (bp_i, bp_j))
        pred_bp_id = getattr(rna.BasepairId, '%s%s' % (pred_bp_i, pred_bp_j))

        # stacking
        if pred_bp_pos_i - bp_pos_i == 1 and bp_pos_j - pred_bp_pos_j == 1:
            # stem end penalty
            if struct.is_stem_end(bp_pos_id):
                penalty_energy = energies.au_stem_end_energy(bp_id)
            else:
                penalty_energy = 0.0
            return stacking.stacking_energy(bp_id, pred_bp_id) + penalty_energy
        # bulge (right)
        elif pred_bp_pos_i - bp_pos_i == 1 and bp_pos_j - pred_bp_pos_j > 1:
            # stem end penalty
            if struct.is_stem_end(bp_pos_id):
                penalty_energy = energies.au_stem_end_energy(bp_id)
            else:
                penalty_energy = 0.0
            return (
                bulge.bulge_energy(
                    bp_pos_j - pred_bp_pos_j, bp_id, pred_bp_id) +
                penalty_energy)
        # bulge (left)
        elif pred_bp_pos_i - bp_pos_i > 1 and bp_pos_j - pred_bp_pos_j == 1:
            # stem end penalty
            if struct.is_stem_end(bp_pos_id):
                penalty_energy = energies.au_stem_end_energy(bp_id)
            else:
                penalty_energy = 0.0
            return (
                bulge.bulge_energy(
                    pred_bp_pos_i - bp_pos_i, bp_id, pred_bp_id) +
                penalty_energy)
        # ml
        elif len(struct.predecessor_pos_ids[bp_pos_id]) > 1:
            # penalties are handled inside the function
            return multiloop.multiloop_energy(bp_pos_id, seq, struct)
        # el
        elif struct.successor_pos_id[bp_pos_id] is None:
            # penalties are handled inside the function
            return exterior.exterior_energy(seq, struct)
        # il
        else:
            # stem end penalty
            if struct.is_stem_end(bp_pos_id):
                penalty_energy = (
                    penalty_energy + energies.au_stem_end_energy(bp_id))
            else:
                penalty_energy = 0.0
            return (
                interior.interior_energy(
                    pred_bp_pos_i - bp_pos_i - 1,
                    bp_pos_j - pred_bp_pos_j - 1,
                    bp_id, pred_bp_id, bp_pos_i, bp_pos_j, seq) +
                penalty_energy)

    ### Calculate energy difference of base variants at position b_pos #########

    str_seq = str(seq)
    avg_e_diffs = []

    struct_0, struct_1 = target_structs

    for mutation_pos in mutation_positions_list:
        seq_new = rna.Sequence(seq)
        e_diffs = []

        # base at b_pos is a free base
        if type(mutation_pos) == int:
            e_0 = get_b_part_energy(mutation_pos, seq, struct_0)
            e_1 = get_b_part_energy(mutation_pos, seq, struct_1)

            b = seq[mutation_pos]
            for b_new_id, b_new in enumerate(rna.BASES[:-1]):
                if b == b_new:
                    e_diffs.append(0.0)
                    continue
                if seq_constraint[mutation_pos][b_new_id]:
                    continue
                seq_new[mutation_pos] = b_new
                e_new_0 = get_b_part_energy(mutation_pos, seq_new, struct_0)
                e_new_1 = get_b_part_energy(mutation_pos, seq_new, struct_1)
                # TODO: speed up -- how?
                # TODO: use different 'combinations' (currently only average is
                # implemented)
                e_diffs.append(
                    helper.min_sum(
                        (helper.min_sum(e_0, e_new_0),
                         helper.min_sum(e_1, e_new_1))) / 2.0)
        # b_pos is a subgraph of dependent basepairs
        else:
            mutation_pos = list(mutation_pos.iter_nodes())
            max_bp_pos_id_0 = (
                max((struct_0.pos_2_bp_pos_id[pos] for pos in mutation_pos)))
            max_bp_pos_id_1 = (
                max((struct_1.pos_2_bp_pos_id[pos] for pos in mutation_pos)))
            e_0 = get_bp_part_energy(max_bp_pos_id_0, seq, struct_0)
            e_1 = get_bp_part_energy(max_bp_pos_id_1, seq, struct_1)
            # e_0 = (
            #     RNA.energy_of_struct(
            #         str_seq, str(target_structs[0])))
            # e_1 = (
            #     RNA.energy_of_struct(
            #         str_seq, str(target_structs[1])))

            for _ in mutate_seq(seq_new, mutation_pos):
                str_seq_new = str(seq_new)
                if str_seq == str_seq_new:
                    e_diffs.append(0.0)
                    continue
                e_new_0 = get_bp_part_energy(max_bp_pos_id_0, seq_new, struct_0)
                e_new_1 = get_bp_part_energy(max_bp_pos_id_1, seq_new, struct_1)
                # e_new_0 = RNA.energy_of_struct(str_seq_new, struct_0)
                # e_new_1 = RNA.energy_of_struct(str_seq_new, struct_1)
                # TODO: speed up
                # TODO: use different 'combinations' (currently only average is
                # implemented)
                e_diffs.append(
                    helper.min_sum(
                        (helper.min_sum(e_0, e_new_0),
                         helper.min_sum(e_1, e_new_1))) / 2.0)

        avg_e_diffs.append((sum(e_diffs) / len(e_diffs), mutation_pos))

    for _, mutation_pos in sorted(avg_e_diffs,
                                  key=lambda x: x[0], reverse=True):
        yield mutation_pos
