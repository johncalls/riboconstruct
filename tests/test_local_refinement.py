#!/usr/bin/env python

import unittest

from riboconstruct import rna

from riboconstruct.inverse_folding import structure
from riboconstruct.inverse_folding import two_target_local_refinement as local_refinement
from riboconstruct.inverse_folding import two_target_inverse_fold


class TestDependencyGraph(unittest.TestCase):
    def setUp(self):
        self.struct_pairs = []
        self.struct_pairs.append(
            (structure.Structure("..(...)"),
             structure.Structure("(.....)")))
        self.struct_pairs.append(
            (structure.Structure("(...(...)...)"),
             structure.Structure("(...)...(...)")))
        self.struct_pairs.append(
            (structure.Structure("(...(...)....)"),
             structure.Structure("(...)....(...)")))
        self.struct_pairs.append(
            (structure.Structure("(...(...).....)"),
             structure.Structure("(...).....(...)")))
        self.struct_pairs.append(
            (structure.Structure("((.((...)).((...))))"),
             structure.Structure(".(.((.(((...))).))).")))
        self.struct_pairs.append(
            (structure.Structure("((((...))))"),
             structure.Structure("((((...))))")))

        self.dep_graphs = []
        for structs in self.struct_pairs:
            dep_graph = local_refinement.generate_dependency_subgraphs(structs)
            self.dep_graphs.append(dep_graph)

    def test_graph_generation(self):
        expected_dep_graphs = []
        expected_dep_graphs.append(((0, 2, 6),))
        expected_dep_graphs.append(((8, 0, 12, 4),))
        expected_dep_graphs.append(((8, 0, 4, 13, 9),))
        expected_dep_graphs.append(((8, 0, 10, 4, 14),))
        expected_dep_graphs.append(
            ((0, 19), (1, 18), (11, 9, 3, 17), (8, 12, 4, 16), (6, 14),
             (13, 7)))
        expected_dep_graphs.append(((0, 10), (9, 1), (8, 2), (3, 7)))

        for dep_graph, exp_dep_graph in zip(self.dep_graphs,
                                            expected_dep_graphs):
            for subgraph, exp_subgraph in zip(dep_graph, exp_dep_graph):
                self.assertTupleEqual(
                    tuple(subgraph.iter_nodes()), exp_subgraph)

    def test_mutated_bp_positions(self):
        pass
        # TODO: implement test


class TestLocalRefinementSimple(unittest.TestCase):
    def setUp(self):
        self.structs = (
            structure.Structure("((((...))))"),
            structure.Structure("((((...))))"))

        local_refinement.start_seq = rna.Sequence("AGGGAAACCCU")
        local_refinement.target_structs = self.structs
        local_refinement.seq_constraint = (
            two_target_inverse_fold.calculate_sequence_constraint(
                self.structs,
                rna.BASES[rna.BaseId.UNSPEC] * len(self.structs[0])))

        local_refinement.OPTIMIZATION_CRITERIA = (
            local_refinement.OptimizationCriteriaId.min_pr_defect)
        local_refinement.MAX_STEPS = 1

    def test_mutation_pos_list(self):
        dep_b_positions = list()
        subgraphs = local_refinement.generate_dependency_subgraphs(self.structs)
        for subgraph in subgraphs:
            dep_b_positions.extend(subgraph.iter_nodes())
        self.assertListEqual(list([0, 10, 9, 1, 8, 2, 3, 7]), dep_b_positions)

        start_seq = local_refinement.start_seq
        indep_b_positions = list(
            set(xrange(len(start_seq))).difference(dep_b_positions))
        self.assertListEqual([4, 5, 6], indep_b_positions)

        mutation_positions_list = []
        subgraphs = local_refinement.generate_dependency_subgraphs(self.structs)
        mutation_positions_list.extend(subgraphs)
        mutation_positions_list.extend(indep_b_positions)

        all_mutation_positions = []
        for mutation_positions in mutation_positions_list:
            if type(mutation_positions) == int:
                all_mutation_positions.append(mutation_positions)
            else:
                all_mutation_positions.extend(
                    mutation_positions.iter_nodes())
        all_mutation_positions.sort()
        self.assertListEqual(
            list(range(len(start_seq))), all_mutation_positions)

    def test_fls_pf(self):
        seq_start = local_refinement.start_seq
        seq_new, _, _ = local_refinement.local_search_fls_pf()
        self.assertEqual("CGGGAAACCCG", str(seq_new))

        cost_start = local_refinement.get_cost(seq_start)
        cost_new = local_refinement.get_cost(seq_new)
        self.assertAlmostEqual(3.088, cost_start, 3)
        self.assertAlmostEqual(0.845, cost_new, 3)


if __name__ == "__main__":
    unittest.main()
