#!/usr/bin/env python

import unittest

from riboconstruct.inverse_folding import dependency_graph as dg
from riboconstruct.inverse_folding import structure


def get_subgraphs_iter(structs):
    dep_graph = dg.DependencyGraph()
    for bp_pos in structs[0].bp_positions:
        dep_graph.add_edge(bp_pos)
    for bp_pos in structs[1].bp_positions:
        dep_graph.add_edge(bp_pos)
    return dep_graph.iter_subgraphs()


class TestDependencyGraph(unittest.TestCase):
    def test_graph_1(self):
        structs = (structure.Structure("((...)...)"),
                   structure.Structure(".....(...)"))
        # only one subgraph
        subgraph = get_subgraphs_iter(structs).next()
        self.assertTupleEqual((0, 1, 5, 9),
                              tuple(sorted(subgraph.iter_nodes())))
        self.assertTupleEqual(((0, 9), (5, 1), (9, 5)),
                              tuple(sorted(subgraph.iter_edges())))


if __name__ == "__main__":
    unittest.main()
