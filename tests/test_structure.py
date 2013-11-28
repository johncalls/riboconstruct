#!/usr/bin/env python

import unittest

from riboconstruct.inverse_folding import structure


class TestStructure(unittest.TestCase):
    def setUp(self):
        self.structs = (
            structure.Structure('..((...))'),
            structure.Structure('(((...)..((...).).(...)..))...(...)'))

    def test_predecessors(self):
        struct = self.structs[0]
        self.assertTupleEqual(
            ((None,), (0,), (1,)), struct.predecessor_pos_ids)
        struct = self.structs[1]
        self.assertTupleEqual(
            ((None,), (None,), (None,), (2,), (None,), (1, 3, 4), (5,), (0, 6)),
            struct.predecessor_pos_ids)

    def test_successors(self):
        struct = self.structs[0]
        self.assertTupleEqual(
            (1, 2, None), struct.successor_pos_id)
        struct = self.structs[1]
        self.assertTupleEqual(
            (7, 5, 3, 5, 5, 6, 7, None), struct.successor_pos_id)

    def test_el_stem_check(self):
        struct = self.structs[0]
        self.assertTrue(struct.successor_pos_id[1] == len(struct.bp_positions))
        struct = self.structs[1]
        self.assertTrue(struct.successor_pos_id[0] == len(struct.bp_positions))
        self.assertTrue(struct.successor_pos_id[6] == len(struct.bp_positions))

    def test_pos_mapping(self):
        struct = self.structs[0]
        expected_bp_positions = (None, None, 1, 0, None, None, None, 0, 1)
        self.assertTupleEqual(expected_bp_positions, struct.pos_2_bp_pos_id)
        struct = self.structs[1]
        expected_bp_positions = (
            6, 5,
            4, None, None, None, 4,
            None, None,
            3, 2, None, None, None, 2, None, 3,
            None,
            1, None, None, None, 1,
            None, None, 5, 6,
            None, None, None,
            0, None, None, None, 0)
        self.assertTupleEqual(expected_bp_positions, struct.pos_2_bp_pos_id)

    def test_stem_end_check(self):
        struct = self.structs[0]
        expected_stem_ends = (False, True, False)
        for bp_pos_id in xrange(len(struct.bp_positions) + 1):
            self.assertEqual(
                expected_stem_ends[bp_pos_id], struct.is_stem_end(bp_pos_id))
        struct = self.structs[1]
        expected_stem_ends = (True, True, False, True, True, False, True, False)
        for bp_pos_id in xrange(len(struct.bp_positions) + 1):
            self.assertEqual(
                expected_stem_ends[bp_pos_id], struct.is_stem_end(bp_pos_id))

    def test_substruct_iteration(self):
        struct = self.structs[0]
        expected_substructs = (
            structure.HairpinSubstructure,
            structure.StackSubstructure,
            structure.FinalSubstructure)
        for i, substruct in enumerate(struct.iter_substructs()):
            self.assertEqual(expected_substructs[i], type(substruct))
        struct = self.structs[1]
        expected_substructs = (
            structure.HairpinSubstructure,
            structure.HairpinSubstructure,
            structure.HairpinSubstructure,
            structure.BulgeSubstructure,
            structure.HairpinSubstructure,
            structure.MultiloopSubstructure,
            structure.StackSubstructure,
            structure.FinalSubstructure)
        for i, substruct in enumerate(struct.iter_substructs()):
            self.assertEqual(expected_substructs[i], type(substruct))


if __name__ == '__main__':
    unittest.main()
