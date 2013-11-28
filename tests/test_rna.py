#!/usr/bin/env python

import unittest

from riboconstruct import rna


class TestStructure(unittest.TestCase):
    def setUp(self):
        self.s = (
            rna.Structure("((.(((...))).(((....))).((((...)).)).))..(.(...).)"))

    def test_trivial_1(self):
        struct = "(....)"
        struct = rna.Structure(struct)
        self.assertEqual("(....)", str(struct))
        struct = rna.Structure(struct)
        self.assertEqual("(....)", str(struct))

    def test_basepairs(self):
        basepairs = self.s.basepairs
        self.assertEqual(basepairs[0], 38)
        self.assertEqual(basepairs[1], 37)
        self.assertEqual(basepairs[2], None)
        self.assertEqual(basepairs[3], 11)
        self.assertEqual(basepairs[4], 10)
        self.assertEqual(basepairs[5], 9)
        self.assertEqual(basepairs[6], None)
        self.assertEqual(basepairs[7], None)
        self.assertEqual(basepairs[8], None)
        self.assertEqual(basepairs[9], 5)
        self.assertEqual(basepairs[10], 4)
        self.assertEqual(basepairs[11], 3)
        self.assertEqual(basepairs[12], None)
        self.assertEqual(basepairs[13], 22)
        self.assertEqual(basepairs[14], 21)
        self.assertEqual(basepairs[15], 20)
        self.assertEqual(basepairs[16], None)
        self.assertEqual(basepairs[17], None)
        self.assertEqual(basepairs[18], None)
        self.assertEqual(basepairs[19], None)
        self.assertEqual(basepairs[20], 15)
        self.assertEqual(basepairs[21], 14)
        self.assertEqual(basepairs[22], 13)
        self.assertEqual(basepairs[23], None)
        self.assertEqual(basepairs[24], 35)
        self.assertEqual(basepairs[25], 34)
        self.assertEqual(basepairs[26], 32)
        self.assertEqual(basepairs[27], 31)
        self.assertEqual(basepairs[28], None)
        self.assertEqual(basepairs[29], None)
        self.assertEqual(basepairs[30], None)
        self.assertEqual(basepairs[31], 27)
        self.assertEqual(basepairs[32], 26)
        self.assertEqual(basepairs[33], None)
        self.assertEqual(basepairs[34], 25)
        self.assertEqual(basepairs[35], 24)
        self.assertEqual(basepairs[36], None)
        self.assertEqual(basepairs[37], 1)
        self.assertEqual(basepairs[38], 0)
        self.assertEqual(basepairs[39], None)
        self.assertEqual(basepairs[40], None)
        self.assertEqual(basepairs[41], 49)
        self.assertEqual(basepairs[42], None)
        self.assertEqual(basepairs[43], 47)
        self.assertEqual(basepairs[44], None)
        self.assertEqual(basepairs[45], None)
        self.assertEqual(basepairs[46], None)
        self.assertEqual(basepairs[47], 43)
        self.assertEqual(basepairs[48], None)
        self.assertEqual(basepairs[49], 41)

    def test_bp_positions(self):
        self.assertEqual(
            self.s.bp_positions,
            ((43, 47), (41, 49), (27, 31), (26, 32), (25, 34), (24, 35),
             (15, 20), (14, 21), (13, 22), (5, 9), (4, 10), (3, 11), (1, 37),
             (0, 38)))

    def test_bad_structures(self):
        with self.assertRaises(ValueError):
            # hairpin size > 30
            rna.Structure(
                "((.(((...))).(((....)))))..(...............................)")
        with self.assertRaises(ValueError):
            # hairpin size < 3
            rna.Structure("((.(((...))).(((....)))))..(..)")
        with self.assertRaises(ValueError):
            # missing matching closing bracket
            rna.Structure("((.(((...))).(((....))))..")
        with self.assertRaises(ValueError):
            # missing matching opening bracket
            rna.Structure(".(((...))).(((....))))..")


class TestIUPACSequence(unittest.TestCase):
    def test_init(self):
        rna.IUPACSequence('ACGUNRYSWKMBDHV')

    def test_init_fails(self):
        with self.assertRaises(TypeError):
            rna.IUPACSequence('ACGUNRYSWKMBDHVL')


if __name__ == "__main__":
    unittest.main()
