#!/usr/bin/env python

import unittest

from riboconstruct import helper
from riboconstruct import rna

from riboconstruct.inverse_folding import inverse_fold
from riboconstruct.inverse_folding import structure


class TestExtendedStructure(unittest.TestCase):
    def setUp(self):
        self.s = (
            structure.Structure(
                "((.(((...))).(((....))).((((...)).)).))..(.(...).)"))

    def test_trivial(self):
        struct = "(....)"
        struct = structure.Structure(struct)
        self.assertEqual("(....)", str(struct))
        struct = structure.Structure(struct)
        self.assertEqual("(....)", str(struct))

    def test_predecessor_pos_ids(self):
        self.assertEqual(
            self.s.predecessor_pos_ids,
            ((None,), (0,), (None,), (2,), (3,), (4,), (None,), (6,), (7,),
             (None,), (9,), (10,), (5, 8, 11), (12,), (1, 13)))


class TestSubstructure(unittest.TestCase):
    def setUp(self):
        structs = []
        structs.append(
            structure.Structure(
                "...((..(((...))).(((....))).((((...)).)).))..(...(...).)."))
        structs.append(
            structure.Structure("((...)(...)...(...))"))

        self.substructs = []
        self.substructs.append([])
        for substruct in structs[0].iter_substructs():
            self.substructs[0].append(substruct)
        self.substructs.append([])
        for substruct in structs[1].iter_substructs():
            self.substructs[1].append(substruct)

        self.types = []
        self.types.append(
            (rna.StructType.HAIRPIN, rna.StructType.INTERIOR,

             rna.StructType.HAIRPIN, rna.StructType.STACKING,
             rna.StructType.BULGE, rna.StructType.STACKING,

             rna.StructType.HAIRPIN, rna.StructType.STACKING,
             rna.StructType.STACKING,

             rna.StructType.HAIRPIN, rna.StructType.STACKING,
             rna.StructType.STACKING,

             rna.StructType.MULTILOOP, rna.StructType.STACKING,

             rna.StructType.EXTERIOR))
        self.types.append(
            (rna.StructType.HAIRPIN, rna.StructType.HAIRPIN,
             rna.StructType.HAIRPIN,

             rna.StructType.MULTILOOP,

             rna.StructType.EXTERIOR))

        self.num_free_bases = []
        self.num_free_bases.append(
            (3, (3, 1), 3, 0, (0, 1), 0, 4, 0, 0, 3, 0, 0, (1, 1, 1, 2),
             0, (1, 2, 3)))
        self.num_free_bases.append((3, 3, 3, (0, 3, 0, 0), (0, 0)))

        self.stem_ends = []
        self.stem_ends.append((5, 8, 11, 1, 13))
        self.stem_ends.append((0, 1, 2, 3))

    def test_structure_types(self):
        for i in xrange(len(self.substructs)):
            for substruct, substruct_type in zip(self.substructs[i],
                                                 self.types[i]):
                self.assertEqual(substruct.struct_type, substruct_type)

    def test_structure_num_free_bases(self):
        for i in xrange(len(self.substructs)):
            for substruct, num_free_bases in zip(self.substructs[i],
                                                 self.num_free_bases[i]):
                self.assertEqual(substruct.num_free_bases, num_free_bases)

    def test_structure_is_stem_end(self):
        for i in xrange(len(self.substructs)):
            for bp_pos_id, substruct in enumerate(self.substructs[i]):
                if bp_pos_id in self.stem_ends[i]:
                    self.assertTrue(substruct.is_stem_end)
                else:
                    self.assertFalse(substruct.is_stem_end)


class TestSequenceConstraint(unittest.TestCase):
    def test_constraint_1(self):
        struct = rna.Structure("(...)")
        iupac_seq = "MNNNN"
        seq_constraint = (
            inverse_fold.calculate_sequence_constraint(struct, iupac_seq))
        self.assertTupleEqual((False, False, True, True), seq_constraint[0])
        self.assertTupleEqual((True, True, False, False), seq_constraint[-1])

    def test_constraint_2(self):
        struct = rna.Structure("(...)")
        iupac_seq = "GNNNU"
        seq_constraint = (
            inverse_fold.calculate_sequence_constraint(struct, iupac_seq))
        self.assertTupleEqual((True, True, False, True), seq_constraint[0])
        self.assertTupleEqual((True, True, True, False), seq_constraint[-1])

    def test_constraint_3(self):
        struct = rna.Structure("(...)")
        iupac_seq = "RNNNC"
        seq_constraint = (
            inverse_fold.calculate_sequence_constraint(struct, iupac_seq))
        self.assertTupleEqual((True, True, False, True), seq_constraint[0])
        self.assertTupleEqual((True, False, True, True), seq_constraint[-1])


class TestEnergyCalculations(unittest.TestCase):
    def setUp(self):
        # energy values checked against INFO-RNA-2.1.2
        self.struct_energy_list = (
            # hairpin
            (("(...)", 5.7),
             ("(....)", 0.2),
             ("(.....)", 2.7),
             # stacking
             ("((...))", 2.3),
             ("(((....)))", -6.4),
             (".(((.....)))..", -6.0),
             # bulge
             ("(.(...))", 6.1),
             ("((....)..)", 3.0),
             ("((...)........)", 10.4),
             # interior (bug in INFO-RNA-2.1.2 for sizes>2)
             ("(.(...).)", 3.6),
             ("(..(...).)", 6.1),
             ("(.(...)..)", 6.1),
             ("(..(...)..)", 0.8),
             # exterior loop
             ("(...).", 4.0),
             (".(...)", 5.2),
             ("(....).", -1.5),
             ("(.....).", 1.0),
             (".(.....)...", 0.7),
             ("(...)(...)(...)", 17.1),
             ("(...)(...)(...).", 15.4),
             ("(...)(...)(...)..", 15.4),
             ("(...)(...).(...)", 15.4),
             ("(...)(...)..(...)", 14.9),
             ("(...)(...)...(...)", 14.9),
             (".(...)(...)(...)", 16.6),
             ("..(...)(...)(...)", 16.6),
             # multiloop
             ("((...)(...)(...))", 22.1),
             ("((...)(...)(...).)", 20.4),
             ("((...)(...)(...)..)", 19.9),
             ("((...)(...)(...)...)", 19.9),
             ("((...)(...).(...))", 20.4),
             ("((...)(...)..(...))", 19.9),
             ("((...)(...)...(...))", 19.9),
             ("(.(...)(...)(...))", 20.4),
             ("(..(...)(...)(...))", 19.9),
             ("(...(...)(...)(...))", 19.9),
             ("(....(...).(...)..(...))", 16.2),
             # longer structures
             ("...((..(((...))).(((....))).((((...)).)).))...(...(...).).",
              -11.2)
             ))

    def test_energy_calculations(self):
        for struct, energy in self.struct_energy_list:
            struct = structure.Structure(struct)
            seq_constraint = (
                inverse_fold.calculate_sequence_constraint(
                    struct, 'N' * len(struct)))
            subsol = inverse_fold.inverse_fold(struct, seq_constraint)
            min_energy = helper.MAX_FLOAT
            for bp_id in xrange(rna.BasepairId.count):
                for filling in subsol.iter_bp_id_fillings(bp_id):
                    if filling.energy < min_energy:
                        min_energy = filling.energy
            self.assertTrue(helper.float_equal(min_energy, energy))


if __name__ == "__main__":
    unittest.main()
