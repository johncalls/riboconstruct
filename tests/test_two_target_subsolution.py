#!/usr/bin/env python

import unittest

from riboconstruct import rna

from riboconstruct.inverse_folding import structure
from riboconstruct.inverse_folding import two_target_inverse_fold as inverse_fold

from riboconstruct.inverse_folding.two_target_subsolution import bulge
from riboconstruct.inverse_folding.two_target_subsolution import exterior
from riboconstruct.inverse_folding.two_target_subsolution import filling as fl
from riboconstruct.inverse_folding.two_target_subsolution import hairpin
from riboconstruct.inverse_folding.two_target_subsolution import interior
from riboconstruct.inverse_folding.two_target_subsolution import multiloop
from riboconstruct.inverse_folding.two_target_subsolution import stacking


def get_subsolution(struct_id, substruct, task):
    def create_subsolution(predecessor=None):
        if substruct.struct_type == rna.StructType.HAIRPIN:
            return hairpin.TwoTargetHairpinSubsolution(
                struct_id, substruct, task)
        elif substruct.struct_type == rna.StructType.BULGE:
            return bulge.TwoTargetBulgeSubsolution(
                struct_id, substruct, task, predecessor)
        elif substruct.struct_type == rna.StructType.STACKING:
            return stacking.TwoTargetStackSubsolution(
                struct_id, substruct, task, predecessor)
        elif substruct.struct_type == rna.StructType.INTERIOR:
            return interior.TwoTargetInteriorSubsolution(
                struct_id, substruct, task, predecessor)
        elif substruct.struct_type == rna.StructType.MULTILOOP:
            return multiloop.TwoTargetMultiloopSubsolution(
                struct_id, substruct, task, predecessor)
        elif substruct.struct_type == rna.StructType.EXTERIOR:
            return exterior.TwoTargetFinalSubsolution(
                struct_id, substruct, task, predecessor)

    if substruct.struct_type == rna.StructType.MULTILOOP:
        pred = []
        for stem_bp_pos_id in substruct.predecessor_pos_ids:
            pred.append(
                task.subsolution_collector[struct_id][stem_bp_pos_id])
            del task.subsolution_collector[struct_id][stem_bp_pos_id]
        subsol = create_subsolution(pred)
        task.subsolution_collector[struct_id][substruct.bp_pos_id] = subsol
    elif substruct.struct_type == rna.StructType.EXTERIOR:
        pred = []
        for stem_bp_pos_id in substruct.predecessor_pos_ids:
            pred.append(
                task.subsolution_collector[struct_id][stem_bp_pos_id])
            del task.subsolution_collector[struct_id][stem_bp_pos_id]
        subsol = create_subsolution(pred)
        task.subsolution_collector[struct_id][substruct.bp_pos_id] = subsol
    elif substruct.struct_type == rna.StructType.HAIRPIN:
        subsol = create_subsolution()
        task.subsolution_collector[struct_id][substruct.bp_pos_id] = subsol
    else:
        pred = (
            task.subsolution_collector[struct_id][substruct.bp_pos_id - 1])
        del task.subsolution_collector[struct_id][substruct.bp_pos_id - 1]
        subsol = create_subsolution(pred)
        task.subsolution_collector[struct_id][substruct.bp_pos_id] = subsol

    return subsol


def get_base(subsol, filling, b_pos):
    filling_b_pos, bp_pos_id = subsol.b_pos_mapping[b_pos]
    if filling.bp_pos_id == bp_pos_id:
        return filling.base_ids[filling_b_pos]
    if (type(filling) == fl.TwoTargetExteriorFilling or
        type(filling) == fl.TwoTargetMultiloopStemFilling):
        _, bp_pos_j = (
            subsol.struct.bp_positions[filling.stem_filling.bp_pos_id])
        if b_pos <= bp_pos_j:
            return get_base(subsol, filling.stem_filling, b_pos)
        if filling.predecessor:
            pred_bp_id_i, _ = (
                subsol.struct.bp_positions[
                    filling.predecessor.stem_filling.bp_pos_id])
            if b_pos < pred_bp_id_i:
                return filling.base_ids[filling_b_pos]
            return get_base(subsol, filling.predecessor, b_pos)
        return filling.base_ids[filling_b_pos]
    return get_base(subsol, filling.predecessor, b_pos)


class TestFillings(unittest.TestCase):
    def setUp(self):
        self.target_structs = (
            structure.Structure("......(((((...)))))............."),
            structure.Structure(".(((((...))))).................."))
        self.iupac_seq = rna.IUPACSequence("CAAAACNNNGUUUUGNNNNNNNNNNNNNNNNN")
        self.seq_constraint = (
            inverse_fold.calculate_sequence_constraint(
                self.target_structs, self.iupac_seq))
        self.task = inverse_fold.Task(self.target_structs, self.seq_constraint)
        self.subsolutions = ([], [])
        for bp_belonging, substruct in self.task.iter_substructs():
            if bp_belonging == 2:
                subsol = get_subsolution(0, substruct, self.task)
                self.subsolutions[0].append(subsol)
                subsol = get_subsolution(1, substruct, self.task)
                self.subsolutions[1].append(subsol)
            else:
                subsol = get_subsolution(bp_belonging, substruct, self.task)
                if type(subsol) == exterior.TwoTargetFinalSubsolution:
                    exterior.TwoTargetExteriorStem(subsol, 0)
                self.subsolutions[bp_belonging].append(subsol)

    def test_fillings(self):
        filling = fl.TwoTargetHairpinFilling(0, 1, (), 3, 0.0)

        filling = fl.TwoTargetStackFilling(1, 0, 0.0, filling)
        bp_pos_i, bp_pos_j = (
            self.target_structs[1].bp_positions[filling.bp_pos_id])
        self.assertTupleEqual((4, 10), (bp_pos_i, bp_pos_j))
        subseq = ''.join(self.iupac_seq[bp_pos_i:bp_pos_j + 1])
        self.assertEqual(subseq, str(filling))
        filling = fl.TwoTargetStackFilling(2, 0, 0.0, filling)
        filling = fl.TwoTargetStackFilling(3, 0, 0.0, filling)

        filling = fl.TwoTargetStackFilling(4, 0, 0.0, filling)
        bp_pos_i, bp_pos_j = (
            self.target_structs[1].bp_positions[filling.bp_pos_id])
        self.assertTupleEqual((1, 13), (bp_pos_i, bp_pos_j))
        subseq = ''.join(self.iupac_seq[bp_pos_i:bp_pos_j + 1])
        self.assertEqual(subseq, str(filling))

        filling = fl.TwoTargetExteriorFilling(0, 5, (2,), 18, 0.0, filling)
        self.assertEqual(str(self.iupac_seq)[1:], str(filling))

        filling = fl.TwoTargetFinalFilling(5, 1, 1, 1, 0.0, filling)
        self.assertEqual(str(self.iupac_seq), str(filling))

    def test_fillings_2(self):
        b_u = rna.BASES[rna.BaseId.UNSPEC]

        filling = fl.TwoTargetHairpinFilling(0, 1, (b_u, b_u, b_u), 3, 0.0)
        subsol = self.subsolutions[1][filling.bp_pos_id]
        base_ids = (get_base(subsol, filling, b_pos) for b_pos in [5, 9])
        self.assertListEqual(['C', 'G'],
                             [rna.BASES[b_id] for b_id in base_ids])

        filling = fl.TwoTargetStackFilling(1, 0, 0.0, filling)
        filling = fl.TwoTargetStackFilling(2, 0, 0.0, filling)
        filling = fl.TwoTargetStackFilling(3, 0, 0.0, filling)
        filling = fl.TwoTargetStackFilling(4, 0, 0.0, filling)
        subsol = self.subsolutions[1][filling.bp_pos_id]
        base_ids = (get_base(subsol, filling, b_pos)
                    for b_pos in [1, 2, 3, 4, 5, 9, 10, 11, 12, 13])
        self.assertListEqual(['A', 'A', 'A', 'A', 'C', 'G', 'U', 'U', 'U', 'U'],
                             [rna.BASES[b_id] for b_id in base_ids])

        filling = fl.TwoTargetStackFilling(4, 0, 0.0, filling)
        subsol = self.subsolutions[1][filling.bp_pos_id]
        base_ids = (get_base(subsol, filling, b_pos)
                    for b_pos in [1, 2, 3, 4, 5, 9, 10, 11, 12, 13])
        self.assertListEqual(['A', 'A', 'A', 'A', 'C', 'G', 'U', 'U', 'U', 'U'],
                             [rna.BASES[b_id] for b_id in base_ids])

        filling = fl.TwoTargetExteriorFilling(0, 5, (2,), 18, 0.0, filling)
        subsol = self.subsolutions[1][filling.bp_pos_id]
        base_ids = (get_base(subsol, filling, b_pos)
                    for b_pos in [1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 14])
        self.assertListEqual(['A', 'A', 'A', 'A', 'C', 'G', 'U', 'U', 'U', 'U',
                              'G'],
                             [rna.BASES[b_id] for b_id in base_ids])

    def test_dependencies(self):
        b_u = rna.BASES[rna.BaseId.UNSPEC]

        filling = fl.TwoTargetHairpinFilling(0, 1, (b_u, b_u, b_u), 3, 0.0)
        filling = fl.TwoTargetStackFilling(1, 0, 0.0, filling)
        filling = fl.TwoTargetStackFilling(2, 0, 0.0, filling)
        filling = fl.TwoTargetStackFilling(3, 0, 0.0, filling)
        filling = fl.TwoTargetStackFilling(4, 0, 0.0, filling)
        filling = fl.TwoTargetExteriorFilling(0, 4, (2,), 13, 0.0, filling)

        subsol = self.subsolutions[0][-1]
        # self.assertTrue(
        #     subsol.opp_struct_b_constrained(filling, 5, rna.BaseId.A))
        # self.assertTrue(
        #     subsol.opp_struct_b_constrained(filling, 5, rna.BaseId.C))
        # self.assertTrue(
        #     subsol.opp_struct_b_constrained(filling, 5, rna.BaseId.G))
        # self.assertFalse(
        #     subsol.opp_struct_b_constrained(filling, 5, rna.BaseId.U))


class TestFillings2(unittest.TestCase):
    def setUp(self):
        self.target_structs = (
            structure.Structure(".((...))..((.(((...).))...(...)(...)).)."),
            structure.Structure(".(((((...)))))..........................."))
        self.iupac_seq = rna.IUPACSequence('N' * len(self.target_structs[0]))
        self.seq_constraint = (
            inverse_fold.calculate_sequence_constraint(
                self.target_structs, self.iupac_seq))
        self.task = inverse_fold.Task(self.target_structs, self.seq_constraint)
        self.subsolutions = ([], [])
        for bp_belonging, substruct in self.task.iter_substructs():
            if bp_belonging == 2:
                subsol = get_subsolution(0, substruct, self.task)
                self.subsolutions[0].append(subsol)
                subsol = get_subsolution(1, substruct, self.task)
                self.subsolutions[1].append(subsol)
            else:
                subsol = get_subsolution(bp_belonging, substruct, self.task)
                if bp_belonging == 0:
                    if type(subsol) == exterior.TwoTargetFinalSubsolution:
                        ext_stem_subsol = (
                            exterior.TwoTargetExteriorStem(subsol, 0))
                        exterior.TwoTargetExteriorStem(
                            subsol, 1, ext_stem_subsol)
                    if type(subsol) == multiloop.TwoTargetMultiloopSubsolution:
                        ml_stem_subsol = (
                            multiloop.TwoTargetMultiloopStem(subsol, 0))
                        ml_stem_subsol = (
                            multiloop.TwoTargetMultiloopStem(
                                subsol, 1, ml_stem_subsol))
                        multiloop.TwoTargetMultiloopStem(
                            subsol, 2, ml_stem_subsol)
                self.subsolutions[bp_belonging].append(subsol)

    def test_get_base(self):
        b_u = rna.BASES[rna.BaseId.UNSPEC]

        filling = fl.TwoTargetHairpinFilling(0, 0, (b_u, b_u, b_u), 3, 0.0)
        filling_0 = (
            fl.TwoTargetMultiloopStemFilling(0, 5, (None,), 0, 0.0, filling))
        filling = fl.TwoTargetHairpinFilling(1, 0, (b_u, b_u, b_u), 3, 0.0)
        filling_1 = (
            fl.TwoTargetMultiloopStemFilling(
                0, 5, (None,), 0, 0.0, filling, filling_0))
        filling = fl.TwoTargetHairpinFilling(2, 0, (b_u, b_u, b_u), 3, 0.0)
        filling = fl.TwoTargetBulgeFilling(3, 1, (0, 1), 0.0, filling)
        filling = fl.TwoTargetStackFilling(4, 2, 0.0, filling)
        filling_2 = (
            fl.TwoTargetMultiloopStemFilling(
                2, 5, (3, 1), 3, 0.0, filling, filling_1))
        filling = (
            fl.TwoTargetMultiloopFilling(5, 3, (3,), 1, 0.0, filling_2))
        filling = fl.TwoTargetBulgeFilling(6, 4, (0, 1), 0.0, filling)
        filling_0 = fl.TwoTargetExteriorFilling(4, 9, (1,), 1, 0.0, filling)
        filling = fl.TwoTargetHairpinFilling(7, 0, (b_u, b_u, b_u), 3, 0.0)
        filling = fl.TwoTargetStackFilling(8, 5, 0.0, filling)
        filling_0 = (
            fl.TwoTargetExteriorFilling(
                5, 9, (3, 3), 2, 0.0, filling, filling_0))
        filling = fl.TwoTargetFinalFilling(9, 5, 0, 1, 0.0, filling_0)
        self.assertEqual("AUANNNUGUUUGUGCANNNUNGCUNCANNNUANNNUUNAC",
                         str(filling))
        seq_list = ['N'] * len(self.target_structs[0])
        filling.get_seq_list(seq_list, self.target_structs[0])
        self.assertEqual("AUANNNUGUUUGUGCANNNUNGCUNCANNNUANNNUUNAC",
                         ''.join(seq_list))
        subsol = self.subsolutions[0][-1]
        base_ids = (subsol.get_base(filling, pos) for pos in
                    (0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 19, 21, 22,
                     23, 25, 26, 30, 31, 35, 36, 38, 39))
        self.assertListEqual(['A', 'U', 'A', 'U', 'G', 'U', 'U', 'U', 'G', 'U',
                              'G', 'C', 'A', 'U', 'G', 'C', 'U', 'C', 'A', 'U',
                              'A', 'U', 'U', 'A', 'C'],
                             [rna.BASES[b_id] for b_id in base_ids])

if __name__ == '__main__':
    unittest.main()
