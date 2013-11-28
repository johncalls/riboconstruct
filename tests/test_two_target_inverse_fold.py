#!/usr/bin/env python

import unittest

from riboconstruct import rna

from riboconstruct.inverse_folding import structure
from riboconstruct.inverse_folding import two_target_inverse_fold as inverse_fold

from riboconstruct.inverse_folding.two_target_subsolution import bulge
from riboconstruct.inverse_folding.two_target_subsolution import exterior
from riboconstruct.inverse_folding.two_target_subsolution import hairpin
from riboconstruct.inverse_folding.two_target_subsolution import interior
from riboconstruct.inverse_folding.two_target_subsolution import multiloop
from riboconstruct.inverse_folding.two_target_subsolution import stacking


def create_subsolution(struct_id, substruct, task, predecessor=None):
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


def get_subsolution(struct_id, substruct, task):
    if substruct.struct_type == rna.StructType.MULTILOOP:
        pred = []
        for stem_bp_pos_id in substruct.predecessor_pos_ids:
            pred.append(
                task.subsolution_collector[struct_id][stem_bp_pos_id])
            del task.subsolution_collector[struct_id][stem_bp_pos_id]
        subsol = create_subsolution(struct_id, substruct, task, pred)
        task.subsolution_collector[struct_id][substruct.bp_pos_id] = subsol
    elif substruct.struct_type == rna.StructType.EXTERIOR:
        pred = []
        for stem_bp_pos_id in substruct.predecessor_pos_ids:
            pred.append(
                task.subsolution_collector[struct_id][stem_bp_pos_id])
            del task.subsolution_collector[struct_id][stem_bp_pos_id]
        subsol = create_subsolution(struct_id, substruct, task, pred)
        task.subsolution_collector[struct_id][substruct.bp_pos_id] = subsol
    elif substruct.struct_type == rna.StructType.HAIRPIN:
        subsol = create_subsolution(struct_id, substruct, task)
        task.subsolution_collector[struct_id][substruct.bp_pos_id] = subsol
    else:
        pred = (
            task.subsolution_collector[struct_id][substruct.bp_pos_id - 1])
        del task.subsolution_collector[struct_id][substruct.bp_pos_id - 1]
        subsol = create_subsolution(struct_id, substruct, task, pred)
        task.subsolution_collector[struct_id][substruct.bp_pos_id] = subsol

    return subsol


class TestDependencies(unittest.TestCase):
    def test_dependencies(self):
        struct_0 = structure.Structure(".....(....).....")
        struct_1 = structure.Structure(".(...)..(......)")
        # struct_0 = structure.Structure("(...(...)..(...)..)")
        # struct_1 = structure.Structure("(.................)")
        seq_constraint = (
            inverse_fold.calculate_sequence_constraint(
                (struct_0, struct_1), 'N' * len(struct_0)))
        task = inverse_fold.Task((struct_0, struct_1), seq_constraint)
        for belonging, substruct in task.iter_substructs():
            if belonging == 2:
                raise ValueError("Reached invalid case.")
            else:
                subsol = get_subsolution(belonging, substruct, task)
                if subsol.substruct.struct_type == rna.StructType.HAIRPIN:
                    if subsol.struct_id == 0:
                        self.assertEqual(
                            [5, 6, 7, 8, 9], sorted(subsol.dep_b_positions))
                    else:  # struct_id == 1
                        if subsol.bp_pos_id == 0:
                            self.assertEqual(
                                [8, 9, 10, 11], sorted(subsol.dep_b_positions))
                        else:  # bp_pos_id == 1
                            self.assertEqual(
                                [4, 5], sorted(subsol.dep_b_positions))
                elif subsol.substruct.struct_type == rna.StructType.EXTERIOR:
                    if subsol.struct_id == 0:
                        self.assertEqual(
                            [0, 1, 5, 6, 7, 8, 9, 14, 15],
                            sorted(subsol.dep_b_positions))
                    else:  # struct_id == 1
                        self.assertEqual(
                            [4, 5, 6, 7, 8, 9, 10, 11],
                            sorted(subsol.dep_b_positions))
                else:
                    raise ValueError("Reached invalid case.")


if __name__ == "__main__":
    unittest.main()
