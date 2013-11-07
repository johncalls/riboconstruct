from . import subsolution
from .fillings import bulge
from .fillings import exterior
from .fillings import hairpin
from .fillings import interior
from .fillings import multiloop
from .fillings import stacking

from .. import rna


def calculate_sequence_constraint(struct, iupac_seq):
    # check if IUPAC sequence is a valid constraint for the given structure
    for bp_pos_i, bp_pos_j in struct.bp_positions:
        possible = False
        for bp_id_i, bp_id_j in rna.BASEPAIRS:
            iupac_id_i = getattr(rna.IUPAC_Id, iupac_seq[bp_pos_i])
            iupac_id_j = getattr(rna.IUPAC_Id, iupac_seq[bp_pos_j])
            if (rna.base_valid_IUPAC(iupac_id_i, bp_id_i) and
                rna.base_valid_IUPAC(iupac_id_j, bp_id_j)):
                possible = True
                break
        if not possible:
            raise ValueError(
                'Constraint (%s, %s) at base pair (%i, %i) not '
                'compatible.' %
                (iupac_seq[bp_pos_i], iupac_seq[bp_pos_j],
                 bp_pos_i, bp_pos_j))

    seq_constraint = {}
    for pos, iu in enumerate(iupac_seq):
        iupac_id = getattr(rna.IUPAC_Id, iu)
        bp_pos_i = struct.basepairs[pos]
        if bp_pos_i is None:
            pos_constraint = []
            # fill sequence constraint
            for b_id in xrange(rna.BaseId.count):
                invalid = not rna.base_valid_IUPAC(iupac_id, b_id)
                pos_constraint.append(invalid)
        else:
            if bp_pos_i < pos:
                pos_constraint_i = [True] * rna.BaseId.count
                pos_constraint = [True] * rna.BaseId.count
                iupac_id_i = getattr(rna.IUPAC_Id, iupac_seq[bp_pos_i])
                for bp_id_i, bp_id_j in rna.BASEPAIRS:
                    invalid = not (
                        rna.base_valid_IUPAC(iupac_id_i, bp_id_i) and
                        rna.base_valid_IUPAC(iupac_id, bp_id_j))
                    pos_constraint_i[bp_id_i] &= invalid
                    pos_constraint[bp_id_j] &= invalid
                seq_constraint[bp_pos_i] = tuple(pos_constraint_i)
            else:
                continue
        seq_constraint[pos] = tuple(pos_constraint)
    return tuple((seq_constraint[i] for i in xrange(len(iupac_seq))))


def inverse_fold(target_struct, seq_constraint):
    def get_fillings_generator(substruct, seq_constraint, predecessor=None):
        if substruct.struct_type == rna.StructType.HAIRPIN:
            return hairpin.HairpinFillingsGenerator(
                subsolution.HairpinSubsolution(substruct, seq_constraint))
        if substruct.struct_type == rna.StructType.STACKING:
            return stacking.StackFillingsGenerator(
                subsolution.StackSubsolution(
                    substruct, seq_constraint, predecessor))
        elif substruct.struct_type == rna.StructType.BULGE:
            return bulge.BulgeFillingsGenerator(
                subsolution.BulgeSubsolution(
                    substruct, seq_constraint, predecessor))
        elif substruct.struct_type == rna.StructType.INTERIOR:
            return interior.InteriorFillingsGenerator(
                subsolution.InteriorSubsolution(
                    substruct, seq_constraint, predecessor))
        elif substruct.struct_type == rna.StructType.MULTILOOP:
            return multiloop.MultiloopFillingsGenerator(
                subsolution.MultiloopSubsolution(
                    substruct, seq_constraint, predecessor))
        elif substruct.struct_type == rna.StructType.EXTERIOR:
            return exterior.FinalFillingsGenerator(
                subsolution.FinalSubsolution(
                    substruct, seq_constraint, predecessor))

    ############################################################################

    substruct_iter = target_struct.iter_substructs()
    substruct = substruct_iter.next()
    fillings_generator = get_fillings_generator(substruct, seq_constraint)
    fillings_generator.evaluate()

    pred_subsol = fillings_generator.subsolution

    stem_ends = {}
    if substruct.is_stem_end:
        stem_ends[substruct.bp_pos_id] = pred_subsol

    for bp_pos_id, substruct in enumerate(substruct_iter, 1):
        if substruct.struct_type == rna.StructType.MULTILOOP:
            pred = []
            for stem_bp_pos_id in substruct.predecessor_pos_ids:
                pred.append(stem_ends[stem_bp_pos_id])
                del stem_ends[stem_bp_pos_id]
            fillings_generator = (
                get_fillings_generator(substruct, seq_constraint, pred))
            if substruct.is_stem_end:
                stem_ends[substruct.bp_pos_id] = fillings_generator.subsolution
        elif substruct.struct_type == rna.StructType.EXTERIOR:
            pred = []
            pred.extend(stem_ends.values())
            fillings_generator = (
                get_fillings_generator(substruct, seq_constraint, pred))
        elif substruct.struct_type == rna.StructType.HAIRPIN:
            fillings_generator = (
                get_fillings_generator(substruct, seq_constraint))
            if substruct.is_stem_end:
                stem_ends[substruct.bp_pos_id] = fillings_generator.subsolution
        else:
            fillings_generator = (
                get_fillings_generator(substruct, seq_constraint, pred_subsol))
            if substruct.is_stem_end:
                stem_ends[substruct.bp_pos_id] = fillings_generator.subsolution

        fillings_generator.evaluate()
        pred_subsol = fillings_generator.subsolution

    return pred_subsol
