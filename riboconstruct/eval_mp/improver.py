import multiprocessing as mp
import os.path

from . import settings
from .. import rna
from ..inverse_folding import structure
from ..inverse_folding import two_target_inverse_fold as inverse_fold
from ..inverse_folding import two_target_local_refinement as loc_ref
from ..riboswitch import riboswitch as rs


class Loader(mp.Process):
    def __init__(self, q_out, output_dir):
        super(Loader, self).__init__()
        self.q_out = q_out
        self.output_dir = output_dir

    def run(self):
        def parent_ids():
            for group_id in os.listdir(self.output_dir):
                group_folder = os.path.join(self.output_dir, group_id)
                for parent_id in os.listdir(group_folder):
                    yield int(parent_id)

        for parent_id in sorted(parent_ids()):
            group_id = parent_id % settings.NUM_PARENT_GROUPS
            parent_folder = os.path.join(self.output_dir,
                                         str(group_id), str(parent_id))
            siblings_file = os.path.join(parent_folder, "siblings")
            with open(siblings_file) as fh:
                for line in fh:
                    riboswitch_id, riboswitch_str = line.rstrip().split(' ', 1)
                    riboswitch_id = int(riboswitch_id)
                    riboswitch_seqs = os.path.join(parent_folder,
                                                   "seqs_%i" % riboswitch_id)
                    seqs = []
                    try:
                        with open(riboswitch_seqs) as fh2:
                            for line in fh2:
                                seq = line.rstrip().split(' ', 1)[0]
                                seqs.append(seq)
                    except IOError:
                        continue
                    if len(seqs):
                        self.q_out.put(
                            (parent_id, riboswitch_id, riboswitch_str, seqs))


class Improver(mp.Process):
    def __init__(self, q_in, q_out):
        super(Improver, self).__init__()
        self.q_in = q_in
        self.q_out = q_out

    def run(self):
        while True:
            task = self.q_in.get()
            if task is None:
                break
            parent_id, riboswitch_id, riboswitch_str, seqs = task
            riboswitch = rs.get_riboswitch_from_str(riboswitch_str)
            targets, seq_constraint = riboswitch.get_constraints_riboswitch()
            targets = (structure.Structure(targets[0]),
                       structure.Structure(targets[1]))
            seq_constraint = (
                inverse_fold.calculate_sequence_constraint(
                    targets, rna.IUPACSequence(seq_constraint)))
            loc_ref.MAX_STEPS = 10 * len(targets[0])
            seqs = [loc_ref.local_search(seq, targets, seq_constraint)
                    for seq in seqs]
            self.q_out.put((parent_id, riboswitch_id, seqs))


class ImprovementsWriter(mp.Process):
    def __init__(self, q_in, output_dir):
        super(ImprovementsWriter, self).__init__()
        self.q_in = q_in
        self.output_dir = output_dir

    def run(self):
        while True:
            task = self.q_in.get()
            if task is None:
                break
            parent_id, riboswitch_id, seqs = task
            print "Writing improved sequences of %i (%i)" % (riboswitch_id,
                                                             parent_id)
            group_id = parent_id % settings.NUM_PARENT_GROUPS
            riboswitch_seqs = os.path.join(self.output_dir,
                                           str(group_id), str(parent_id),
                                           "seqs_%i" % riboswitch_id)
            with open(riboswitch_seqs, 'w') as fh:
                for seq, cost, steps in seqs:
                    fh.write("%s %f\n" % (seq, cost))
