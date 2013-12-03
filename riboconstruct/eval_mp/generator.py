import multiprocessing as mp
import os.path

from . import settings
from .. import riboswitch as rs
from ..inverse_folding import two_target_inverse_fold as inverse_fold


def load_eval_dir(output_dir, q_out):
    def parent_ids():
        for group_id in os.listdir(output_dir):
            for parent_id in os.listdir(os.path.join(output_dir, group_id)):
                yield int(parent_id)

    for p_id in sorted(parent_ids()):
        group_id = p_id % settings.NUM_PARENT_GROUPS
        siblings_file = os.path.join(output_dir, str(group_id), str(p_id),
                                     "siblings")
        with open(siblings_file) as fh:
            for line in fh:
                r_id, r_str = line.rstrip().split(' ', 1)
                q_out.put((p_id, int(r_id), r_str))


class Generator(mp.Process):
    def __init__(self, q_in, q_out):
        super(Generator, self).__init__()
        self.q_in = q_in
        self.q_out = q_out

    def run(self):
        while True:
            task = self.q_in.get()
            if task is None:
                break
            parent_id, riboswitch_id, riboswitch_str = task
            riboswitch = rs.get_riboswitch_from_str(riboswitch_str)
            riboswitch_fullstr = riboswitch.get_constraints_riboswitch()
            seqs = [seq for seq, _, _ in
                    inverse_fold.generate_sequences(*riboswitch_fullstr,
                                                    local_refinement=False)]
            self.q_out.put((parent_id, riboswitch_id, seqs))


class GeneratorWriter(mp.Process):
    def __init__(self, q_in, output_dir):
        super(GeneratorWriter, self).__init__()
        self.q_in = q_in
        self.output_dir = output_dir

    def run(self):
        while True:
            task = self.q_in.get()
            if task is None:
                break
            parent_id, riboswitch_id, seqs = task
            group_id = parent_id % settings.NUM_PARENT_GROUPS
            if len(seqs):
                riboswitches_eval = os.path.join(self.output_dir,
                                                 str(group_id),
                                                 str(parent_id),
                                                 "seqs_%i" % riboswitch_id)
                with open(riboswitches_eval, 'w') as fh:
                    for seq in seqs:
                        fh.write("%s\n" % seq)
