import multiprocessing as mp
import os.path

from . import settings
from ..inverse_folding import two_target_inverse_fold as inverse_fold
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
                    self.q_out.put(
                        (parent_id, int(riboswitch_id), riboswitch_str))


class Generator(mp.Process):
    def __init__(self, in_, q_out):
        super(Generator, self).__init__()
        self.in_ = in_
        self.q_out = q_out

    def run(self):
        while True:
            task = self.in_.get()
            if task is None:
                break
            parent_id, riboswitch_id, riboswitch_str = task
            riboswitch = rs.get_riboswitch_from_str(riboswitch_str)
            riboswitch_fullstr = riboswitch.get_constraints_riboswitch()
            seqs = [seq for seq, _, _ in
                    inverse_fold.generate_sequences(*riboswitch_fullstr,
                                                    local_refinement=False)]
            self.q_out.put((parent_id, riboswitch_id, seqs))


class SequencesWriter(mp.Process):
    def __init__(self, q_in, output_dir):
        super(SequencesWriter, self).__init__()
        self.q_in = q_in
        self.output_dir = output_dir

    def run(self):
        while True:
            task = self.q_in.get()
            if task is None:
                break
            parent_id, riboswitch_id, seqs = task
            print "Writing sequences for %i (%i)" % (riboswitch_id, parent_id)
            group_id = parent_id % settings.NUM_PARENT_GROUPS
            parent_folder = os.path.join(self.output_dir,
                                         str(group_id), str(parent_id))
            if len(seqs):
                riboswitches_eval = os.path.join(parent_folder,
                                                 "seqs_%i" % riboswitch_id)
                with open(riboswitches_eval, 'w') as fh:
                    for seq in seqs:
                        fh.write("%s\n" % seq)
