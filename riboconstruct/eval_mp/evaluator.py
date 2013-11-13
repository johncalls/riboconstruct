import multiprocessing as mp
import os.path
import shutil

from . import settings
from .. import main
from .. import helper
from .. import params
from .. import riboswitch as rs
from .. import rna_f
from ..riboswitch import element as rs_e


rna_f = rna_f.RNAf()


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
                    riboswitch_eval = (
                        os.path.join(parent_folder,
                                     "rnaf_eval_%i" % riboswitch_id))
                    if os.path.exists(riboswitch_eval):
                        shutil.rmtree(riboswitch_eval)
                    seqs = []
                    try:
                        with open(riboswitch_seqs) as fh2:
                            for i, line in enumerate(fh2):
                                seq, cost = line.rstrip().split(' ', 1)
                                seqs.append((i, seq, cost))
                    except IOError:
                        continue
                    if len(seqs):
                        seqs = [(i, seq)
                                for i, seq, cost
                                in sorted(seqs, key=lambda x: x[2])]
                        self.q_out.put(
                            (parent_id, riboswitch_id, riboswitch_str,
                             (seqs[0],)))


class Evaluator(mp.Process):
    def __init__(self, q_in, q_out, output_dir):
        super(Evaluator, self).__init__()
        self.q_in = q_in
        self.q_out = q_out
        self.output_dir = output_dir

    def run(self):
        while True:
            task = self.q_in.get()
            if task is None:
                break
            parent_id, riboswitch_id, riboswitch_str, seqs = task
            group_id = parent_id % settings.NUM_PARENT_GROUPS
            parent_folder = os.path.join(self.output_dir,
                                         str(group_id), str(parent_id))
            riboswitch_eval = os.path.join(parent_folder,
                                           "rnaf_eval_%i" % riboswitch_id)
            riboswitch = rs.get_riboswitch_from_str(riboswitch_str)
            os.mkdir(riboswitch_eval)
            scores = []
            for seq_id, seq in seqs:
                instance_folder = os.path.join(riboswitch_eval, str(seq_id))
                os.mkdir(instance_folder)
                seq = rs_e.SequenceConstraint(riboswitch.pos_riboswitch, seq)
                riboswitch.add(seq)
                evaluation_folders = (
                    riboswitch.create_evaluation_files(instance_folder))
                # evaluate
                evaluation_output = (
                    rna_f.evaluate(evaluation_folders[0],
                                   params.PATTERN_FILE_NAME),
                    rna_f.evaluate(evaluation_folders[1],
                                   params.PATTERN_FILE_NAME))
                # extract features
                features = main.get_features(evaluation_output)
                # calculate score
                scores.append((main.calculate_score(features), features))
                riboswitch.remove(seq)
            self.q_out.put((parent_id, riboswitch_id, seqs, scores))
            shutil.rmtree(riboswitch_eval)


class EvaluatorRNAplfold(mp.Process):
    def __init__(self, q_in, output_dir):
        super(EvaluatorRNAplfold, self).__init__()
        self.q_in = q_in
        self.output_dir = output_dir

    def run(self):
        while True:
            task = self.q_in.get()
            if task is None:
                break
            parent_id, riboswitch_id, riboswitch_str, seqs = task
            group_id = parent_id % settings.NUM_PARENT_GROUPS
            parent_folder = os.path.join(self.output_dir,
                                         str(group_id), str(parent_id))
            riboswitch_eval = os.path.join(parent_folder,
                                           "eval_plfold_%i" % riboswitch_id)
            riboswitch = rs.get_riboswitch_from_str(riboswitch_str)
            with open(riboswitch_eval, 'w') as fh:
                for seq_id, seq in seqs:
                    riboswitch.add(
                        rs_e.SequenceConstraint(riboswitch.pos_riboswitch, seq))
                    _, seq = riboswitch.get_constraints()
                    (positions,
                     lengths,
                     l) = riboswitch.get_accessibility_positions()
                    accessibilities = helper.rnaplfold(seq, positions, lengths,
                                                       L=l)
                    fh.write("%i %f %f\n" % (seq_id,
                                             accessibilities[0],
                                             accessibilities[1]))


class EvaluationWriter(mp.Process):
    def __init__(self, q_in, output_dir):
        super(EvaluationWriter, self).__init__()
        self.q_in = q_in
        self.output_dir = output_dir

    def run(self):
        while True:
            task = self.q_in.get()
            if task is None:
                break
            parent_id, riboswitch_id, seqs, scores = task
            print "Writing evaluated riboswitch %i (%i)" % (riboswitch_id,
                                                            parent_id)
            group_id = parent_id % settings.NUM_PARENT_GROUPS
            riboswitch_eval = os.path.join(self.output_dir,
                                           str(group_id), str(parent_id),
                                           "eval_%i" % riboswitch_id)
            with open(riboswitch_eval, 'w') as fh:
                for seq, (score, features) in zip(seqs, scores):
                    fh.write("%f %s %s\n" % (score, features, seq))
