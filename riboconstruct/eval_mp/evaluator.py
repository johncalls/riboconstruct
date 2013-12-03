import multiprocessing as mp
import operator
import os.path
import shutil

from . import settings
from .. import helper
from .. import riboswitch as rs
from .. import rna_f
from ..riboswitch import element as rs_e


rna_f = rna_f.RNAf()


def load_eval_dir(output_dir, q_out, get_best_num_seqs=None):
    def parent_ids():
        for group_id in os.listdir(output_dir):
            for parent_id in os.listdir(os.path.join(output_dir, group_id)):
                yield int(parent_id)

    for parent_id in sorted(parent_ids()):
        group_id = parent_id % settings.NUM_PARENT_GROUPS
        parent_folder = os.path.join(output_dir, str(group_id), str(parent_id))
        siblings_file = os.path.join(parent_folder, "siblings")
        with open(siblings_file) as fh:
            for line in fh:
                riboswitch_id, riboswitch_str = line.rstrip().split(' ', 1)
                riboswitch_id = int(riboswitch_id)
                riboswitch_seqs = os.path.join(parent_folder,
                                               "seqs_%i" % riboswitch_id)
                riboswitch_eval = os.path.join(parent_folder,
                                               "rnaf_eval_%i" % riboswitch_id)
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
                    # sort by cost
                    seqs = [(i, seq)
                            for i, seq, cost
                            in sorted(seqs, key=lambda x: x[2])]
                    if get_best_num_seqs is None:
                        q_out.put((parent_id, riboswitch_id, riboswitch_str,
                                   seqs))
                    else:
                        q_out.put((parent_id, riboswitch_id, riboswitch_str,
                                   seqs[:get_best_num_seqs]))


def get_features(evaluation_output):
    data = evaluation_output[rs_e.State.bound].split(';')
    b_accuracy_a = float(data[2]), float(data[4])
    b_accuracy_h_c = float(data[7]), float(data[9])
    b_accuracy_h_uc = float(data[14]), float(data[16])

    b_access_rs_c = 1.0 - float(data[11])  # fix accessibility
    b_access_rs_uc = 1.0 - float(data[18])  # fix accessibility

    data = evaluation_output[rs_e.State.unbound].split(';')
    ub_accuracy_h = float(data[7]), float(data[11])
    ub_accuracy_a = float(data[9]), float(data[13])
    ub_access_rs = 1.0 - float(data[15])  # fix accessibility
    ub_access_a = 1.0 - float(data[17])  # fix accessibility

    return (b_accuracy_a[0],
            b_accuracy_h_c[0], b_accuracy_h_uc[0],
            b_access_rs_c, b_access_rs_uc,

            ub_accuracy_h[0],
            ub_accuracy_a[0],
            ub_access_rs,
            ub_access_a)


def calculate_score(features):
    b_accuracy_h_c = features[1]  # maximize
    b_access_rs_c = 1.0 - features[3]  # minimize

    ub_accuracy_h = features[5]  # maximize
    ub_accuracy_a = features[6]  # maximize
    ub_access_rs = features[7]  # maximize
    ub_access_a = features[8]  # maximize

    objectives = (b_accuracy_h_c, b_access_rs_c, ub_accuracy_h, ub_accuracy_a,
                  ub_access_rs, ub_access_a)

    # score: multiply feature values
    return reduce(operator.mul, objectives)


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
            riboswitch_eval = os.path.join(self.output_dir,
                                           str(group_id), str(parent_id),
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
                                   helper.PATTERN_FILE_NAME),
                    rna_f.evaluate(evaluation_folders[1],
                                   helper.PATTERN_FILE_NAME))
                # extract features
                features = get_features(evaluation_output)
                # calculate score
                scores.append((calculate_score(features), features))
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
            riboswitch_eval = os.path.join(self.output_dir,
                                           str(group_id), str(parent_id),
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


class EvaluatorWriter(mp.Process):
    def __init__(self, q_in, output_dir):
        super(EvaluatorWriter, self).__init__()
        self.q_in = q_in
        self.output_dir = output_dir

    def run(self):
        while True:
            task = self.q_in.get()
            if task is None:
                break
            parent_id, riboswitch_id, seqs, scores = task
            group_id = parent_id % settings.NUM_PARENT_GROUPS
            riboswitch_eval = os.path.join(self.output_dir,
                                           str(group_id), str(parent_id),
                                           "eval_%i" % riboswitch_id)
            with open(riboswitch_eval, 'w') as fh:
                for seq, (score, features) in zip(seqs, scores):
                    fh.write("%f %s %s\n" % (score, features, seq))
