import itertools
import operator
import os
import os.path

from . import settings


def threshold_all(data, limits, comparison=operator.gt):
    for entry in data:
        features = transform_features(entry[1])
        if all(comparison(f, l) for f, l in itertools.izip(features, limits)):
            yield entry


def threshold_any(data, limits, comparison=operator.gt):
    for entry in data:
        features = transform_features(entry[1])
        if any(comparison(f, l) for f, l in itertools.izip(features, limits)):
            yield entry


def best_x(data, scoring, x=20):
    if x is None:
        x = len(data)
    return sorted(data, key=scoring, reverse=True)[:x]


def transform_features(features):
    b_accuracy_h_c = features[1]  # maximize
    b_access_rs_c = 1.0 - features[3]  # minimize

    ub_accuracy_h = features[5]  # maximize
    ub_accuracy_a = features[6]  # maximize
    ub_access_rs = features[7]  # maximize
    ub_access_a = features[8]  # maximize

    return (b_accuracy_h_c, b_access_rs_c,
            ub_accuracy_h, ub_accuracy_a, ub_access_rs, ub_access_a)


def linear_combination(entry):
    return reduce(operator.mul, transform_features(entry[1]))


def iter_eval_data(eval_folder):
    def parent_ids():
        for group_id in os.listdir(eval_folder):
            group_folder = os.path.join(eval_folder, group_id)
            for parent_id in os.listdir(group_folder):
                yield int(parent_id)

    for parent_id in sorted(parent_ids()):
        group_id = parent_id % settings.NUM_PARENT_GROUPS
        parent_folder = os.path.join(eval_folder, str(group_id),
                                     str(parent_id))
        siblings_file = os.path.join(parent_folder, "siblings")
        with open(siblings_file) as fh:
            for line in fh:
                riboswitch_id, _ = line.rstrip().split(' ', 1)
                riboswitch_id = int(riboswitch_id)
                try:
                    riboswitch_eval = os.path.join(parent_folder,
                                                   "eval_%i" % riboswitch_id)
                    with open(riboswitch_eval) as fh:
                        line = fh.next().rstrip().split(' ')
                        features = [float(line[1][1:-1])]
                        features.extend(float(f[:-1]) for f in line[2:10])
                        seq_id = int(line[10][1:-1])
                        seq = line[11][1:-2]
                    riboswitch_eval = (
                        os.path.join(parent_folder,
                                     "eval_plfold_%i" % riboswitch_id))
                    with open(riboswitch_eval) as fh:
                        line = fh.next().rstrip().split(' ')
                        # seq_id = int(line[0])
                        features_plfold = float(line[1]), float(line[2])
                except IOError:
                    continue
                yield ((parent_id, riboswitch_id), features, features_plfold,
                       (seq_id, seq))


def get_riboswitch_str(eval_folder, parent_id, riboswitch_id):
    group_id = parent_id % settings.NUM_PARENT_GROUPS
    siblings_file = os.path.join(eval_folder, str(group_id), str(parent_id),
                                 "siblings")
    with open(siblings_file) as fh:
        for line in fh:
            r_id, r_str = line.rstrip().split(' ', 1)
            if riboswitch_id == int(r_id):
                return r_str


def get_eval_seqs(eval_folder, parent_id, riboswitch_id):
    group_id = parent_id % settings.NUM_PARENT_GROUPS
    eval_file = os.path.join(eval_folder, str(group_id), str(parent_id),
                             "eval_%i" % riboswitch_id)
    with open(eval_file) as fh:
        for line in fh:
            _, seq_id, seq = line.rstrip().rsplit(' ', 2)
            yield int(seq_id[1:-1]), seq[1:-2]
