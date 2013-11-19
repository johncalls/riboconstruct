import collections
import itertools
import operator
import os
import os.path
import shutil

from . import helper
from . import riboswitch as rs

from .inverse_folding import two_target_inverse_fold as inverse_fold

from .riboswitch import element as rs_e
from .riboswitch import siblings_generator as sg

from . import rna_f


OUTPUT_DIR = None


rna_f = rna_f.RNAf()


def set_output_dir(output_dir):
    if not os.path.isdir(output_dir):
        raise ValueError("Ouput folder '%s' does not exist." % output_dir)
    global OUTPUT_DIR
    OUTPUT_DIR = output_dir


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


def evaluate_riboswitch(riboswitch_id, riboswitch, parent_folder):
    out = []
    # generate riboswitch folder
    riboswitch_folder = os.path.join(parent_folder, str(riboswitch_id))
    os.mkdir(riboswitch_folder)
    # generate the hairpin's sequence
    constraints_fs = riboswitch.get_constraints_riboswitch()
    seqs_info = inverse_fold.generate_sequences(*constraints_fs)
    # evaluate the current riboswitch for the calculated sequences
    for seq_id, (sequence, cost, steps) in enumerate(seqs_info):
        sequence = str(sequence)
        # set the expression platform sequence
        ep_seq = rs_e.SequenceConstraint(riboswitch.pos_instance,
                                         sequence)
        riboswitch.add(ep_seq)
        # evaluate
        instance_folder = os.path.join(riboswitch_folder, str(seq_id))
        os.mkdir(instance_folder)
        evaluation_folders = riboswitch.create_evaluation_files(instance_folder)
        evaluation_output = (
            rna_f.evaluate(evaluation_folders[0], helper.PATTERN_FILE_NAME),
            rna_f.evaluate(evaluation_folders[1], helper.PATTERN_FILE_NAME))
        # extract features
        features = get_features(evaluation_output)
        # calculate score
        score = calculate_score(features)
        # remove expression platform sequence for next iteration
        riboswitch.remove(ep_seq)
        # add results
        out.append((score, features, sequence))
    # delete evaluation folders
    shutil.rmtree(riboswitch_folder)
    return out


def evaluate():
    def _eval(known_riboswitches, unevaluated_siblings, folder, start_id=0):
        riboswitch_id = start_id
        while len(unevaluated_siblings):
            parent_id, parent_str = unevaluated_siblings.popleft()
            parent = rs.get_riboswitch_from_str(parent_str)
            parent_folder = os.path.join(folder, str(parent_id))
            os.mkdir(parent_folder)
            with open(os.path.join(parent_folder, "parent"), 'w') as fh:
                fh.write("%s\n" % parent_str)
            with open(os.path.join(parent_folder, "riboswitches"), 'w') as fh:
                # improve the current hairpin model
                siblings_generator = sg.SimpleSiblingsGenerator(parent)
                for riboswitch in siblings_generator.iter_siblings():
                    riboswitch_str = str(riboswitch)
                    if riboswitch_str in known_riboswitches:
                        continue
                    riboswitch_id += 1
                    print ("Evaluating riboswitch %i (parent %i)" % (
                           riboswitch_id, parent_id))
                    fh.write("%i %s\n" % (riboswitch_id, riboswitch_str))
                    fh.flush()
                    output = evaluate_riboswitch(riboswitch_id, riboswitch,
                                                 parent_folder)
                    for score, features, ep_seq in output:
                        fh.write("\t%f %s %s\n" % (score, str(features),
                                                   str(ep_seq)))
                    known_riboswitches.add(riboswitch_str)
                    unevaluated_siblings.append(
                        (riboswitch_id, riboswitch_str))

    ####################################################################
    config_file = os.path.join(OUTPUT_DIR, "config.dat")
    if not os.path.isfile(config_file):
        raise ValueError("Output folder does not contain the configuration "
                         "file.")

    riboswitches_output_dir = os.path.join(OUTPUT_DIR, "riboswitches")
    if not os.path.exists(riboswitches_output_dir):
        os.makedirs(riboswitches_output_dir)

    known_riboswitches = set()
    unevaluated_siblings = collections.deque()

    ####################################################################
    riboswitch_id = 0
    riboswitch = helper.get_riboswitch_from_config_file(config_file)
    riboswitch_str = str(riboswitch)

    parent_id = str(riboswitch_id - 1)
    parent_folder = os.path.join(riboswitches_output_dir, parent_id)
    os.mkdir(parent_folder)

    unevaluated_siblings.append((riboswitch_id, riboswitch_str))
    known_riboswitches.add(riboswitch_str)

    with open(os.path.join(parent_folder, "riboswitches"), 'w') as fh:
        fh.write("%i %s\n" % (riboswitch_id, riboswitch_str))
        fh.flush()
        output = evaluate_riboswitch(riboswitch_id, riboswitch, parent_folder)
        for score, features, ep_seq in output:
            fh.write("\t%f %s %s\n" % (score, str(features), str(ep_seq)))

    _eval(known_riboswitches, unevaluated_siblings, riboswitches_output_dir)
