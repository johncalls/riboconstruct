import collections
import itertools
import operator
import os
import os.path
import shutil

from . import helper
from . import params
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
    constraints_fs = riboswitch.get_constraints_functional_site()
    seqs_info = inverse_fold.generate_sequences(*constraints_fs)
    # evaluate the current riboswitch for the calculated sequences
    for seq_id, (sequence, cost, steps) in enumerate(seqs_info):
        sequence = str(sequence)
        # set the expression platform sequence
        ep_seq = rs_e.SequenceConstraint(riboswitch.pos_functional_site,
                                         sequence)
        riboswitch.add(ep_seq)
        # evaluate
        instance_folder = os.path.join(riboswitch_folder, str(seq_id))
        os.mkdir(instance_folder)
        evaluation_folders = riboswitch.create_evaluation_files(instance_folder)
        evaluation_output = (
            rna_f.evaluate(evaluation_folders[0], params.PATTERN_FILE_NAME),
            rna_f.evaluate(evaluation_folders[1], params.PATTERN_FILE_NAME))
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


def evaluate(known_riboswitches, unevaluated_siblings, folder, start_id=0):
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
                print "Evaluating riboswitch %i (parent %i)" % (riboswitch_id,
                                                                parent_id)
                fh.write("%i %s\n" % (riboswitch_id, riboswitch_str))
                fh.flush()
                output = evaluate_riboswitch(riboswitch_id, riboswitch,
                                             parent_folder)
                for score, features, ep_seq in output:
                    fh.write("\t%f %s %s\n" % (score,
                                               str(features), str(ep_seq)))
                known_riboswitches.add(riboswitch_str)
                unevaluated_siblings.append((riboswitch_id, riboswitch_str))


def main():
    ############################################################################
    config_file = os.path.join(OUTPUT_DIR, "config.dat")
    if not os.path.isfile(config_file):
        raise ValueError("Output folder does not contain the configuration "
                         "file.")

    riboswitches_output_dir = os.path.join(OUTPUT_DIR, "riboswitches")
    if not os.path.exists(riboswitches_output_dir):
        os.makedirs(riboswitches_output_dir)

    known_riboswitches = set()
    unevaluated_siblings = collections.deque()

    ############################################################################
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

    evaluate(known_riboswitches, unevaluated_siblings, riboswitches_output_dir)


def iter_parent_eval(eval_data):
    evaluation_output = []
    riboswitch_id = riboswitch_str = None
    for line in eval_data:
        if line[0] == '\t':
            evaluation_output.append(line)
        else:
            if riboswitch_id is not None:
                if len(evaluation_output):
                    yield (True, riboswitch_id, riboswitch_str,
                           evaluation_output)
                    evaluation_output = []
                else:
                    yield False, riboswitch_id, riboswitch_str, None
            riboswitch_id, riboswitch_str = line.rstrip().split(None, 1)
            riboswitch_id = int(riboswitch_id)
    if riboswitch_id is not None:
        if len(evaluation_output):
            yield True, riboswitch_id, riboswitch_str, evaluation_output
        else:
            yield False, riboswitch_id, riboswitch_str, None


def fix_eval_data(folder):
    reexpand = []
    known_riboswitches = set()
    highest_r_id = 0

    print "Load re-evaluate data; re-evaluate missing entries"
    for parent_id in (int(i) for i in os.listdir(folder)):  # unsorted iteration

        parent_folder = os.path.join(folder, str(parent_id))
        riboswitches_file = os.path.join(parent_folder, "riboswitches")
        # check if parent has been expanded, i.e. a riboswitches file exists
        if ((not os.path.exists(riboswitches_file) or
             os.stat(riboswitches_file).st_size == 0) and
            parent_id != -1):
            # save for re-evaluation if the parent has not been expanded
            parent_file = os.path.join(parent_folder, "parent")
            with open(parent_file) as fh:
                parent_str = fh.readline().rstrip()
            reexpand.append((parent_id, parent_str))
        # check siblings
        cleaned_data = []

        with open(riboswitches_file) as fh:
            for complete, r_id, r_str, eval_data in iter_parent_eval(fh):
                if r_str in known_riboswitches:
                    raise ValueError("Riboswitch %i (parent: %i) seen more "
                                     "than once." % (r_id, parent_id))
                if r_id > highest_r_id:
                    highest_r_id = r_id
                known_riboswitches.add(r_str)

                cleaned_data.append("%i %s\n" % (r_id, r_str))
                if complete:
                    cleaned_data.extend(eval_data)
                else:  # re-evaluate
                    # remove the sibling folder (which should not exist anyway)
                    riboswitch_folder = os.path.join(parent_folder, str(r_id))
                    if os.path.exists(riboswitch_folder):
                        shutil.rmtree(riboswitch_folder)
                    # re-evaluate
                    riboswitch = rs.get_riboswitch_from_str(r_str)
                    output = evaluate_riboswitch(r_id, riboswitch,
                                                 parent_folder)
                    for score, features, ep_seq in output:
                        cleaned_data.append("\t%f %s %s\n" % (
                                            score, str(features), str(ep_seq)))
        # overwrite with cleaned data
        if len(cleaned_data):
            with open(riboswitches_file, 'w') as fh:
                for line in cleaned_data:
                    fh.write(line)
    # re-expand parent
    print "Re-expand riboswitches"
    for parent_id, parent_str in reexpand:
        parent = rs.get_riboswitch_from_str(parent_str)
        parent_folder = os.path.join(folder, str(parent_id))
        with open(os.path.join(parent_folder, "parent"), 'w') as fh:
            fh.write("%s\n" % parent_str)
        with open(os.path.join(parent_folder, "riboswitches"), 'w') as fh:
            siblings_generator = sg.SimpleSiblingsGenerator(parent)
            for riboswitch in siblings_generator.iter_siblings():
                riboswitch_str = str(riboswitch)
                if riboswitch_str in known_riboswitches:
                    continue
                highest_r_id += 1

                fh.write("%i %s\n" % (highest_r_id, riboswitch_str))
                fh.flush()
                output = evaluate_riboswitch(highest_r_id, riboswitch,
                                             parent_folder)
                for score, features, ep_seq in output:
                    fh.write("\t%f %s %s\n" % (score, str(features),
                                               str(ep_seq)))


def load_eval_data(folder, initial_riboswitch):
    known_riboswitches = set()
    unevaluated_siblings = dict()
    not_found = collections.deque()
    highest_r_id = 0

    parent_ids = sorted(int(i) for i in os.listdir(folder))
    if parent_ids[0] != -1:
        raise ValueError("Initial riboswitch missing")
    initial_riboswitch_str = str(initial_riboswitch)
    known_riboswitches.add(initial_riboswitch_str)
    unevaluated_siblings[initial_riboswitch_str] = 0

    print "Load data"
    for parent_id in itertools.islice(parent_ids, 1, None):
        parent_folder = os.path.join(folder, str(parent_id))
        with open(os.path.join(parent_folder, "parent")) as fh:
            parent_str = fh.readline().rstrip()
        del unevaluated_siblings[parent_str]
        parent = rs.get_riboswitch_from_str(parent_str)
        # list of expected siblings
        siblings = set(str(s)
                       for s
                       in sg.SimpleSiblingsGenerator(parent).iter_siblings())
        riboswitches_file = os.path.join(parent_folder, "riboswitches")
        if not os.path.exists(riboswitches_file):
            print parent_id
        with open(riboswitches_file) as fh:
            for _, r_id, r_str, _ in iter_parent_eval(fh):
                if r_id > highest_r_id:
                    highest_r_id = r_id
                known_riboswitches.add(r_str)
                unevaluated_siblings[r_str] = r_id
                # remove the riboswitch since it has been found
                siblings.remove(r_str)
        # mark all siblings that have not been found (but might be found later)
        for riboswitch_str in siblings:
            not_found.append((parent_id, riboswitch_str))
    # check the siblings that have not been found
    for parent_id, riboswitch_str in not_found:
        if riboswitch_str in known_riboswitches:
            continue
        known_riboswitches.add(riboswitch_str)
        highest_r_id += 1
        parent_folder = os.path.join(folder, str(parent_id))
        riboswitch = rs.get_riboswitch_from_str(riboswitch_str)
        with open(os.path.join(parent_folder, "riboswitches"), 'a') as fh:
            fh.write("%i %s\n" % (highest_r_id, riboswitch_str))
            fh.flush()
            output = evaluate_riboswitch(highest_r_id, riboswitch,
                                         parent_folder)
            for score, features, ep_seq in output:
                fh.write("\t%f %s %s\n" % (score, str(features),
                                           str(ep_seq)))
    # save all riboswitches that have not been expanded
    unevaluated_siblings = (
        collections.deque(
            sorted((r_id, r_str)
                   for r_str, r_id in unevaluated_siblings.iteritems())))

    return highest_r_id, known_riboswitches, unevaluated_siblings


def main_continued():
    ###########################################################################
    config_file = os.path.join(OUTPUT_DIR, "config.dat")
    if not os.path.isfile(config_file):
        raise ValueError("Output folder does not contain the configuration "
                         "file.")

    riboswitches_output_dir = os.path.join(OUTPUT_DIR, "riboswitches")
    if not os.path.exists(riboswitches_output_dir):
        os.makedirs(riboswitches_output_dir)

    initial_riboswitch = helper.get_riboswitch_from_config_file(config_file)

    ###########################################################################
    fix_eval_data(riboswitches_output_dir)
    riboswitch_id, known_riboswitches, unevaluated_siblings = (
        load_eval_data(riboswitches_output_dir, initial_riboswitch))

    print "Evaluate"
    evaluate(known_riboswitches, unevaluated_siblings, riboswitches_output_dir,
             riboswitch_id)
