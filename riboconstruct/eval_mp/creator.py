#!/usr/bin/env python

import os
import os.path

from . import settings


def build_eval_dir(instance_iter, output_dir, q_out=None):
    """
    Builds the evaluation directory within *output_dir*. This directory
    is organized in a parent-sibling way, that is, all sibling
    riboswitches derived from a certain parent riboswitch are stored
    within the same subdirectory identified by the parent's id.

    To be able to handle the vast amount of parent ids, these are
    organized previously in groups. The maximum number of groups is
    defined by :attr:`settings.NUM_PARENT_GROUPS`.

    If *q_out* is not ``None``, a tuple containing the parent id, the
    riboswitch id and the riboswitch string representation is put into
    the queue. *q_out* is expected to be a
    :class:`multiprocessing.Queue`.
    """
    siblings = []
    previous_p_id = -1
    for p_id, r_id, r in instance_iter:
        if p_id != previous_p_id:
            group_id = p_id % settings.NUM_PARENT_GROUPS
            parent_folder = os.path.join(output_dir, str(group_id),
                                         str(previous_p_id))
            os.makedirs(parent_folder)
            with open(os.path.join(parent_folder, "siblings"), 'w') as fh:
                for s_id, s in siblings:
                    fh.write("%i %s\n" % (s_id, s))
            if q_out is not None:
                q_out.put((p_id, s_id, str(s)))
            siblings = []
            previous_p_id = p_id
        siblings.append((r_id, r))
    if q_out is not None:
        q_out.put(None)
