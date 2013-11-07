#!/usr/bin/env python

import collections
import os
import os.path

from . import settings
from ..riboswitch import riboswitch as rs
from ..riboswitch import siblings_generator_2 as sg


def create_instances(initial_riboswitch, output_dir):
    parents = collections.deque()
    known_riboswitches = set()

    parent_id = -1
    group_id = parent_id % settings.NUM_PARENT_GROUPS
    parent_folder = os.path.join(output_dir, str(group_id), str(parent_id))
    os.makedirs(parent_folder)
    riboswitch_id = 0
    riboswitch_str = str(initial_riboswitch)
    with open(os.path.join(parent_folder, "siblings"), 'w') as fh:
        fh.write("%i %s\n" % (0, riboswitch_str))
    parents.append((riboswitch_id, riboswitch_str))
    known_riboswitches.add(riboswitch_str)

    while len(parents):
        parent_id, parent_str = parents.popleft()
        print "Expanding %i" % parent_id
        parent = rs.get_riboswitch_from_str(parent_str)
        # generate folders
        group_id = parent_id % settings.NUM_PARENT_GROUPS
        parent_folder = os.path.join(output_dir,
                                     str(group_id), str(parent_id))
        os.makedirs(parent_folder)
        # write parent infos
        with open(os.path.join(parent_folder, "parent"), 'w') as fh:
            fh.write("%s\n" % parent_str)
        # generate siblings
        siblings_generator = sg.SimpleSiblingsGenerator(parent)
        with open(os.path.join(parent_folder, "siblings"), 'w') as fh:
            for riboswitch in siblings_generator.iter_siblings():
                riboswitch_str = str(riboswitch)
                if riboswitch_str in known_riboswitches:
                    continue
                known_riboswitches.add(riboswitch_str)
                riboswitch_id += 1
                parents.append((riboswitch_id, riboswitch_str))
                fh.write("%i %s\n" % (riboswitch_id, riboswitch_str))


def update_instances(output_dir):
    def parent_ids():
        for group_id in os.listdir(output_dir):
            group_folder = os.path.join(output_dir, group_id)
            for parent_id in os.listdir(group_folder):
                yield int(parent_id)

    known_riboswitches = set()
    parents = collections.deque()
    max_r_id = 0

    # load instances
    for parent_id in sorted(parent_ids()):
        group_id = parent_id % settings.NUM_PARENT_GROUPS
        parent_folder = os.path.join(output_dir, str(group_id), str(parent_id))
        try:
            with open(os.path.join(parent_folder, "parent")) as fh:
                parent_str = fh.next().rstrip()
                known_riboswitches.add(parent_str)
                parents.append((parent_id, parent_str))
        except IOError:
            pass
        max_r_id = parent_id

    # reexpand
    while len(parents):
        parent_id, parent_str = parents.popleft()
        group_id = parent_id % settings.NUM_PARENT_GROUPS
        parent_folder = os.path.join(output_dir,
                                     str(group_id), str(parent_id))
        siblings_file = os.path.join(parent_folder, "siblings")
        # check if parent has been seen before; otherwise write parent infos
        if not os.path.exists(parent_folder):
            os.makedirs(parent_folder)
            with open(os.path.join(parent_folder, "parent"), 'w') as fh:
                fh.write("%s\n" % parent_str)
            with open(siblings_file, 'w'):
                pass
        # generate siblings
        parent = rs.get_riboswitch_from_str(parent_str)
        siblings_generator = sg.SimpleSiblingsGenerator(parent)
        with open(siblings_file, 'a') as fh:
            for riboswitch in siblings_generator.iter_siblings():
                riboswitch_str = str(riboswitch)
                if riboswitch_str in known_riboswitches:
                    continue
                known_riboswitches.add(riboswitch_str)
                max_r_id += 1
                parents.append((max_r_id, riboswitch_str))
                fh.write("%i %s\n" % (max_r_id, riboswitch_str))
                print "New sibling %i (%i)" % (max_r_id, parent_id)
