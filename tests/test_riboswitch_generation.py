#!/usr/bin/env python

from riboconstruct import riboswitch as rs
from riboconstruct.riboswitch import generation as rs_gen


initial_rs = (
    rs.get_riboswitch_from_config_file(
        "/home/john/projects/riboswitch_construction/riboconstruct_git/tests/"
        "initial.rs"))
ss_instance_space = rs_gen.SpliceSiteInstanceSpace()
ss_instance_space_iter = rs_gen.InstanceIterator(ss_instance_space, initial_rs)

for riboswitch in ss_instance_space_iter:
    print riboswitch
