#!/usr/bin/env python

from riboconstruct import riboswitch as rs
from riboconstruct.riboswitch import generation as rs_gen


initial_riboswitch_str = '[0, 103];[0, 103];(A_ub (50, 103) "......(((((....)))))..(((((...........))))).........." "AACAUACCAGAGAAAUCUGGAGAGGUGAAGAAUACGACCACCUAGGNNNNNNN");(A_b (38, 103) "((((((((((........(((((....)))))..(((((...........)))))))))))))))" "NNNNNNNCCUAAAACAUACCAGAGAAAUCUGGAGAGGUGAAGAAUACGACCACCUAGGNNNNNNN");(H_ub (17, 34) "(((((((...)))))))");(H_b (0, 17) "(((((((...)))))))");(TS (0, 10) "AUGAGUAUGU");(AC (50, 57))'
initial_riboswitch = rs.get_riboswitch_from_str(initial_riboswitch_str)

# use the splice site instance space
i_s = rs_gen.SpliceSiteInstanceSpace()

# check whether the initial riboswitch is valid within the splice site
# instance space
if not i_s.validate(initial_riboswitch):
    print i_s.get_detailed_validation(initial_riboswitch)

# iterate the instances
instance_iter = rs_gen.InstanceIterator(i_s, initial_riboswitch)
for i, (p_id, r_id, riboswitch) in enumerate(instance_iter):
    # explicitly validate each instance again in each step
    # NOTE: this is _not_ necessary since the instance iterator only returns
    #       valid instances
    if not i_s.validate(riboswitch):
        print "ERROR: riboswitch not valid:",
        print p_id, r_id, riboswitch
        print i_s.get_detailed_validation(riboswitch)
        break
    else:
        print p_id, r_id, riboswitch
