from .. import riboswitch as rs
from .. import rna_f
from ..inverse_folding import two_target_local_refinement as loc_ref


OUTPUT_DIR = None
NUM_PARENT_GROUPS = 1000

rna_f.save_mode = 's'
rna_f.hd_mode = 'n'

rs.FUNCTIONAL_SITE_MAX_LENGTH = 50

HAIRPIN_LOOP_MIN_SIZE = 3
HAIRPIN_LOOP_MAX_SIZE = 5
HAIRPIN_STEM_MIN_SIZE = 7
HAIRPIN_STEM_MAX_SIZE = 20
UNBOUND_H_APT_DIFF = 4
EXPRESSION_PLATFORM_MAX_LEN = 50

loc_ref.OPTIMIZATION_CRITERIA = loc_ref.OptimizationCriteriaId.struct_prob_variance
loc_ref.NEIGHBOUR_SELECTION = loc_ref.NeighbourSelection.propability_dependent
loc_ref.SLS_ACCEPT_PROB = 0.025
loc_ref.MAX_STEPS = 10 * rs.FUNCTIONAL_SITE_MAX_LENGTH
loc_ref.MAX_SEARCH_TIME = 1.5
