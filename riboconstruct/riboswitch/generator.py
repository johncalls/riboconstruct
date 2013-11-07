import collections

from . import riboswitch as rs
from . import element as rs_e

from .. import rna


EXPRESSION_PLATFORM_MAX_LEN = 50

HAIRPIN_LOOP_MAX_SIZE = 7
HAIRPIN_LOOP_MIN_SIZE = 3
HAIRPIN_STEM_MIN_SIZE = 5
HAIRPIN_STEM_MAX_SIZE = 30


class Generator(object):
    """
    Basic riboswitch iterator. Based on the initial riboswitch the whole
    riboswitch space within the defined constraints is iterated. In each step
    a new :class:`riboconstruct.riboswitch.riboswitch.Riboswitch` is returned
    and used to generate new riboswitches. Thereby, riboswitches already
    generated are discarded.

    How to generate a new riboswitch based on the current riboswitch is
    specificied in the specific subclasses by defining the set of actions used
    to alter the single riboswitch elements and the set of constraints that
    have to be fulfilled by each riboswitch.
    """

    def __init__(self, initial_riboswitch):
        self.validate(initial_riboswitch)
        # for the riboswitch generation
        self.open_list = collections.deque()
        self.closed_list = set()
        self.counter = 0
        # for the enforcement of constraints
        self.elements = dict()
        self.actions = ActionContainer()
        # add the initial riboswitch
        self.open_list.append((-1, self.counter, str(initial_riboswitch)))

    def __iter__(self):
        return self

    def next(self):
        parent_id, riboswitch_id, riboswitch_str = self.open_list.popleft()
        riboswitch = rs.get_riboswitch_from_str(riboswitch_str)
        self.rearrange_riboswitch_elements(riboswitch)
        for sibling in self.iter_siblings(riboswitch):
            if str(sibling) in self.closed_list:
                continue
            self.counter += 1
            self.open_list.append((riboswitch_id, self.counter, str(sibling)))
            self.closed_list.put(str(sibling))
        return parent_id, riboswitch_id, riboswitch

    def iter_siblings(self, riboswitch):
        """
        Returns an iterator over all sibling riboswitches based on alternations
        of the current one.

        The current riboswitch is represented by its elements stored in
        ``Generator.elements``. The elements are specificied in the respective
        subclasses via :func:`rearrange_riboswitch_elements`.
        """
        def constraints_fulfilled(target_ident, action):
            return all(constraint(*c_args) if not inverse else
                       not constraint(*c_args)
                       for constraint, c_args, inverse
                       in self.actions.iter_constraints(target_ident, action))

        for target_ident in self.actions.iter_target_identifiers():
            for action, a_args in self.actions.iter_actions(target_ident):
                # check constraints for each target_ident and action
                if constraints_fulfilled(target_ident, action):
                    # generate sibling; alter this sibling in the following
                    sibling = riboswitch.copy()
                    for single_target_ident in target_ident:
                        # perform the action and replace the old riboswitch
                        # element by the new one
                        old = self.elements[single_target_ident]
                        new = action(old, *a_args)
                        sibling.replace(old, new)
                    yield sibling

    def within_dist(self, ident_a, ident_b, max_dist, num_pos=1):
        """
        Returns whether there are at most *max_dist* positions between
        riboswitch elements *a* and *b* after *a* has been shifted upstream by
        *num_pos* positions, if ``num_pos > 0``. If ``num_pos < 0``, the shift
        is downstream.
        """
        a = self.elements[ident_a]
        b = self.elements[ident_b]
        new_a_pos_i = a.pos[0] + num_pos
        new_a_pos_j = a.pos[1] + num_pos
        if new_a_pos_j < b.pos[0]:  # check if a is downstream of b
            return b.pos[0] - new_a_pos_j <= max_dist
        elif new_a_pos_i > b.pos[1]:  # check if a is upstream of b
            return new_a_pos_i - b.pos[1] <= max_dist
        # else: overlapping; do not handle this case
        return True

    def overlapping(self, ident_a, ident_b, num_pos=1):
        """
        Returns whether the riboswitch element *a* is overlapping with *b*
        after *a* has been shifted upstream by *num_pos* positions, if
        ``num_pos > 0``. If ``num_pos < 0``, the shift is downstream.
        """
        a = self.elements[ident_a]
        b = self.elements[ident_b]
        return a.pos[1] + num_pos > b.pos[0] or a.pos[0] - num_pos < b.pos[1]

    def matching_shift(self, ident_a, ident_b, num_pos=1):
        """
        Checks whether the structure of element *a* can be formed by sequence
        of element *b* after *a* has been shifted upstream by *num_pos*
        positions, if ``num_pos > 0``. If ``num_pos < 0``, the shift is
        downstream.
        """
        def seq_overlapping(struct_pos_i, struct_pos_j, seq_pos):
            return seq_pos[0] <= struct_pos_i and seq_pos[1] >= struct_pos_j

        a = self.elements[ident_a]
        b = self.elements[ident_b]
        if not self.overlapping(a, b):
            return True
        a_i = a.pos[0]
        b_i = b.pos[0]
        for bp_id, (bp_pos_i, bp_pos_j) in enumerate(a.struct.bp_positions):
            bp_pos_i += a_i + num_pos
            bp_pos_j += a_i + num_pos
            # here it is assured that in each step b_start is <= bp_pos_i
            if not seq_overlapping(bp_pos_i, bp_pos_j, b.pos):
                return True
            base_i = b.seq[bp_pos_i - b_i + num_pos]
            base_j = b.seq[bp_pos_j - b_i + num_pos]
            if not rna.valid_bp(base_i, base_j):
                return False


class ActionContainer(object):
    """
    Container to manage action, the targets they aim to modify and the
    constraints they have to fulfill.
    """

    def __init__(self):
        self.container = dict()

    def add(self, target_ident, action, constraint):
        """
        Add an action which is only performed for a specific target if the
        specified constraint is fulfilled to the container.

        *target_ident* identifies the specific riboswitch element(s) at runtime
        for which the action is executed for. Has to be a tuple and can contain
        several identifiers for different riboswitch elements.

        *action* specifies how the targeted riboswitch element is altered. Has
        to be a tuple of the form ``(action_fct, action_arguments)``, where
        the latter is a tuple specifying the arguments used for the action.

        *constraint* defines the constraint that has to be fulfilled before the
        action is executed. Has to be a tuple of the form
        ``(constraint_fct, constraint_arguments, inverse)``.
        *constraint_arguments* is a tuple specifying the arguments used for the
        constraint, whereas *inverse* specifies whether the constraint has to
        be met (``inverse = False``, standard value) or not
        (``inverse = True``). Several constraints can be added for the same
        target and action, and each is checked before the action is performed.
        """
        try:
            self.container[target_ident][action].append(constraint)
        except KeyError:
            self.container[target_ident] = dict()
            self.container[target_ident][action] = [constraint]

    def iter_target_identifiers(self):
        """Iterate all target identifiers."""
        return self.container.iterkeys()

    def iter_actions(self, target_ident):
        """Iterate all actions added for a specific target."""
        return self.container[target_ident].iterkeys()

    def iter_constraints(self, target_ident, action):
        """Iterate all constraints added for a specific target and action."""
        return self.container[target_ident][action].itervalues()


class SpliceSiteRiboswitchGenerator(Generator):
    """
    Generates riboswitches controlling the 5' splice site (cf. with evaluation
    in the master thesis) by defining actions to iterate over all possible
    riboswitches within specificied constraints.
    """

    def __init__(self,
                 initial_riboswitch,
                 hairpin_stem_size=(7, 20),
                 ub_hairpin_aptamer_dist=4,
                 expression_platform_max_len=50):
        super(SpliceSiteRiboswitchGenerator, self).__init__(initial_riboswitch)
        # store the parameters
        self.hairpin_stem_size = hairpin_stem_size
        self.ub_hairpin_aptamer_dist = ub_hairpin_aptamer_dist
        self.expression_platform_max_len = expression_platform_max_len
        # specify target identifiers
        h_b = "%s_%s" % (rs_e.Hairpin.ident, rs_e.State.bound)
        h_ub = "%s_%s" % (rs_e.Hairpin.ident, rs_e.State.unbound)
        a_b = "%s_%s" % (rs_e.Aptamer.ident, rs_e.State.bound)
        a_ub = "%s_%s" % (rs_e.Aptamer.ident, rs_e.State.unbound)
        ts = str(rs_e.TargetSite.ident)
        cf = str(rs_e.ContextFront.ident)
        cb = str(rs_e.ContextBack.ident)
        b1_2_ac = "%s_%s" % ("b1_2", str(rs_e.AccessConstraint.ident))
        # define actions and constraints
        # =====================================================================
        # all constraints are in the following format:
        #
        # self.actions.add((h_b,),                         # target_ident
        #                  (rs_e.Hairpin.shift, (1,)),     # action
        #                  (self.matching_shift, (h_b, ts, 1)))  # constraint
        # =====================================================================
        # ----------------------------------------------------------------------
        # sequestor - shift up
        # ----------------------------------------------------------------------
        self.actions.add((h_b,),
                         (rs_e.Hairpin.shift, (1,)),
                         (self.matching_shift, (h_b, ts, 1)))
        self.actions.add((h_b,),
                         (rs_e.Hairpin.shift, (1,)),
                         (self.overlapping, (h_b, ts, 1)))
        self.actions.add((h_b,),
                         (rs_e.Hairpin.shift, (1,)),
                         (self.overlapping, (h_b, a_b, 1), True))
        # constraint: sequestor.pos[0] < anti-sequestor.pos[1]
        # --> indirectly enforced by anti-sequestor.pos[0] > target_site.pos[1]
        # and target_site.pos[0] >= sequestor.pos[0]
        # self.actions.add((h_b,),
        #                  (rs_e.Hairpin.shift, (1,)),
        #                  (lambda h_b, h_ub: h_b.pos[0] + 1 <= h_ub.pos[0],
        #                   (h_b, h_ub)))
        # ----------------------------------------------------------------------
        # anti-sequestor - shift up
        # ----------------------------------------------------------------------
        self.actions.add((h_ub,),
                         (rs_e.Hairpin.shift, (1,)),
                         (self.matching_shift, (h_ub, ts, 1)))
        self.actions.add((h_ub,),
                         (rs_e.Hairpin.shift, (1,)),
                         (self.overlapping, (h_ub, a_ub, 1), True))
        self.actions.add((h_ub,),
                         (rs_e.Hairpin.shift, (1,)),
                         (self.within_dist, (h_ub, h_b, 0)))
        self.actions.add((h_ub,),
                         (rs_e.Hairpin.shift, (1,)),
                         (self.within_dist, (h_ub, h_b, 0, 1)))
        # should not overlap with the B1-2 region
        self.actions.add((h_ub,),
                         (rs_e.Hairpin.shift, (1,)),
                         (self.overlapping, (h_ub, b1_2_ac, 1), True))
        # ----------------------------------------------------------------------
        # sequestor - shift down
        # ----------------------------------------------------------------------
        self.actions.add((h_b,),
                         (rs_e.Hairpin.shift, (-1,)),
                         (self.matching_shift, (h_b, ts, -1)))
        self.actions.add((h_b,),
                         (rs_e.Hairpin.shift, (-1,)),
                         (self.within_dist, (h_b, h_ub, 0, -1)))
        # sequestor has to stay within expression platform
        self.actions.add((h_b,),
                         (rs_e.Hairpin.shift, (-1,)),
                         (lambda h_b, cf: h_b.pos[0] - 1 >= cf.pos[1],
                          (h_b, cf)))
        # ----------------------------------------------------------------------
        # anti-sequestor - shift down
        # ----------------------------------------------------------------------
        self.actions.add((h_ub,),
                         (rs_e.Hairpin.shift, (-1,)),
                         (self.matching_shift, (h_ub, ts, -1)))
        self.actions.add((h_ub,),
                         (rs_e.Hairpin.shift, (-1,)),
                         (self.overlapping, (h_ub, ts, -1), True))
        self.actions.add((h_ub,),
                         (rs_e.Hairpin.shift, (-1,)),
                         (self.within_dist,
                          (h_ub, a_b, self.ub_hairpin_aptamer_dist, -1)))
        # constraint: sequestor.pos[0] < anti-sequestor.pos[1]
        # --> indirectly enforced by anti-sequestor.pos[0] > target_site.pos[1]
        # and target_site.pos[0] >= sequestor.pos[0]
        # self.actions.add((h_ub,),
        #                  (rs_e.Hairpin.shift, (-1,)),
        #                  (lambda h_ub, h_b: h_b.pos[0] <= h_ub.pos[0] - 1,
        #                   (h_ub, h_b)))
        # ----------------------------------------------------------------------
        # target site - shift up
        # ----------------------------------------------------------------------
        self.actions.add((ts, cf),
                         (rs_e.TargetSite.shift, (1,)),
                         (self.overlapping, (ts, h_ub, 1), True))
        # context front (!) and hairpin should not overlap
        self.actions.add((ts, cf),
                         (rs_e.TargetSite.shift, (1,)),
                         (self.overlapping, (cf, h_b, 1), True))
        self.actions.add((ts, cf),
                         (rs_e.TargetSite.shift, (1,)),
                         (self.matching_shift, (h_b, ts, -1)))  # shift ts up eqv. to shift h down
        self.actions.add((ts, cf),
                         (rs_e.TargetSite.shift, (1,)),
                         (self.matching_shift, (h_ub, ts, -1)))  # shift ts up eqv. to shift h down
        # ----------------------------------------------------------------------
        # target site - shift down
        # ----------------------------------------------------------------------
        # size of expression platform is limited
        self.actions.add((ts, cf),
                         (rs_e.TargetSite.shift, (-1,)),
                         (lambda ts, cf, a_ub:
                             a_ub.pos[0] - (ts.pos[0] - 1) <=
                             self.expression_platform_max_len,
                          (ts, cf, a_ub)))
        # target site and sequestor should still overlap
        self.actions.add((ts, cf),
                         (rs_e.TargetSite.shift, (-1,)),
                         (self.overlapping, (h_b, ts, 1)))  # shift ts down eqv. to shift h up
        self.actions.add((ts, cf),
                         (rs_e.TargetSite.shift, (-1,)),
                         (self.matching_shift, (h_b, ts, 1)))  # shift ts down eqv. to shift h up
        self.actions.add((ts, cf),
                         (rs_e.TargetSite.shift, (-1,)),
                         (self.matching_shift, (h_ub, ts, 1)))  # shift ts down eqv. to shift h up
        # ----------------------------------------------------------------------
        # sequestor - increase stem
        # ----------------------------------------------------------------------
        self.actions.add((h_b,),
                         (rs_e.Hairpin.insert_bp_before, (0,)),
                         (self.matching_shift, (h_b, ts, 1)))

    def rearrange_riboswitch_elements(self, riboswitch):
        # rearrange riboswitch elements by their identifiers
        for element in riboswitch.elements:
            if element.state is None:
                ident = str(element.ident)
            else:
                ident = "%s_%s" % (element.ident, element.state)
            self.elements[ident] = element

    @classmethod
    def validate(self):
        """
        Validate the current riboswitch according to the defined constraints.
        """
        pass
