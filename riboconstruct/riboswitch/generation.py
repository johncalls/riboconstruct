import collections
import os.path
import sys

from . import element as rs_e
from .. import riboswitch as rs
from .. import rna


BB_FILE_EXT = "bb"
BB_TERMINATOR = "<<<"
PATTERN_FILE_NAME = "pattern.p"


def write_to_file(file_, *out):
    """
    Writes a variable number of function arguments to *file_*.

    Raises an :exc:`~exceptions.IOError` if *file_* exists.
    """
    if os.path.isfile(file_):
        raise IOError("File '%s' already exists." % file_)
    with open(file_, 'w') as f:
        for line in out:
            f.write("%s\n" % str(line))


def write_instances(instance_iterator, output=sys.stdout):
    """
    Write all instances returned by the *instance_iterator* to *output*
    (standard is :attr:`sys.stdout`).
    """
    for p_id, r_id, r in instance_iterator:
        print >> output, p_id, r_id, r


class RiboswitchInstanceSpace(object):
    """
    Spans the space of :class:`~riboconstruct.riboswitch.Riboswitch`
    instances by defining target-action-constraint triplets which are
    used to generate new riboswitch instances from existing ones.

    On the basis of these triplets, :class:`RiboswitchInstanceSpace`
    provides different functionalities to validate and alter riboswitch
    instances, and to produce evaluation files necessary for the
    evaluation of riboswitch instances with ``RNAf``.

    See :func:`~ActionContainer.add` how the triplets have to be defined
    syntactically. Each riboswitch
    :class:`~riboconstruct.riboswitch.element` can be accessed within
    the constraints via :attr:`elements` using the specified
    identifiers. Note that the identified elements are only available on
    runtime.

    See the source of :class:`SpliceSiteInstanceSpace` for an exemplary
    use case.
    """

    __all__ = ["validate", "get_detailed_validation", "elements", "actions"
               "create_evaluation_files", "iter_alterations"]

    def __init__(self):
        self._current_instance = None
        self.elements = dict()
        """
        A :class:`dict` mapping identifiers to riboswitch
        :class:`~riboconstruct.riboswitch.element`\ s of the current
        :class:`~riboconstruct.riboswitch.Riboswitch` instance.
        """
        self.actions = ActionContainer()
        """
        Is the :class:`ActionContainer` storing the
        target-action-constraint triplets. The target is identified
        via the respective identifier in :attr:`elements`.
        """

    def validate(self, instance):
        """
        Validates the riboswitch *instance* according to the defined
        constraints.

        Overwrite this method in the specific
        :class:`RiboswitchInstanceSpace` subclass, if the validation
        defined in this method is not sufficient.
        """
        self._register(instance)
        if any(constraint(*c_args) if inverse else not constraint(*c_args)
               for constraint, c_args, _, inverse
               in self.actions.iter_all_constraints()):
            return False
        cf = self.elements[rs_e.ContextFront.ident]
        if cf.pos[0] != instance.pos[0]:
            return False
        if cf.pos[0] == cf.pos[1] and cf.pos[0] != instance.pos_riboswitch[0]:
            return False
        cb = self.elements[rs_e.ContextBack.ident]
        if cb.pos[1] != instance.pos[1]:
            return False
        if cb.pos[0] == cb.pos[1] and cb.pos[0] != instance.pos_riboswitch[1]:
            return False
        return True

    def get_detailed_validation(self, instance):
        """
        Returns a :class:`str` containing detailed information which
        of the defined constraints failed on the given riboswitch
        *instance*.
        """
        s = []
        self._register(instance)
        s.append("The following constraints failed:")
        for c in set(self.actions.iter_all_constraints()):
            constraint, c_args, descr, inverse = c
            result = (constraint(*c_args) if not inverse else
                      not constraint(*c_args))
            if not result:
                s.append('- %s' % descr)
        return '\n'.join(s)

    def create_evaluation_files(self, instance, folder):
        """
        Creates the evaluation files needed by ``RNAf`` to evaluate the
        riboswitch *instance*.

        The functionality has to be implemented in the specific
        :class:`RiboswitchInstanceSpace` subclass.
        """
        raise NotImplementedError

    def iter_alterations(self, instance):
        """
        Alters the riboswitch *instance* using the defined actions and
        iterates over the resulting instances. Only valid instances with
        respect to the constraints defined for each target-action pair
        are yielded.
        """
        def constraints_fulfilled(target_ident, action):
            return all(constraint(*c_args) if not inverse else
                       not constraint(*c_args)
                       for constraint, c_args, _, inverse
                       in self.actions.iter_constraints(target_ident, action))

        for target_ident in self.actions.iter_target_identifiers():
            for action_ident in self.actions.iter_actions(target_ident):
                self._register(instance, overwrite=True)
                alternation = instance.copy()
                action, a_args = action_ident
                # alter the current instance
                for single_target_ident in target_ident:
                    # perform the action and replace the old riboswitch
                    # element by the new one
                    old = self.elements[single_target_ident]
                    new = getattr(old, action)(*a_args)
                    alternation.replace(old, new)
                    self.elements[single_target_ident] = new
                    self._current_instance = None
                # check constraints for each target_ident and action
                if constraints_fulfilled(target_ident, action_ident):
                    yield alternation

    def _register(self, instance, overwrite=False):
        """
        Rearranges the elements of the riboswitch *instance* such that
        they are accessible by an unique identifier.

        Overwrite this method in the specific
        :class:`RiboswitchInstanceSpace` subclass, if the identifiers
        defined in this method are not sufficient.
        """
        if (not overwrite and
            self._current_instance is not None and
            instance == self._current_instance):
            return
        self._current_instance = instance
        self.elements = dict()
        # rearrange riboswitch elements by their identifiers
        for element in instance.elements:
            if element.state is None:
                ident = element.ident
            else:
                ident = "%s_%s" % (element.ident,
                                   rs_e.State.get_str(element.state))
            self.elements[ident] = element
        # if there are no context elements, add dummies
        if rs_e.ContextFront.ident not in self.elements:
            pos = (instance.pos_riboswitch[0], instance.pos_riboswitch[0])
            cf = rs_e.ContextFront(pos, rna.Sequence(""))
            self.elements[cf.ident] = cf
            instance.add(cf)
        if rs_e.ContextBack.ident not in self.elements:
            pos = (instance.pos_riboswitch[1], instance.pos_riboswitch[1])
            cb = rs_e.ContextBack(pos, rna.Sequence(""))
            self.elements[cb.ident] = cb
            instance.add(cb)

    def _within_dist(self, ident_a, ident_b, max_dist):
        """
        Returns whether there are at most *max_dist* base positions
        between riboswitch elements *a* and *b* in the current instance.
        """
        a = self.elements[ident_a]
        b = self.elements[ident_b]
        if a.pos[1] < b.pos[0]:  # check if a is downstream of b
            return b.pos[0] - a.pos[1] <= max_dist
        elif a.pos[0] > b.pos[1]:  # check if a is upstream of b
            return a.pos[0] - b.pos[1] <= max_dist
        # else: overlapping; do not handle this case
        return True

    def _overlapping(self, ident_a, ident_b):
        """
        Returns whether the riboswitch element *a* is overlapping with
        *b* in the current instance.
        """
        a = self.elements[ident_a]
        b = self.elements[ident_b]
        if a.pos[1] <= b.pos[0] or a.pos[0] >= b.pos[1]:
            return False
        return True

    def _matching(self, ident_a, ident_b):
        """
        Returns whether the structure of element *a* can be formed by
        sequence of element *b* in the current instance.
        """
        def seq_overlapping(struct_pos_i, struct_pos_j, seq_pos):
            # NOTE: seq_pos[1] has to be '>' since it's index is +1 (the 2nd
            #       index of riboswitch elements is always +1)
            return seq_pos[0] <= struct_pos_i and seq_pos[1] > struct_pos_j

        if not self._overlapping(ident_a, ident_b):
            return True
        a = self.elements[ident_a]
        b = self.elements[ident_b]
        a_start = a.pos[0]
        b_start = b.pos[0]
        for bp_id, (bp_pos_i, bp_pos_j) in enumerate(a.struct.bp_positions):
            bp_pos_i += a_start  # get 'global' index
            bp_pos_j += a_start  # get 'global' index
            # here it is assured that the structure and sequence overlap for
            # the current bp
            if not seq_overlapping(bp_pos_i, bp_pos_j, b.pos):
                return True
            # get the bases from b.seq that match the positions of the bp
            base_i = b.seq[bp_pos_i - b_start]
            base_j = b.seq[bp_pos_j - b_start]
            # check whether these bases can in fact form a valid bp
            if not rna.valid_bp(base_i, base_j):
                return False


class ActionContainer(object):
    """
    A container to manage target-action-constraint triplets. (See
    :func:`add` for a more detailed description.)
    """

    __all__ = ["add", "iter_target_identifiers", "iter_actions",
               "iter_constraints", "iter_all_constraints"]

    def __init__(self):
        self._container = dict()

    def add(self, target_ident, action, constraint, descr=None):
        """
        Add a target-action-constraint triplet to the container. These
        triplets define actions which are performed on the specified
        target if and only if the specified constraint is fulfilled.

        *target_ident* identifies the specific riboswitch elements the
        action is performed on within
        :attr:`~RiboswitchInstanceSpace.elements` at runtime. Has
        to be a tuple and can contain several identifiers for different
        riboswitch elements.

        *action* specifies how the targeted riboswitch element is
        altered. Has to be a tuple of the form ``(action_fct,
        action_arguments)``, where the first entry is a :class:`str`
        defining an alternation action and the second is a tuple
        specifying the arguments used for the action.

        *constraint* defines the constraint that has to be fulfilled to
        make the action valid. Has to be a tuple of the form
        ``(constraint_fct, constraint_args, constraint_descr,
        [inverse,])``:
            * ``constraint_fct`` is a function defining the constraint
            * ``constraint_args`` is a tuple specifying the
              arguments used for the constraint
            * ``constraint_descr`` is a short human readable description
              of the constraint.
            * ``inverse`` specifies whether the constraint should hold
              (``inverse = False``, standard value) or should not hold
              (``inverse = True``)

        Example: ::

            from riboconstruct.riboswitch import elements as rs_e
            from riboconstruct.riboswitch import generation as rs_gen

            actions = rs_gen.ActionContainer()

            h_ub = "%s_%s" % (rs_e.Hairpin.ident, "ub")
            h_b = "%s_%s" % (rs_e.Hairpin.ident, "b")

            target_ident = (h_ub,)
            action = ("shift", (-1,))
            # ensures that the bound hairpin is always left of the unbound
            # hairpin when shifting the unbound hairpin downstream
            constraint_fct = (lambda h_ub, h_b: elements[h_b].pos[0] <
                                                elements[h_ub].pos[0])
            constraint = (constraint_fct, (h_ub, h_b), "descr")
            actions.add(target_ident, action, constraint)

        """
        if len(constraint) == 3:
            constraint = constraint + (False,)
        try:
            self._container[target_ident][action].append(constraint)
        except KeyError:
            try:
                self._container[target_ident][action] = [constraint]
            except KeyError:
                self._container[target_ident] = dict()
                self._container[target_ident][action] = [constraint]

    def iter_target_identifiers(self):
        """Iterate all target identifiers."""
        return self._container.iterkeys()

    def iter_actions(self, target_ident):
        """Iterate all actions added for a specific target."""
        return self._container[target_ident].iterkeys()

    def iter_constraints(self, target_ident, action):
        """Iterate all constraints for a specific target-action pair."""
        return iter(self._container[target_ident][action])

    def iter_all_constraints(self):
        """
        Iterate all defined target-action pairs and yield the
        constraints defined for each pair.
        """
        for target_ident in self.iter_target_identifiers():
            for action in self.iter_actions(target_ident):
                for constraint in self.iter_constraints(target_ident, action):
                    yield constraint


class InstanceIterator(object):
    """
    Uses a :class:`RiboswitchInstanceSpace` and an initial
    :class:`~riboconstruct.riboswitch.Riboswitch` instance to iterate
    over all valid instances defined by the instance space.

    In each iteration step, the id of the parent the current riboswitch
    is derived from, the id of the current riboswitch, and the current
    riboswitch instance are returned.

    :func:`~RiboswitchInstanceSpace.iter_alterations` is used to
    alternate the current riboswitch instance in each iteration. Only
    instances which have not been generated before are returned.
    """

    __all__ = ["next"]

    def __init__(self, instance_space, initial_instance):
        self._instance_space = instance_space
        self._closed_list = set()
        self._open_list = collections.deque()
        self._id = 0

        self._open_list.append((-1, self._id, repr(initial_instance)))

    def __iter__(self):
        return self

    def next(self):
        """Iterates over the riboswitch instance space."""
        if not len(self._open_list):
           raise StopIteration("No more new riboswitches.")
        parent_id, riboswitch_id, riboswitch_str = self._open_list.popleft()
        riboswitch = rs.get_riboswitch_from_str(riboswitch_str)
        for sibling in self._instance_space.iter_alterations(riboswitch):
            if sibling in self._closed_list:
                continue
            self._id += 1
            self._open_list.append((riboswitch_id, self._id, repr(sibling)))
            self._closed_list.add(sibling)
        return parent_id, riboswitch_id, riboswitch


class SpliceSiteInstanceSpace(RiboswitchInstanceSpace):
    """
    Spans the space of :class:`~riboconstruct.riboswitch.Riboswitch`
    instances controlling the 5' splice site (cf. with the evaluation
    section in the master thesis) by defining target-action-constraint
    triplets which are used to generate new riboswitch instances from
    existing ones.
    """

    __all__ = ["validate", "get_detailed_validation", "elements", "actions"
               "create_evaluation_files", "iter_alterations"]

    def __init__(self,
                 hairpin_stem_size=(7, 20),
                 ub_hairpin_aptamer_dist=4,
                 expression_platform_max_len=50):
        """
        *hairpin_stem_size* is a tuple containing the minimum and
        maximum size the sequestor/anti-sequestor hairpin can have.

        *ub_hairpin_aptamer_dist* is the maximum distance between the
        anti-sequestor hairpin and the B1-2 region of the aptamer.

        *expression_platform_max_len* is the maximum length of the
        expression platform of the riboswitch.
        """
        super(SpliceSiteInstanceSpace, self).__init__()
        # specify target identifiers
        b_str  = rs_e.State.get_str(rs_e.State.bound)
        ub_str = rs_e.State.get_str(rs_e.State.unbound)
        self._h_b = "%s_%s" % (rs_e.Hairpin.ident, b_str)
        self._h_ub = "%s_%s" % (rs_e.Hairpin.ident, ub_str)
        self._a_b = "%s_%s" % (rs_e.Aptamer.ident, b_str)
        self._a_ub = "%s_%s" % (rs_e.Aptamer.ident, ub_str)
        self._ts = str(rs_e.TargetSite.ident)
        self._cf = str(rs_e.ContextFront.ident)
        self._cb = str(rs_e.ContextBack.ident)
        self._b1_2_ac = "b1_2_%s" % str(rs_e.AccessConstraint.ident)
        # ==============================================================
        # define actions and constraints
        # ==============================================================
        # all constraints are in the following format:
        #
        # self.actions.add(
        #   (self._h_b,),                                        # target_ident
        #   ("shift", (1,)),                                     # action
        #   (self._matching, (self._h_b, self._ts, 1), "descr"), # constraint
        # --------------------------------------------------------------
        # sequestor - shift up
        # --------------------------------------------------------------
        self.actions.add((self._h_b,),
                         ("shift", (1,)),
                         (self._matching, (self._h_b, self._ts),
                          "h_b, ts, formed bps are valid"))
        self.actions.add((self._h_b,),
                         ("shift", (1,)),
                         (self._overlapping, (self._h_b, self._ts),
                          "h_b, ts: overlapping"))
        self.actions.add((self._h_b,),
                         ("shift", (1,)),
                         (self._overlapping, (self._h_b, self._a_b),
                          "h_b, a_b: not overlapping", True))
        # --------------------------------------------------------------
        # anti-sequestor - shift up
        # --------------------------------------------------------------
        self.actions.add((self._h_ub,),
                          ("shift", (1,)),
                          (self._matching, (self._h_ub, self._ts),
                           "h_ub, ts: formed bps are valid"))
        self.actions.add((self._h_ub,),
                          ("shift", (1,)),
                          (self._overlapping, (self._h_ub, self._a_ub),
                           "h_ub, a_ub: not overlapping", True))
        self.actions.add((self._h_ub,),
                          ("shift", (1,)),
                          (self._within_dist, (self._h_ub, self._h_b, 0),
                           "h_ub, h_b: within distance"))
        self.actions.add((self._h_ub,),
                          ("shift", (1,)),
                          (self._within_dist, (self._h_ub, self._h_b, 0),
                           "h_ub, h_b: within distance"))
        # should not overlap with the B1-2 region
        self.actions.add((self._h_ub,),
                          ("shift", (1,)),
                          (self._overlapping, (self._h_ub, self._b1_2_ac),
                           "h_ub, b1_2_ac: not overlapping", True))
        # --------------------------------------------------------------
        # sequestor - shift down
        # --------------------------------------------------------------
        self.actions.add((self._h_b,),
                          ("shift", (-1,)),
                          (self._matching, (self._h_b, self._ts),
                           "h_b, ts: formed bps are valid"))
        self.actions.add((self._h_b,),
                          ("shift", (-1,)),
                          (self._within_dist, (self._h_b, self._h_ub, 0),
                           "h_b, h_ub: within distance"))
        # sequestor has to stay within expression platform
        self.actions.add((self._h_b,),
                          ("shift", (-1,)),
                          (lambda h_b, cf: self.elements[h_b].pos[0] >=
                                           self.elements[cf].pos[1],
                           (self._h_b, self._cf),
                           "h_b, cf: cf not to the right of h_b"))
        self.actions.add((self._h_b,),
                          ("shift", (-1,)),
                          (lambda h_b, h_ub: self.elements[h_b].pos[0] <
                                             self.elements[h_ub].pos[0],
                           (self._h_b, self._h_ub),
                           "h_b, h_ub: h_b to the left of h_ub"))
        # --------------------------------------------------------------
        # anti-sequestor - shift down
        # --------------------------------------------------------------
        self.actions.add((self._h_ub,),
                          ("shift", (-1,)),
                          (self._matching, (self._h_ub, self._ts),
                           "h_ub, ts: formed bps are valid"))
        self.actions.add((self._h_ub,),
                          ("shift", (-1,)),
                          (self._overlapping, (self._h_ub, self._ts),
                           "h_ub, ts: not overlapping", True))
        self.actions.add((self._h_ub,),
                          ("shift", (-1,)),
                          (self._within_dist,
                           (self._h_ub, self._a_b, ub_hairpin_aptamer_dist),
                           "h_ub, a_b: within distance"))
        self.actions.add((self._h_ub,),
                          ("shift", (-1,)),
                          (lambda h_ub, h_b: self.elements[h_b].pos[0] <
                                             self.elements[h_ub].pos[0],
                           (self._h_ub, self._h_b),
                           "h_ub, h_b: h_b to the right of h_ub"))
        # --------------------------------------------------------------
        # target site - shift up ---> decreases expression platform
        # --------------------------------------------------------------
        self.actions.add((self._ts, self._cf),
                          ("shift", (1,)),
                          (self._overlapping, (self._ts, self._h_ub),
                           "ts, h_ub: not overlapping", True))
        # context front (!) and hairpin should not overlap
        self.actions.add((self._ts, self._cf),
                          ("shift", (1,)),
                          (self._overlapping, (self._cf, self._h_b),
                           "cf, h_b: not overlapping", True))
        self.actions.add((self._ts, self._cf),
                          ("shift", (1,)),
                          (self._matching, (self._h_b, self._ts),
                           "ts, h_b: formed bps are valid"))
        # --------------------------------------------------------------
        # target site - shift down ---> increases expression platform
        # --------------------------------------------------------------
        # size of expression platform is limited
        self.actions.add((self._ts, self._cf),
                          ("shift", (-1,)),
                          (lambda ts, cf, a_ub:
                              (self.elements[a_ub].pos[0] -
                                 (self.elements[ts].pos[0])) <=
                              expression_platform_max_len,
                           (self._ts, self._cf, self._a_ub),
                           "size of expression platform within limits"))
        # target site and sequestor should still overlap
        self.actions.add((self._ts, self._cf),
                          ("shift", (-1,)),
                          (self._overlapping, (self._h_b, self._ts),
                           "ts, h_b: overlapping"))
        self.actions.add((self._ts, self._cf),
                          ("shift", (-1,)),
                          (self._matching, (self._h_b, self._ts),
                           "ts, h_b: formed bps are valid"))
        # --------------------------------------------------------------
        # sequestor - increase stem
        # --------------------------------------------------------------
        self.actions.add((self._h_b,),
                          ("insert_bp_before", (0,)),
                          (self._matching, (self._h_b, self._ts),
                           "h_b, ts: formed bps are valid"))
        self.actions.add((self._h_b,),
                          ("insert_bp_before", (0,)),
                          (self._overlapping, (self._h_b, self._a_b),
                           "h_b, a_b: not overlapping", True))
        # sequestor has to stay within expression platform
        self.actions.add((self._h_b,),
                          ("insert_bp_before", (0,)),
                          (self._overlapping, (self._h_b, self._cf),
                           "h_b, cf: not overlapping", True))
        # there is a maximum stem size
        self.actions.add((self._h_b,),
                          ("insert_bp_before", (0,)),
                          (lambda h_b:
                              len(self.elements[h_b].struct.bp_positions) <
                              hairpin_stem_size[1],
                           (self._h_b,), "h_b: stem size ok"))
        # --------------------------------------------------------------
        # anti-sequestor - increase stem
        # --------------------------------------------------------------
        self.actions.add((self._h_ub,),
                          ("insert_bp_before", (0,)),
                          (self._overlapping, (self._h_ub, self._a_ub),
                           "h_ub, a_ub: not overlapping", True))
        self.actions.add((self._h_ub,),
                          ("insert_bp_before", (0,)),
                          (self._overlapping, (self._h_ub, self._ts),
                           "h_ub, ts: not overlapping", True))
        self.actions.add((self._h_ub,),
                          ("insert_bp_before", (0,)),
                          (lambda h_ub, h_b: self.elements[h_b].pos[0] <
                                             self.elements[h_ub].pos[0],
                           (self._h_ub, self._h_b),
                           "h_b, h_ub: h_b to the left of h_ub"))
        # there is a maximum stem size
        self.actions.add((self._h_ub,),
                          ("insert_bp_before", (0,)),
                          (lambda h_ub:
                              len(self.elements[h_ub].struct.bp_positions) <
                              hairpin_stem_size[1],
                           (self._h_ub,), "h_ub: stem size ok"))
        # --------------------------------------------------------------
        # sequestor - decrease stem
        # --------------------------------------------------------------
        self.actions.add((self._h_b,),
                          ("remove_bp", (0,)),
                          (self._matching, (self._h_b, self._ts),
                           "h_b, ts: formed bps are valid"))
        self.actions.add((self._h_b,),
                          ("remove_bp", (0,)),
                          (lambda h_b, h_ub: self.elements[h_b].pos[0] <
                                             self.elements[h_ub].pos[0],
                           (self._h_b, self._h_ub),
                           "h_b, h_ub: h_b to the left of h_ub"))
        self.actions.add((self._h_b,),
                          ("remove_bp", (0,)),
                          (self._within_dist, (self._h_b, self._h_ub, 0),
                           "h_b, h_ub: within distance"))
        # there is a minimum stem size
        self.actions.add((self._h_b,),
                          ("remove_bp", (0,)),
                          (lambda h_b:
                             len(self.elements[h_b].struct.bp_positions) >=
                             hairpin_stem_size[0],
                           (self._h_b,), "h_b: stem size ok"))
        self.actions.add((self._h_b,),
                          ("remove_bp", (0,)),
                          (self._overlapping, (self._h_b, self._ts),
                           "h_b, ts: overlapping"))
        # --------------------------------------------------------------
        # anti-sequestor - decrease stem
        # --------------------------------------------------------------
        self.actions.add((self._h_ub,),
                          ("remove_bp", (0,)),
                          (self._within_dist, (self._h_ub, self._h_b, 0),
                           "h_ub, h_b: within distance"))
        # there is a minimum stem size
        self.actions.add((self._h_ub,),
                          ("remove_bp", (0,)),
                          (lambda h_ub:
                              len(self.elements[h_ub].struct.bp_positions) >=
                              hairpin_stem_size[0],
                           (self._h_ub,), "h_ub: stem size ok"))
        self.actions.add((self._h_ub,),
                          ("remove_bp", (0,)),
                          (self._within_dist,
                           (self._h_ub, self._a_b, ub_hairpin_aptamer_dist),
                           "h_ub, a_b: within distance"))

    def create_evaluation_files(self, instance, folder):
        """
        Creates the evaluation files needed by ``RNAf`` in
        *folder* to evaluate the riboswitch *instance*.
        """

        self._register(instance)
        if not os.path.exists(folder):
            raise OSError("Folder '%s' does not exist." % folder)
        if not self.validate(instance):
            raise ValueError("Riboswitch instance not valid.")

        b = rs_e.State.bound
        ub = rs_e.State.unbound
        b_str = rs_e.State.get_str(b)
        ub_str = rs_e.State.get_str(ub)

        folders = ["", ""]
        folders[b] = os.path.join(folder, b_str)
        folders[ub] = os.path.join(folder, ub_str)

        structs, seq = instance.get_constraints()
        b_struct = structs[b]
        ub_struct = structs[ub]

        pos = instance.pos
        pos_riboswitch = instance.pos_riboswitch
        start = pos[0]

        aptamer_start_b = self.elements[self._a_b].pos[0]
        aptamer_start_ub = self.elements[self._a_ub].pos[0]

        # --------------------------------------------------------------
        # write sequestor building block (bb)
        # --------------------------------------------------------------
        i = pos_riboswitch[0] - start
        j = aptamer_start_b - start
        write_to_file(os.path.join(folders[b],
                                   "%s.%s" % (self._h_b, BB_FILE_EXT)),
                      ">%s" % self._h_b,
                      seq[i:j],
                      b_struct[i:j],
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write anti-sequestor bb
        # --------------------------------------------------------------
        j = aptamer_start_ub - start
        write_to_file(os.path.join(folders[ub],
                                   "%s.%s" % (self._h_ub, BB_FILE_EXT)),
                      ">%s" % self._h_ub,
                      seq[i:j],
                      ub_struct[i:j],
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write bound aptamer bb
        # --------------------------------------------------------------
        i = aptamer_start_b - start
        j = pos_riboswitch[1] - start
        write_to_file(os.path.join(folders[b],
                                   "%s.%s" % (self._a_b, BB_FILE_EXT)),
                      ">%s" % self._a_b,
                      seq[i:j],
                      b_struct[i:j],
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write unbound aptamer bb
        # --------------------------------------------------------------
        i = aptamer_start_ub - start
        write_to_file(os.path.join(folders[ub],
                                   "%s.%s" % (self._a_ub, BB_FILE_EXT)),
                      ">%s" % self._a_ub,
                      seq[i:j],
                      ub_struct[i:j],
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write target site bb
        # --------------------------------------------------------------
        ts = self.elements[self._ts]
        # NOTE: +1 because of different index counting
        # HACK: +3 omits the ATG/AUG at the start
        i = ts.pos[0] - pos_riboswitch[0] + 1 + 3
        j = ts.pos[1] - pos_riboswitch[0]
        pattern_ts = "%c(%i,%i)" % (ord('%'), i, j)
        write_to_file(os.path.join(folders[0],
                                   "%s.%s" % (self._ts, BB_FILE_EXT)),
                      ">%s" % self._ts,
                      str(ts.seq),
                      BB_TERMINATOR)
        write_to_file(os.path.join(folders[1],
                                   "%s.%s" % (self._ts, BB_FILE_EXT)),
                      ">%s" % self._ts,
                      str(ts.seq),
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write accessibility constraint bb
        # --------------------------------------------------------------
        ac = self.elements[self._b1_2_ac]
        # NOTE: +1 because of different index counting
        i = ac.pos[0] - self.elements[self._a_ub].pos[0] + 1
        j = ac.pos[1] - self.elements[self._a_ub].pos[1]
        pattern_ac = "%c(%i,%i)" % (ord('%'), i, j)
        write_to_file(os.path.join(folders[ub],
                                   "%s.%s" % (self._b1_2_ac, BB_FILE_EXT)),
                      ">%s" % self._b1_2_ac,
                      str(ac.seq),
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write context front bb (if existent)
        # --------------------------------------------------------------
        if pos_riboswitch[0] != pos[0]:
            i = 0
            j = pos_riboswitch[0]
            pattern_cf = "%s%s" % (self._cf, unichr(167).encode("utf-8"))
            seq = self.elements[self._cf].seq
            write_to_file(os.path.join(folders[0],
                                       "%s.%s" % (self._cf, BB_FILE_EXT)),
                          ">%s" % self._cf,
                          seq[i:j],
                          BB_TERMINATOR)
            write_to_file(os.path.join(folders[1],
                                       "%s.%s" % (self._cf, BB_FILE_EXT)),
                          ">%s" % self._cf,
                          seq[i:j],
                          BB_TERMINATOR)
        # --------------------------------------------------------------
        # write context back bb (if existent)
        # --------------------------------------------------------------
        if pos_riboswitch[1] != pos[1]:
            i = pos_riboswitch[1]
            j = pos[1]
            pattern_cb = "%s%s" % (self._cf, unichr(167).encode("utf-8"))
            seq = self.elements[self._cb].seq
            write_to_file(os.path.join(folders[0],
                                       "%s.%s" % (self._cb, BB_FILE_EXT)),
                          ">%s" % self._cb,
                          seq[i:j],
                          BB_TERMINATOR)
            write_to_file(os.path.join(folders[1],
                                       "%s.%s" % (self._cb, BB_FILE_EXT)),
                          ">%s" % self._cb,
                          seq[i:j],
                          BB_TERMINATOR)
        # --------------------------------------------------------------
        # write pattern file
        # --------------------------------------------------------------
        evaluation_patterns = [None, None]
        # bound case
        evaluation_patterns[rs_e.State.bound] = (
            "({}|%%s%%s#%%s%s&%%s%%s|1)" % unichr(167).encode("utf-8"))
        # unbound case
        evaluation_patterns[rs_e.State.unbound] = (
            "({}|%%s%%s#%%s%s%%s#%%s%%s|1)" % unichr(167).encode("utf-8"))
        # fill the patterns with the identifiers and write them to file
        write_to_file(os.path.join(folders[b], PATTERN_FILE_NAME),
                      ">pattern_%s" % b_str,
                      evaluation_patterns[b] % (pattern_cf, pattern_ts,
                                                self._h_b, self._a_b,
                                                pattern_cb))
        write_to_file(os.path.join(folders[ub], PATTERN_FILE_NAME),
                      ">pattern_%s" % ub_str,
                      evaluation_patterns[ub] % (pattern_cf, pattern_ts,
                                                 self._h_ub, pattern_ac,
                                                 self._a_ub, pattern_cb))

        return folders

    def _register(self, instance, overwrite=False):
        super(SpliceSiteInstanceSpace, self)._register(instance, overwrite)
        try:
            access_constraint = self.elements.pop(rs_e.AccessConstraint.ident)
            self.elements[self._b1_2_ac] = access_constraint
        except KeyError:
            pass
