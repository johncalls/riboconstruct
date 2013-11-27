import collections
import os.path

from . import element as rs_e
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


class InstanceSpace(object):
    """
    Defines the boundaries an altered riboswitch instance has to stay
    within.

    The current instance that is going to be altered is registered
    using :func:`register` and the resulting instances are iterated
    using :func:`iter_alterations`. The rules how to alter the current
    instance and which constraints they have to fulfill are stored in an
    :class:`ActionContainer`. They have to be defined in the specific
    subclasses of :class:`InstanceSpace`.
    """

    __all__ = ["validate", "register", "create_evaluation_files",
               "iter_alterations"]

    def __init__(self):
        self._current_instance = None
        self._elements = dict()
        self._actions = ActionContainer()

    def validate(self, instance):
        """
        Validates the riboswitch *instance* according to the defined
        constraints.

        Overwrite this method in the specific :class:`InstanceSpace`
        subclass, if the validation defined in this method is not
        sufficient.
        """
        self.register(instance)
        return all(constraint(*c_args) if not inverse else
                   not constraint(*c_args)
                   for (constraint, c_args, inverse), _
                   in self._actions.iter_all_constraints())

    def print_validation(self, instance):
        print "The following constraints failed:"
        for c, descr in set(self._actions.iter_all_constraints()):
            constraint, c_args, inverse = c
            result = (constraint(*c_args) if not inverse else
                      not constraint(*c_args))
            if not result:
                print '-', descr

    def register(self, instance):
        """
        Rearranges the elements of the riboswitch *instance* such that
        they are accessible by an unique identifier.

        Overwrite this method in the specific :class:`InstanceSpace`
        subclass, if the identifiers defined in this method are not
        sufficient.
        """
        if (self._current_instance is not None and
            instance == self._current_instance):
            return
        self._current_instance = instance
        self._elements = dict()
        # rearrange riboswitch elements by their identifiers
        for element in instance.elements:
            if element.state is None:
                ident = element.ident
            else:
                ident = "%s_%s" % (element.ident,
                                   rs_e.State.get_str(element.state))
            self._elements[ident] = element
        # if there are no context elements, add dummies
        if rs_e.ContextFront.ident not in instance.elements:
            pos = (instance.pos_riboswitch[0], instance.pos_riboswitch[0])
            cf = rs_e.ContextFront(pos, ())
            self._elements[cf.ident] = cf
        if rs_e.ContextBack.ident not in instance.elements:
            pos = (instance.pos_riboswitch[1], instance.pos_riboswitch[1])
            cb = rs_e.ContextBack(pos, ())
            self._elements[cb.ident] = cb

    def create_evaluation_files(self, instance, folder):
        """
        Creates the evaluation files needed by ``RNAf`` to evaluate the
        riboswitch *instance*.

        The functionality has to be defined in the specific
        :class:`InstanceSpace` subclass.
        """
        raise NotImplementedError

    def iter_alterations(self, instance):
        """
        Alters the riboswitch *instance* using the defined actions and
        iterates over the resulting instances. Only valid instances with
        respect to the constraints defined for each action are returned.
        """
        def constraints_fulfilled(target_ident, action):
            return all(constraint(*c_args) if not inverse else
                       not constraint(*c_args)
                       for (constraint, c_args, inverse), _
                       in self._actions.iter_constraints(target_ident, action))

        self.register(instance)
        for target_ident in self._actions.iter_target_identifiers():
            for action_ident in self._actions.iter_actions(target_ident):
                action, a_args = action_ident
                # alter the current instance
                alternation = instance.copy()
                for single_target_ident in target_ident:
                    # perform the action and replace the old riboswitch
                    # element by the new one
                    old = self._elements[single_target_ident]
                    new = getattr(old, action)(*a_args)
                    alternation.replace(old, new)
                # check constraints for each target_ident and action
                if constraints_fulfilled(target_ident, action_ident):
                    yield alternation

    def _within_dist(self, ident_a, ident_b, max_dist):
        """
        Returns whether there are at most *max_dist* base positions
        between riboswitch elements *a* and *b* in the current instance.
        """
        a = self._elements[ident_a]
        b = self._elements[ident_b]
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
        a = self._elements[ident_a]
        b = self._elements[ident_b]
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
        a = self._elements[ident_a]
        b = self._elements[ident_b]
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

    def element(self, ident):
        return self._elements[ident]


class ActionContainer(object):
    """
    A Container to manage actions, the riboswitch elements they target
    at to modify them and the constraints they have to fulfill.
    """

    __all__ = ["add", "iter_target_identifiers", "iter_actions",
               "iter_constraints", "iter_all_constraints"]

    def __init__(self):
        self._container = dict()

    def add(self, target_ident, action, constraint, descr=None):
        """
        Add an action which is only performed for a specific target if
        the specified constraint is fulfilled to the container.

        *target_ident* identifies the specific riboswitch element(s) at
        runtime for which the action is executed for. Has to be a tuple
        and can contain several identifiers for different riboswitch
        elements.

        *action* specifies how the targeted riboswitch element is
        altered. Has to be a tuple of the form ``(action_fct,
        action_arguments)``, where the latter is a tuple specifying the
        arguments used for the action.

        *constraint* defines the constraint that has to be fulfilled to
        make the action valid. Has to be a tuple of the form
        ``(constraint_fct, constraint_arguments, inverse)``.
        ``constraint_arguments`` is a tuple specifying the arguments
        used for the constraint, whereas *inverse* specifies whether the
        constraint has to be met (``inverse = False``, standard value)
        or not (``inverse = True``). Several constraints can be added
        for the same target and action, and each is checked before the
        action is performed.

        *descr* is a short human readable description of each
        target-action-constraint triple.
        """
        if len(constraint) == 2:
            constraint = constraint + (False,)
        try:
            self._container[target_ident][action].append((constraint, descr))
        except KeyError:
            try:
                self._container[target_ident][action] = [(constraint, descr)]
            except KeyError:
                self._container[target_ident] = dict()
                self._container[target_ident][action] = [(constraint, descr)]

    def iter_target_identifiers(self):
        """Iterate all target identifiers."""
        return self._container.iterkeys()

    def iter_actions(self, target_ident):
        """Iterate all actions added for a specific target."""
        return self._container[target_ident].iterkeys()

    def iter_constraints(self, target_ident, action):
        """
        Iterate all constraint-description pairs for a specific target
        and action.
        """
        return iter(self._container[target_ident][action])

    def iter_all_constraints(self):
        """
        Iterate all defined target-action-constraint triples and yield
        the constraint and description for each triple."""
        for target_ident in self.iter_target_identifiers():
            for action in self.iter_actions(target_ident):
                for constraint, descr in self.iter_constraints(target_ident,
                                                               action):
                    yield constraint, descr


class InstanceIterator(object):
    """
    Uses a riboswitch instance space object and an initial riboswitch
    instance to iterate over all valid instances within the constraints
    defined by the instance space.

    In each iteration step, the id of the parent the current riboswitch
    is derived from, the id of the current riboswitch and the current
    riboswitch are returned.
    """

    __all__ = ["next"]

    def __init__(self, instance_space, initial_instance):
        self._instance_space = instance_space
        self._closed_list = set()
        self._open_list = collections.deque()
        self._id = 0

        self._open_list.append((-1, initial_instance))

    def __iter__(self):
        return self

    def next(self):
        """Defines how to iterate over all possible riboswitches."""
        if not len(self._open_list):
           raise StopIteration("No more new riboswitches.")
        parent_id, riboswitch = self._open_list.popleft()
        riboswitch_id = self._id
        for i, sibling in enumerate(self._instance_space.iter_alterations(riboswitch)):
            if sibling in self._closed_list:
                continue
            self._id += 1
            self._closed_list.add(sibling)
            self._open_list.append((riboswitch_id, sibling))
        return parent_id, riboswitch_id, riboswitch


class SpliceSiteInstanceSpace(InstanceSpace):
    """
    Defines the space of riboswitch instances controlling the 5' splice
    site (cf. with the evaluation section in the master thesis).
    """

    __all__ = ["iter_alterations", "validate", "create_evaluation_files",
               "register"]

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
        # self._actions.add(
        #   (self._h_b,),                               # target_ident
        #   ("shift", (1,)),                           # action
        #   (self._matching, (self._h_b, self._ts, 1)),  # constraint
        #   "descr")                                   # descr
        # --------------------------------------------------------------
        # sequestor - shift up
        # --------------------------------------------------------------
        self._actions.add((self._h_b,),
                          ("shift", (1,)),
                          (self._matching, (self._h_b, self._ts)),
                          "h_b, ts, formed bps are valid")
        self._actions.add((self._h_b,),
                          ("shift", (1,)),
                          (self._overlapping, (self._h_b, self._ts)),
                          "h_b, ts: overlapping")
        self._actions.add((self._h_b,),
                          ("shift", (1,)),
                          (self._overlapping, (self._h_b, self._a_b), True),
                          "h_b, a_b: not overlapping")
        # --------------------------------------------------------------
        # anti-sequestor - shift up
        # --------------------------------------------------------------
        self._actions.add((self._h_ub,),
                          ("shift", (1,)),
                          (self._matching, (self._h_ub, self._ts)),
                          "h_ub, ts: formed bps are valid")
        self._actions.add((self._h_ub,),
                          ("shift", (1,)),
                          (self._overlapping, (self._h_ub, self._a_ub), True),
                          "h_ub, a_ub: not overlapping")
        self._actions.add((self._h_ub,),
                          ("shift", (1,)),
                          (self._within_dist, (self._h_ub, self._h_b, 0)),
                          "h_ub, h_b: within distance")
        self._actions.add((self._h_ub,),
                          ("shift", (1,)),
                          (self._within_dist, (self._h_ub, self._h_b, 0)),
                          "h_ub, h_b: within distance")
        # should not overlap with the B1-2 region
        self._actions.add((self._h_ub,),
                          ("shift", (1,)),
                          (self._overlapping, (self._h_ub, self._b1_2_ac),
                           True),
                          "h_ub, b1_2_ac: not overlapping")
        # --------------------------------------------------------------
        # sequestor - shift down
        # --------------------------------------------------------------
        self._actions.add((self._h_b,),
                          ("shift", (-1,)),
                          (self._matching, (self._h_b, self._ts)),
                          "h_b, ts: formed bps are valid")
        self._actions.add((self._h_b,),
                          ("shift", (-1,)),
                          (self._within_dist, (self._h_b, self._h_ub, 0)),
                          "h_b, h_ub: within distance")
        # sequestor has to stay within expression platform
        self._actions.add((self._h_b,),
                          ("shift", (-1,)),
                          (lambda h_b, cf: self.element(h_b).pos[0] >=
                                           self.element(cf).pos[1],
                           (self._h_b, self._cf)),
                          "h_b, cf: cf not to the right of h_b")
        self._actions.add((self._h_b,),
                          ("shift", (-1,)),
                          (lambda h_b, h_ub: self.element(h_b).pos[0] <
                                             self.element(h_ub).pos[0],
                           (self._h_b, self._h_ub)),
                          "h_b, h_ub: h_b to the left of h_ub")
        # --------------------------------------------------------------
        # anti-sequestor - shift down
        # --------------------------------------------------------------
        self._actions.add((self._h_ub,),
                          ("shift", (-1,)),
                          (self._matching, (self._h_ub, self._ts)),
                          "h_ub, ts: formed bps are valid")
        self._actions.add((self._h_ub,),
                          ("shift", (-1,)),
                          (self._overlapping, (self._h_ub, self._ts), True),
                          "h_ub, ts: not overlapping")
        self._actions.add((self._h_ub,),
                          ("shift", (-1,)),
                          (self._within_dist,
                           (self._h_ub, self._a_b, ub_hairpin_aptamer_dist)),
                          "h_ub, a_b: within distance")
        self._actions.add((self._h_ub,),
                          ("shift", (-1,)),
                          (lambda h_ub, h_b: self.element(h_b).pos[0] <
                                             self.element(h_ub).pos[0],
                           (self._h_ub, self._h_b)),
                          "h_ub, h_b: h_b to the right of h_ub")
        # --------------------------------------------------------------
        # target site - shift up ---> decreases expression platform
        # --------------------------------------------------------------
        self._actions.add((self._ts, self._cf),
                          ("shift", (1,)),
                          (self._overlapping, (self._ts, self._h_ub), True),
                          "ts, h_ub: not overlapping")
        # context front (!) and hairpin should not overlap
        self._actions.add((self._ts, self._cf),
                          ("shift", (1,)),
                          (self._overlapping, (self._cf, self._h_b), True),
                          "cf, h_b: not overlapping")
        self._actions.add((self._ts, self._cf),
                          ("shift", (1,)),
                          (self._matching, (self._h_b, self._ts)),
                          "ts, h_b: formed bps are valid")
        # --------------------------------------------------------------
        # target site - shift down ---> increases expression platform
        # --------------------------------------------------------------
        # size of expression platform is limited
        self._actions.add((self._ts, self._cf),
                          ("shift", (-1,)),
                          (lambda ts, cf, a_ub:
                              (self.element(a_ub).pos[0] -
                                 (self.element(ts).pos[0])) <=
                              expression_platform_max_len,
                           (self._ts, self._cf, self._a_ub)),
                          ("size of expression platform within limits"))
        # target site and sequestor should still overlap
        self._actions.add((self._ts, self._cf),
                          ("shift", (-1,)),
                          (self._overlapping, (self._h_b, self._ts)),
                          "ts, h_b: overlapping")
        self._actions.add((self._ts, self._cf),
                          ("shift", (-1,)),
                          (self._matching, (self._h_b, self._ts)),
                          "ts, h_b: formed bps are valid")
        # --------------------------------------------------------------
        # sequestor - increase stem
        # --------------------------------------------------------------
        self._actions.add((self._h_b,),
                          ("insert_bp_before", (0,)),
                          (self._matching, (self._h_b, self._ts)),
                          "h_b, ts: formed bps are valid")
        self._actions.add((self._h_b,),
                          ("insert_bp_before", (0,)),
                          (self._overlapping, (self._h_b, self._a_b), True),
                          "h_b, a_b: not overlapping")
        # sequestor has to stay within expression platform
        self._actions.add((self._h_b,),
                          ("insert_bp_before", (0,)),
                          (self._overlapping, (self._h_b, self._cf), True),
                          "h_b, cf: not overlapping")
        # there is a maximum stem size
        self._actions.add((self._h_b,),
                          ("insert_bp_before", (0,)),
                          (lambda h_b:
                              len(self.element(h_b).struct.bp_positions) <
                              hairpin_stem_size[1],
                           (self._h_b,)),
                          "h_b: stem size ok")
        # --------------------------------------------------------------
        # anti-sequestor - increase stem
        # --------------------------------------------------------------
        self._actions.add((self._h_ub,),
                          ("insert_bp_before", (0,)),
                          (self._overlapping, (self._h_ub, self._a_ub), True),
                          "h_ub, a_ub: not overlapping")
        self._actions.add((self._h_ub,),
                          ("insert_bp_before", (0,)),
                          (self._overlapping, (self._h_ub, self._ts), True),
                          "h_ub, ts: not overlapping")
        self._actions.add((self._h_ub,),
                          ("insert_bp_before", (0,)),
                          (lambda h_ub, h_b: self.element(h_b).pos[0] <
                                             self.element(h_ub).pos[0],
                           (self._h_ub, self._h_b)),
                          "h_b, h_ub: h_b to the left of h_ub")
        # there is a maximum stem size
        self._actions.add((self._h_ub,),
                          ("insert_bp_before", (0,)),
                          (lambda h_ub:
                              len(self.element(h_ub).struct.bp_positions) <
                              hairpin_stem_size[1],
                           (self._h_ub,)),
                          "h_ub: stem size ok")
        # --------------------------------------------------------------
        # sequestor - decrease stem
        # --------------------------------------------------------------
        self._actions.add((self._h_b,),
                          ("remove_bp", (0,)),
                          (self._matching, (self._h_b, self._ts)),
                          "h_b, ts: formed bps are valid")
        self._actions.add((self._h_b,),
                          ("remove_bp", (0,)),
                          (lambda h_b, h_ub: self.element(h_b).pos[0] <
                                             self.element(h_ub).pos[0],
                           (self._h_b, self._h_ub)),
                          "h_b, h_ub: h_b to the left of h_ub")
        self._actions.add((self._h_b,),
                          ("remove_bp", (0,)),
                          (self._within_dist, (self._h_b, self._h_ub, 0)),
                          "h_b, h_ub: within distance")
        # there is a minimum stem size
        self._actions.add((self._h_b,),
                          ("remove_bp", (0,)),
                          (lambda h_b:
                             len(self.element(h_b).struct.bp_positions) >=
                             hairpin_stem_size[0],
                           (self._h_b,)),
                          "h_b: stem size ok")
        # --------------------------------------------------------------
        # anti-sequestor - decrease stem
        # --------------------------------------------------------------
        self._actions.add((self._h_ub,),
                          ("remove_bp", (0,)),
                          (self._within_dist, (self._h_ub, self._h_b, 0)),
                          "h_ub, h_b: formed bps are valid")
        # there is a minimum stem size
        self._actions.add((self._h_ub,),
                          ("remove_bp", (0,)),
                          (lambda h_ub:
                              len(self.element(h_ub).struct.bp_positions) >=
                              hairpin_stem_size[0],
                           (self._h_ub,)),
                          "h_ub: stem size ok")
        self._actions.add((self._h_ub,),
                          ("remove_bp", (0,)),
                          (self._within_dist, (self._h_ub, self._a_b,
                                               ub_hairpin_aptamer_dist)),
                          "h_ub, a_b: within distance")

    def register(self, instance):
        super(SpliceSiteInstanceSpace, self).register(instance)
        try:
            access_constraint = self._elements.pop(rs_e.AccessConstraint.ident)
            self._elements[self._b1_2_ac] = access_constraint
        except KeyError:
            pass

    def create_evaluation_files(self, instance, folder):
        self.register(instance)
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

        aptamer_start_b = self._elements[self._a_b].pos[0]
        aptamer_start_ub = self._elements[self._a_ub].pos[0]

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
        ts = self._elements[self._ts]
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
        ac = self._elements[self._b1_2_ac]
        # NOTE: +1 because of different index counting
        i = ac.pos[0] - self._elements[self._a_ub].pos[0] + 1
        j = ac.pos[1] - self._elements[self._a_ub].pos[1]
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
            seq = self._elements[self._cf].seq
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
            seq = self._elements[self._cb].seq
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
