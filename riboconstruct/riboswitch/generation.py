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
                   for constraint, c_args, inverse
                   in self._actions.iter_all_constraints())

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
                ident = str(element.ident)
            else:
                ident = "%s_%s" % (element.ident, element.state)
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
                       for constraint, c_args, inverse
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
        return a.pos[1] > b.pos[0] or a.pos[0] < b.pos[1]

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

    def add(self, target_ident, action, constraint):
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
        """
        if len(constraint) == 2:
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
        """
        Iterate all constraints added for a specific target and action.
        """
        return iter(self._container[target_ident][action])

    def iter_all_constraints(self):
        """Iterate all defined constraints."""
        for target_ident in self.iter_target_identifiers():
            for action in self.iter_actions(target_ident):
                for constraint in self.iter_constraints(target_ident, action):
                    yield constraint


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
        for sibling in self._instance_space.iter_alterations(riboswitch):
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
        # store the evaluation pattern
        self._evaluation_patterns = [None, None]
        # bound case
        self._evaluation_patterns[rs_e.State.bound] = (
            "({}|%%s%%s#%%s%s&%%s%%s|1)" % unichr(167).encode("utf-8"))
        # unbound case
        self._evaluation_patterns[rs_e.State.unbound] = (
            "({}|%%s%%s#%%s%s%%s#%%s%%s|1)" % unichr(167).encode("utf-8"))
        self._evaluation_patterns = tuple(self._evaluation_patterns)
        # specify target identifiers
        h_b = "%s_%s" % (rs_e.Hairpin.ident, rs_e.State.bound)
        h_ub = "%s_%s" % (rs_e.Hairpin.ident, rs_e.State.unbound)
        a_b = "%s_%s" % (rs_e.Aptamer.ident, rs_e.State.bound)
        a_ub = "%s_%s" % (rs_e.Aptamer.ident, rs_e.State.unbound)
        ts = str(rs_e.TargetSite.ident)
        cf = str(rs_e.ContextFront.ident)
        # cb = str(rs_e.ContextBack.ident)
        b1_2_ac = "%s_%s" % ("b1_2", str(rs_e.AccessConstraint.ident))
        # ==============================================================
        # define actions and constraints
        # ==============================================================
        # all constraints are in the following format:
        #
        # self._actions.add(
        #   (h_b,),                          # target_ident
        #   ("shift", (1,)),      # action
        #   (self._matching, (h_b, ts, 1)))  # constraint
        # --------------------------------------------------------------
        # sequestor - shift up
        # --------------------------------------------------------------
        self._actions.add((h_b,),
                         ("shift", (1,)),
                         (self._matching, (h_b, ts)))
        self._actions.add((h_b,),
                         ("shift", (1,)),
                         (self._overlapping, (h_b, ts)))
        self._actions.add((h_b,),
                         ("shift", (1,)),
                         (self._overlapping, (h_b, a_b), True))
        self._actions.add((h_b,),
                         ("shift", (-1,)),
                         (lambda h_b, h_ub: self.element(h_b).pos[0] <
                                            self.element(h_ub).pos[0],
                          (h_b, h_ub)))
        # --------------------------------------------------------------
        # anti-sequestor - shift up
        # --------------------------------------------------------------
        self._actions.add((h_ub,),
                         ("shift", (1,)),
                         (self._matching, (h_ub, ts)))
        self._actions.add((h_ub,),
                         ("shift", (1,)),
                         (self._overlapping, (h_ub, a_ub), True))
        self._actions.add((h_ub,),
                         ("shift", (1,)),
                         (self._within_dist, (h_ub, h_b, 0)))
        self._actions.add((h_ub,),
                         ("shift", (1,)),
                         (self._within_dist, (h_ub, h_b, 0)))
        # should not overlap with the B1-2 region
        self._actions.add((h_ub,),
                         ("shift", (1,)),
                         (self._overlapping, (h_ub, b1_2_ac), True))
        # --------------------------------------------------------------
        # sequestor - shift down
        # --------------------------------------------------------------
        self._actions.add((h_b,),
                         ("shift", (-1,)),
                         (self._matching, (h_b, ts)))
        self._actions.add((h_b,),
                         ("shift", (-1,)),
                         (self._within_dist, (h_b, h_ub, 0)))
        # sequestor has to stay within expression platform
        self._actions.add((h_b,),
                         ("shift", (-1,)),
                         (lambda h_b, cf: self.element(h_b).pos[0] >=
                                          self.element(cf).pos[1],
                          (h_b, cf)))
        # --------------------------------------------------------------
        # anti-sequestor - shift down
        # --------------------------------------------------------------
        self._actions.add((h_ub,),
                         ("shift", (-1,)),
                         (self._matching, (h_ub, ts)))
        self._actions.add((h_ub,),
                         ("shift", (-1,)),
                         (self._overlapping, (h_ub, ts), True))
        self._actions.add((h_ub,),
                         ("shift", (-1,)),
                         (self._within_dist,
                          (h_ub, a_b, ub_hairpin_aptamer_dist)))
        self._actions.add((h_ub,),
                         ("shift", (-1,)),
                         (lambda h_ub, h_b: self.element(h_b).pos[0] <
                                            self.element(h_ub).pos[0],
                          (h_ub, h_b)))
        # --------------------------------------------------------------
        # target site - shift up ---> decreases expression platform
        # --------------------------------------------------------------
        self._actions.add((ts, cf),
                         ("shift", (1,)),
                         (self._overlapping, (ts, h_ub), True))
        # context front (!) and hairpin should not overlap
        self._actions.add((ts, cf),
                         ("shift", (1,)),
                         (self._overlapping, (cf, h_b), True))
        self._actions.add((ts, cf),
                         ("shift", (1,)),
                         (self._matching, (h_b, ts)))
        # --------------------------------------------------------------
        # target site - shift down ---> increases expression platform
        # --------------------------------------------------------------
        # size of expression platform is limited
        self._actions.add((ts, cf),
                         ("shift", (-1,)),
                         (lambda ts, cf, a_ub:
                             (self.element(a_ub).pos[0] -
                                (self.element(ts).pos[0] - 1)) <=
                             expression_platform_max_len,
                          (ts, cf, a_ub)))
        # target site and sequestor should still overlap
        self._actions.add((ts, cf),
                         ("shift", (-1,)),
                         (self._overlapping, (h_b, ts)))
        self._actions.add((ts, cf),
                         ("shift", (-1,)),
                         (self._matching, (h_b, ts)))
        # --------------------------------------------------------------
        # sequestor - increase stem
        # --------------------------------------------------------------
        self._actions.add((h_b,),
                         ("insert_bp_before", (0,)),
                         (self._matching, (h_b, ts)))
        self._actions.add((h_b,),
                         ("insert_bp_before", (0,)),
                         (self._overlapping, (h_b, a_b), True))
        # sequestor has to stay within expression platform
        self._actions.add((h_b,),
                         ("insert_bp_before", (0,)),
                         (self._overlapping, (h_b, cf), True))
        # there is a maximum stem size
        self._actions.add((h_b,),
                         ("insert_bp_before", (0,)),
                         (lambda h_b:
                             len(self.element(h_b).struct.bp_positions) <
                             hairpin_stem_size[1],
                          (h_b,)))
        # --------------------------------------------------------------
        # anti-sequestor - increase stem
        # --------------------------------------------------------------
        self._actions.add((h_ub,),
                         ("insert_bp_before", (0,)),
                         (self._overlapping, (h_ub, a_ub), True))
        self._actions.add((h_ub,),
                         ("insert_bp_before", (0,)),
                         (self._overlapping, (h_ub, ts), True))
        self._actions.add((h_ub,),
                         ("insert_bp_before", (0,)),
                         (lambda h_ub, h_b: self.element(h_b).pos[0] <
                                            self.element(h_ub).pos[0],
                          (h_ub, h_b)))
        # there is a maximum stem size
        self._actions.add((h_ub,),
                         ("insert_bp_before", (0,)),
                         (lambda h_ub:
                             len(h_ub.struct.bp_positions) <
                             hairpin_stem_size[1],
                          (h_ub,)))
        # --------------------------------------------------------------
        # sequestor - decrease stem
        # --------------------------------------------------------------
        self._actions.add((h_b,),
                         ("remove_bp", (0,)),
                         (self._matching, (h_b, ts)))
        self._actions.add((h_b,),
                         ("remove_bp", (0,)),
                         (lambda h_b, h_ub: self.element(h_b).pos[0] <
                                            self.element(h_ub).pos[0],
                          (h_b, h_ub)))
        self._actions.add((h_b,),
                         ("remove_bp", (0,)),
                         (self._within_dist, (h_b, h_ub, 0)))
        # there is a minimum stem size
        self._actions.add((h_b,),
                         ("remove_bp", (0,)),
                         (lambda h_b:
                            len(self.element(h_b).struct.bp_positions) >
                            hairpin_stem_size[0],
                          (h_b,)))
        # --------------------------------------------------------------
        # anti-sequestor - decrease stem
        # --------------------------------------------------------------
        self._actions.add((h_ub,),
                         ("remove_bp", (0,)),
                         (self._within_dist, (h_ub, h_b, 0)))
        # there is a minimum stem size
        self._actions.add((h_b,),
                         ("remove_bp", (0,)),
                         (lambda h_b:
                             len(self.element(h_b).struct.bp_positions) >
                             hairpin_stem_size[0],
                          (h_b,)))
        self._actions.add((h_ub,),
                         ("remove_bp", (0,)),
                         (self._within_dist,
                          (h_ub, a_b, ub_hairpin_aptamer_dist)))

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
        folders[b] = os.path.join(folder, rs_e.State.get_str(b))
        folders[ub] = os.path.join(folder, rs_e.State.get_str(ub))
        folders = tuple(folders)

        structs, seq = instance.get_constraints()
        b_struct = structs[b]
        ub_struct = structs[ub]

        pos_instance = instance.pos
        pos_riboswitch = instance.pos_riboswitch
        start = pos_instance[0]

        ident_a_b = "%s_%s" % (rs_e.Aptamer.ident, b_str)
        ident_a_ub = "%s_%s" % (rs_e.Aptamer.ident, ub_str)
        aptamer_start_b = self._elements[ident_a_b].pos[0]
        aptamer_start_ub = self._elements[ident_a_ub].pos[0]

        # --------------------------------------------------------------
        # write sequestor building block (bb)
        # --------------------------------------------------------------
        ident = "%s_%s" % (rs_e.Hairpin.ident, b_str)
        i = pos_riboswitch[0] - start
        j = aptamer_start_b - start
        patterns_h_b = ident
        write_to_file(os.path.join(folders[b], "%s.%s" % (ident, BB_FILE_EXT)),
                      ">%s" % ident,
                      seq[i:j],
                      b_struct[i:j],
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write anti-sequestor bb
        # --------------------------------------------------------------
        ident = "%s_%s" % (rs_e.Hairpin.ident, ub_str)
        j = aptamer_start_ub - start
        patterns_h_ub = ident
        write_to_file(os.path.join(folders[ub],
                                   "%s.%s" % (ident, BB_FILE_EXT)),
                      ">%s" % ident,
                      seq[i:j],
                      ub_struct[i:j],
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write bound aptamer bb
        # --------------------------------------------------------------
        i = aptamer_start_b - start
        j = pos_riboswitch[1] - start
        patterns_a_b = ident_a_b
        write_to_file(os.path.join(folders[b],
                                   "%s.%s" % (ident_a_b, BB_FILE_EXT)),
                      ">%s" % ident_a_b,
                      seq[i:j],
                      b_struct[i:j],
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write unbound aptamer bb
        # --------------------------------------------------------------
        i = aptamer_start_ub - start
        patterns_a_ub = ident_a_ub
        write_to_file(os.path.join(folders[ub],
                                   "%s.%s" % (ident_a_ub, BB_FILE_EXT)),
                      ">%s" % ident_a_ub,
                      seq[i:j],
                      ub_struct[i:j],
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write target site bb
        # --------------------------------------------------------------
        ident = rs_e.TargetSite.ident
        ts = self._elements[str(ident)]
        # NOTE: +1 because of different index counting
        # HACK: +3 omits the ATG/AUG at the start
        i = ts.pos[0] - pos_riboswitch[0] + 1 + 3
        j = ts.pos[1] - pos_riboswitch[0]
        pattern_ts = "%c(%i,%i)" % (ord('%'), i, j)
        write_to_file(os.path.join(folders[0],
                                   "%s.%s" % (ident, BB_FILE_EXT)),
                      ">%s" % ident,
                      str(ts.seq),
                      BB_TERMINATOR)
        write_to_file(os.path.join(folders[1],
                                   "%s.%s" % (ident, BB_FILE_EXT)),
                      ">%s" % ident,
                      str(ts.seq),
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write accessibility constraint bb
        # --------------------------------------------------------------
        ident = "%s_%s" % ("b1_2", str(rs_e.AccessConstraint.ident))
        ac = self._elements[ident]
        a_ub = "%s_%s" % (self._elements[rs_e.Aptamer.ident], ub_str)
        # NOTE: +1 because of different index counting
        i = ac.pos[0] - a_ub.pos[0] + 1
        j = ac.pos[1] - a_ub.pos[1]
        pattern_ac = "%c(%i,%i)" % (ord('%'), i, j)
        write_to_file(os.path.join(folders[ub],
                                   "%s.%s" % (ident, BB_FILE_EXT)),
                      ">%s" % ident,
                      str(ac.seq),
                      BB_TERMINATOR)
        # --------------------------------------------------------------
        # write context front bb (if existent)
        # --------------------------------------------------------------
        if pos_riboswitch[0] != pos_instance[0]:
            ident = rs_e.ContextFront.ident
            i = 0
            j = pos_riboswitch[0]
            pattern_cf = "%s%s" % (ident, unichr(167).encode("utf-8"))
            write_to_file(os.path.join(folders[0],
                                       "%s.%s" % (ident, BB_FILE_EXT)),
                          ">%s" % ident,
                          seq[i:j],
                          BB_TERMINATOR)
            write_to_file(os.path.join(folders[1],
                                       "%s.%s" % (ident, BB_FILE_EXT)),
                          ">%s" % ident,
                          seq[i:j],
                          BB_TERMINATOR)
        # --------------------------------------------------------------
        # write context back bb (if existent)
        # --------------------------------------------------------------
        if pos_riboswitch[1] != pos_instance[1]:
            ident = rs_e.ContextBack.ident
            i = pos_riboswitch[1]
            j = pos_instance[1]
            pattern_cb = "%s%s" % (ident, unichr(167).encode("utf-8"))
            write_to_file(os.path.join(folders[0],
                                       "%s.%s" % (ident, BB_FILE_EXT)),
                          ">%s" % ident,
                          seq[i:j],
                          BB_TERMINATOR)
            write_to_file(os.path.join(folders[1],
                                       "%s.%s" % (ident, BB_FILE_EXT)),
                          ">%s" % ident,
                          seq[i:j],
                          BB_TERMINATOR)
        # --------------------------------------------------------------
        # write pattern file
        # --------------------------------------------------------------
        write_to_file(os.path.join(folders[b], PATTERN_FILE_NAME),
                      ">pattern_%s" % b_str,
                      self._evaluation_patterns[b] % (pattern_cf,
                                                      pattern_ts,
                                                      patterns_h_b,
                                                      patterns_a_b,
                                                      pattern_cb))
        write_to_file(os.path.join(folders[ub], PATTERN_FILE_NAME),
                      ">pattern_%s" % ub_str,
                      self._evaluation_patterns[ub] % (pattern_cf,
                                                       pattern_ts,
                                                       patterns_h_ub,
                                                       pattern_ac,
                                                       patterns_a_ub,
                                                       pattern_cb))

        return folders
