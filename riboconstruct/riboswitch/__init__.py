from . import element as rs_e

from .. import rna


def iter_riboswitch(config):
    """
    Yields the parsed riboswitch elements that are encoded in *config*
    in a FASTA-like format.

    Example riboswitch in FASTA-like format: ::

        >Aptamer|38,103|bound
        ((((((((((........(((((....)))))..(((((...........)))))))))))))))
        NNNNNNNCCUAAAACAUACCAGAGAAAUCUGGAGAGGUGAAGAAUACGACCACCUAGGNNNNNNN
        >Aptamer|50,103|unbound
        ......(((((....)))))..(((((...........)))))..........
        AACAUACCAGAGAAAUCUGGAGAGGUGAAGAAUACGACCACCUAGGNNNNNNN
        >Access_constraint|50,56
        AACAUA
        >Hairpin|15,32|unbound
        (((((((...)))))))
        >Hairpin|0,17|bound
        (((((((...)))))))
        >Restriction_site|0,10
        AUGAGUAUGU
        >Context_front|-45,0
        AAGCUAUACCAAGCAUACAAUCAACUCCAAGCUAGAUCUCUUAAG
        >Context_back|103,148
        AUCUAGCGCUGGUACCAUCCCAUUUAACUGUAAGAAGAAUUGCAC

    Read file: ::

        >>> from riboconstruct.riboswitch import iter_riboswitch
        >>> with open("example_file") as config:
        ...     for element in iter_riboswitch(config):
        ...             print element

    :func:`get_riboswitch_from_config_file` can be used to read a
    FASTA-like riboswitch file and return a :class:`Riboswitch`.
    """
    def next():
        seq = config_iter.next().strip()
        while seq[0] == '#':
            seq = config_iter.next().strip()
        if seq[0] == '>':
            raise ValueError("Config has wrong format.")
        return seq

    config_iter = iter(config)
    re_type = None
    for line in config_iter:
        line = line.strip()
        if not len(line):
            continue
        if line[0] == '>' and not re_type:
            line = line[1:].split('|')
            if len(line) < 2:
                raise ValueError("Config has wrong format.")
            re_type = getattr(rs_e.Type, line[0].strip().lower())
            pos = line[1].split(',')
            pos = int(pos[0]), int(pos[1])
            if (re_type == rs_e.Type.aptamer or
                re_type == rs_e.Type.hairpin):
                # catch these two exception to raise a ValueError which is
                # better suited
                try:
                    state = getattr(rs_e.State, line[2].strip().lower())
                except (AttributeError, IndexError):
                    raise ValueError("Config has wrong format.")
        elif line[0] == '#' or line[0] == '':
            continue
        elif re_type is not None:
            if line[0] == '>':
                raise ValueError('Config has wrong format.')
            if re_type == rs_e.Type.aptamer:
                struct = rna.Structure(line)
                seq = rna.IUPACSequence(next().upper())
                if pos[1] - pos[0] != len(struct) != len(seq):
                    raise ValueError("Aptamer element has wrong size.")
                yield rs_e.Aptamer(state, pos, struct, seq)
            elif re_type == rs_e.Type.hairpin:
                struct = rna.Structure(line)
                if pos[1] - pos[0] != len(struct):
                    raise ValueError("Hairpin element has wrong size.")
                yield rs_e.Hairpin(state, pos, struct)
            elif re_type == rs_e.Type.target_site:
                seq = rna.IUPACSequence(line.upper())
                if pos[1] - pos[0] != len(seq):
                    raise ValueError("Restriction site element has wrong size.")
                yield rs_e.TargetSite(pos, seq)
            elif re_type == rs_e.Type.context_front:
                seq = rna.IUPACSequence(line)
                if pos[1] - pos[0] != len(seq):
                    raise ValueError("Context front element has wrong size.")
                yield rs_e.ContextFront(pos, seq)
            elif re_type == rs_e.Type.context_back:
                seq = rna.IUPACSequence(line)
                if pos[1] - pos[0] != len(seq):
                    raise ValueError("Context back element has wrong size.")
                yield rs_e.ContextBack(pos, seq)
            elif re_type == rs_e.Type.seq_constraint:
                seq = rna.IUPACSequence(line.upper())
                if pos[1] - pos[0] != len(seq):
                    raise ValueError("Sequence constraint element has wrong "
                                     "size.")
                yield rs_e.Sequence(pos, seq)
            elif re_type == rs_e.Type.access_constraint:
                yield rs_e.AccessConstraint(pos)
            else:
                raise TypeError("Riboswitch config file contains invalid "
                                "riboswitch element.")
            re_type = None
        else:
            raise ValueError("Config has wrong format.")


def get_riboswitch_from_config_file(config_file):
    """
    Returns a riboswitch instance by parsing the elements from a
    FASTA-like format using :func:`iter_riboswitch`.
    """
    riboswitch = Riboswitch()
    # read initial settings from config
    with open(config_file) as config:
        map(riboswitch.add, iter_riboswitch(config))
    # check if the rs elements 'fit' together
    riboswitch.validate()
    return riboswitch


import copy
import itertools
import os
import sys

from .. import params


EXPRESSION_PLATFORM_MAX_LEN = 50

PATTERNS = [None, None]

# in the bound case:
# due to ligand binding the aptamer is in a certain conformation forcing the
# hairpin to fold into a certain structure, i.e. use the aptamer as structural
# constraint and check if the hairpin folds into the modelled structure
PATTERNS[rs_e.State.bound] = (
    "({}|%%s%%s#%%s%s&%%s%%s|1)" % unichr(167).encode("utf-8"))

# in the unbound case:
# both hairpin and aptamer should fold into the modelled structure, i.e. use NO
# structural constraint and check if the aptamer and the hairpin fold into the
# modelled structures
PATTERNS[rs_e.State.unbound] = (
    "({}|%%s%%s#%%s%s%%s#%%s%%s|1)" % unichr(167).encode("utf-8"))

PATTERNS = tuple(PATTERNS)


def get_riboswitch_from_str(riboswitch_str):
    """
    Parse a riboswitch string representation and return the
    corresponding riboswitch instance.
    """
    riboswitch = Riboswitch()
    elements_iter = iter(riboswitch_str.split(';'))

    # first two elements are the positions
    # positions = elements_iter.next().split(', ')
    # riboswitch.pos_instance = (int(positions[0][8:]), int(positions[1][:-1]))
    positions = elements_iter.next().split(', ')
    riboswitch.pos_riboswitch = [int(positions[0][8:]), int(positions[1][:-1])]

    for element in elements_iter:
        element = element[1:-1].split(' ')  # get rid of '(' and ')' and split
        type_ = element[0].rsplit('_', 1)  # split type from state
        try:
            type_, state = type_
        except ValueError:  # element has no state
            type_ = type_[0]
            pos = int(element[1][1:-1]), int(element[2][:-1])
            seq = rna.Sequence(element[3][1:-1])
            if type_ == rs_e.TargetSite.ident:
                riboswitch.add(rs_e.TargetSite(pos, seq))
            elif type_ == rs_e.ContextFront.ident:
                riboswitch.add(rs_e.ContextFront(pos, seq))
            elif type_ == rs_e.ContextBack.ident:
                riboswitch.add(rs_e.ContextBack(pos, seq))
            elif type_ == rs_e.AccessConstraint.ident:
                riboswitch.add(rs_e.AccessConstraint(pos, seq))
        else:  # element has a state
            state = (rs_e.State.bound
                     if state == rs_e.State.get_str(rs_e.State.bound)
                     else rs_e.State.unbound)
            pos = int(element[1][1:-1]), int(element[2][:-1])
            struct = rna.Structure(element[3][1:-1])
            if type_ == rs_e.Aptamer.ident:
                seq = rna.Sequence(element[4][1:-1])
                riboswitch.add(rs_e.Aptamer(state, pos, struct, seq))
            elif type_ == rs_e.Hairpin.ident:
                riboswitch.add(rs_e.Hairpin(state, pos, struct))
    return riboswitch


def get_constraints((start, end), riboswitch_elements):
    def set_base(pos, new_base):
        if new_base == rna.BASES[rna.BaseId.UNSPEC]:
            return
        if seq[pos] == rna.BASES[rna.BaseId.UNSPEC]:
            seq[pos] = new_base
        # TODO: remove later since this should not happen
        elif seq[pos] != new_base:
            raise ValueError("Overlapping bases at position %i not valid %s "
                             "vs. %s" % (pos, seq[pos], new_base))

    # set new_struct_elem if the elem has not been specified before
    def set_struct_elem(state, pos, new_struct_elem):
        if structs[state][pos] == rna.STRUCT_ELEM_UNSPEC:
            structs[state][pos] = new_struct_elem

    length = end - start
    seq = [rna.BASES[rna.BaseId.UNSPEC]] * length
    structs = ([rna.STRUCT_ELEM_UNSPEC] * length,
               [rna.STRUCT_ELEM_UNSPEC] * length)

    for rs_elem in riboswitch_elements:
        # clip rs_elem position into the start-end range
        rs_elem_start = rs_elem.pos[0] if rs_elem.pos[0] > start else start
        rs_elem_end = rs_elem.pos[1] if rs_elem.pos[1] < end else end
        i, j = rs_elem_start - start, rs_elem_end - start

        if (rs_elem.type == rs_e.Type.target_site or
            rs_elem.type == rs_e.Type.seq_constraint):
            for pos in xrange(i, j):
                set_base(pos, rs_elem.seq[pos - i])
        elif rs_elem.type == rs_e.Type.hairpin:
            for pos in xrange(i, j):
                set_struct_elem(rs_elem.state, pos, rs_elem.struct[pos - i])
        elif rs_elem.type == rs_e.Type.aptamer:
            for pos in xrange(i, j):
                set_base(pos, rs_elem.seq[pos - i])
                set_struct_elem(rs_elem.state, pos, rs_elem.struct[pos - i])
        elif rs_elem.type == rs_e.Type.access_constraint:
            # The sequence of the access constraint should fit to the aptamer
            # sequence it is overlapping with (has been checked before), i.e.
            # do nothing here
            pass
        else:  # context_front or context_back
            for pos in xrange(i, j):
                set_base(pos, rs_elem.seq[pos - i])

    return (''.join(structs[0]), ''.join(structs[1])), ''.join(seq)


# INFO: the currenty riboswitch model only supports ONE hairpin loop (for each
# state) although some functionalities are already implemented to support more
# than one hairpin

class Riboswitch(object):
    def __init__(self):
        self.elements = set()
        self._elements = [[] for _ in xrange(rs_e.Type.count)]

        maxint = sys.maxint
        minint = -maxint - 1
        self.pos = [maxint, minint]
        self.pos_instance = [maxint, minint]
        self.pos_riboswitch = [maxint, minint]

    def __hash__(self):
        try:
            return self.__hash
        except AttributeError:
            # start = self.pos_instance[0]
            # target_site = (
            #     self._elements[rs_e.Type.target_site][0])
            # # consider the length of the functional site
            # self.__hash = self.pos_instance[1] - start
            # # consider each hairpin
            # for hairpin in sorted(self._elements[rs_e.Type.hairpin],
            #                       key=lambda h: (h.state, h.pos)):
            #     self.__hash ^= (
            #         hash((hairpin.pos[0] - start, hairpin.pos[1] - start)))
            #     self.__hash ^= hash(str(hairpin.struct))
            #     self.__hash <<= hairpin.state
            # # consider the pos of the target site
            # self.__hash ^= (
            #     hash((target_site.pos[0] - start, target_site.pos[1] - start)))
            # return self.__hash
            self.__hash = hash(str(self.get_raw_info()))
            return self.__hash

    def __eq__(self, other):
        # raw info contains all information that makes riboswitches
        # distinguishable
        # TODO: improve
        return self.get_raw_info() == other.get_raw_info()

    def __repr__(self):
        return self.get_raw_info()

    def reset_hash(self):
        try:
            del self.__hash
        except AttributeError:
            pass
        try:
            del self.__info
        except AttributeError:
            pass
        try:
            del self.__raw_info
        except AttributeError:
            pass

    def copy(self):
        new = copy.copy(self)

        new._elements = (
            list(list(element)
                 for element in self._elements))
        new.elements = set(self.elements)

        new.pos = list(self.pos)
        new.pos_riboswitch = list(self.pos_riboswitch)
        new.pos_instance = list(self.pos_instance)

        new.reset_hash()
        return new

    def add(self, rs_elem):
        # add the riboswitch element
        self._elements[rs_elem.type].append(rs_elem)
        self.elements.add(rs_elem)

        # set start and end position (without context they equal start and end
        # position of the riboswitch)
        if rs_elem.pos[0] < self.pos[0]:
            if rs_elem.pos[1] > self.pos[1]:
                self.pos = list(rs_elem.pos)
            else:
                self.pos[0] = rs_elem.pos[0]
        elif rs_elem.pos[1] > self.pos[1]:
            self.pos[1] = rs_elem.pos[1]

        # set positions of functional site and the hole riboswitch
        if rs_elem.type in rs_e.EXPRESSION_PLATFORM_TYPES:
            # reset the hash since it is not valid any more
            self.reset_hash()
            if rs_elem.pos[0] < self.pos_instance[0]:
                if rs_elem.pos[1] > self.pos_instance[1]:
                    self.pos_instance = list(rs_elem.pos)
                else:
                    self.pos_instance[0] = rs_elem.pos[0]
                self.pos_riboswitch[0] = rs_elem.pos[0]
            elif rs_elem.pos[1] > self.pos_instance[1]:
                self.pos_instance[1] = rs_elem.pos[1]
        elif rs_elem.type == rs_e.Type.aptamer:
            if rs_elem.pos[1] > self.pos_riboswitch[1]:
                self.pos_riboswitch[1] = rs_elem.pos[1]
            if rs_elem.pos[0] > self.pos_instance[1]:
                self.pos_instance[1] = rs_elem.pos[0]
        elif rs_elem.type == rs_e.Type.context_front:
            if rs_elem.pos[1] < self.pos_instance[0]:
                self.pos_instance[0] = rs_elem.pos[1]
                self.pos_riboswitch[0] = rs_elem.pos[1]
        elif rs_elem.type == rs_e.Type.context_back:
            if rs_elem.pos[0] > self.pos_riboswitch[1]:
                self.pos_riboswitch[1] = rs_elem.pos[0]

        # TODO: remove later
        assert(self.pos_riboswitch[0] == self.pos_instance[0])

    def replace(self, old_rs_elem, new_rs_elem):
        self.remove(old_rs_elem)
        self.add(new_rs_elem)

    def remove(self, rs_elem):
        try:
            # remove old_rs_elem if it exists
            self._elements[rs_elem.type].remove(rs_elem)
        except ValueError:
            pass

    def create_evaluation_files(self, folder):
        # NOTE: expects the container to contain only one hairpin for each state

        def write_to_file(file_, *out):
            if os.path.isfile(file_):
                raise IOError("File '%s' already exists." % file_)
            with open(file_, 'w') as f:
                for line in out:
                    f.write("%s\n" % str(line))

        if not os.path.exists(folder):
            raise OSError("Folder '%s' does not exist." % folder)

        state_ident = ["", ""]
        state_ident[rs_e.State.bound] = "b"
        state_ident[rs_e.State.unbound] = "ub"
        folders = ["", ""]
        folders[rs_e.State.bound] = (
            os.path.join(folder, state_ident[rs_e.State.bound]))
        folders[rs_e.State.unbound] = (
            os.path.join(folder, state_ident[rs_e.State.unbound]))
        folders = tuple(folders)
        os.mkdir(folders[0])
        os.mkdir(folders[1])

        ### generate the RNAf evaluation files #################################

        structs, seq = self.get_constraints_riboswitch()
        len_functional_site = (self.pos_instance[1] -
                               self.pos_instance[0])

        # hairpin
        patterns_h = ["", ""]
        for hairpin in self._elements[rs_e.Type.hairpin]:
            state = hairpin.state
            ident = "%s_%s" % (rs_e.Hairpin.ident, state_ident[state])
            patterns_h[state] = str(ident)
            struct = structs[state]
            # INFO: HACK!!!!!! fix: get_constraints_functional_site
            if state == rs_e.State.bound:
                struct_ = struct[:len_functional_site - 12]
                seq_ = seq[:len_functional_site - 12]
            else:
                struct_ = struct[:len_functional_site]
                seq_ = seq[:len_functional_site]
            write_to_file(
                os.path.join(
                    folders[state], "%s.%s" % (ident, params.BB_FILE_EXT)),
                ">%s" % ident, seq_, struct_, params.BB_TERMINATOR)

        # aptamer
        patterns_a = ["", ""]
        for aptamer in self._elements[rs_e.Type.aptamer]:
            state = aptamer.state
            ident = "%s_%s" % (rs_e.Aptamer.ident, state_ident[state])
            # INFO: HACK!!!
            if aptamer.state == rs_e.State.unbound:
                # TODO: integrate 'has to be unbound' constraint
                patterns_a[state] = ident
                struct_ = struct[len_functional_site:]
                seq_ = seq[len_functional_site:]
            else:
                patterns_a[state] = ident
                struct_ = struct[len_functional_site - 12:]
                seq_ = seq[len_functional_site - 12:]
            write_to_file(
                os.path.join(
                    folders[state], "%s.%s" % (ident, params.BB_FILE_EXT)),
                ">%s" % ident, seq_, struct_, params.BB_TERMINATOR)

        # target site
        rs = self._elements[rs_e.Type.target_site][0]
        start_pos_fs = self.pos_instance[0]
        # NOTE: +1 because of different index counting
        pos_i = rs.pos[0] - start_pos_fs + 1 + 3  # HACK: the last omits the ATG/AUG at the start
        pos_j = rs.pos[1] - start_pos_fs
        pattern_rs = "%c(%i,%i)" % (ord('%'), pos_i, pos_j)
        header = ">%s" % rs_e.TargetSite.ident
        write_to_file(
            os.path.join(folders[0], "%s.%s" % (rs_e.TargetSite.ident,
                                                params.BB_FILE_EXT)),
            header, rs.seq, params.BB_TERMINATOR)
        write_to_file(
            os.path.join(folders[1], "%s.%s" % (rs_e.TargetSite.ident,
                                                params.BB_FILE_EXT)),
            header, rs.seq, params.BB_TERMINATOR)

        # INFO: HACK! --> only one accessibility constraint expected which is
        #                 used for the unbound case
        # accessibility constraints
        try:
            ac = self._elements[rs_e.Type.access_constraint][0]
        except IndexError:
            pattern_ac = ""
        else:
            aptamers = self._elements[rs_e.Type.aptamer]
            aptamer_ub = (aptamers[0]
                          if aptamers[0].state == rs_e.State.unbound
                          else aptamers[1])
            # NOTE: +1 because of different index counting
            pos_i = ac.pos[0] - aptamer_ub.pos[0] + 1
            pos_j = ac.pos[1] - aptamer_ub.pos[0]
            pattern_ac = "%c(%i,%i)" % (ord('%'), pos_i, pos_j)
            header = ">%s" % rs_e.AccessConstraint.ident
            write_to_file(
                os.path.join(
                    folders[rs_e.State.unbound],
                    "%s.%s" % (rs_e.AccessConstraint.ident,
                               params.BB_FILE_EXT)),
                header, ac.seq)

        # context front (if existent)
        try:
            cf = self._elements[rs_e.Type.context_front][0]
        except IndexError:
            pattern_cf = ""
        else:
            header = ">%s" % rs_e.ContextFront.ident
            pattern_cf = "%s%s" % (rs_e.ContextFront.ident,
                                   unichr(167).encode("utf-8"))
            write_to_file(
                os.path.join(folders[0], "%s.%s" % (rs_e.ContextFront.ident,
                                                    params.BB_FILE_EXT)),
                header, cf.seq, cf.struct, params.BB_TERMINATOR)
            write_to_file(
                os.path.join(folders[1], "%s.%s" % (rs_e.ContextFront.ident,
                                                    params.BB_FILE_EXT)),
                header, cf.seq, cf.struct, params.BB_TERMINATOR)

        # context back (if existent)
        try:
            cb = self._elements[rs_e.Type.context_back][0]
        except IndexError:
            pattern_cb = ""
        else:
            header = ">%s" % rs_e.ContextBack.ident
            pattern_cb = "%s%s" % (unichr(167).encode("utf-8"),
                                   rs_e.ContextBack.ident)
            write_to_file(
                os.path.join(folders[0], "%s.%s" % (rs_e.ContextBack.ident,
                                                    params.BB_FILE_EXT)),
                header, cb.seq, cb.struct, params.BB_TERMINATOR)
            write_to_file(
                os.path.join(folders[1], "%s.%s" % (rs_e.ContextBack.ident,
                                                    params.BB_FILE_EXT)),
                header, cb.seq, cb.struct, params.BB_TERMINATOR)

        # pattern file
        state = rs_e.State.bound
        write_to_file(
            os.path.join(folders[state], params.PATTERN_FILE_NAME),
            ">pattern_%s" % state_ident[state],
            PATTERNS[state] % (
                pattern_cf, pattern_rs, patterns_h[state], patterns_a[state],
                pattern_cb))
        state = rs_e.State.unbound
        write_to_file(
            os.path.join(folders[state], params.PATTERN_FILE_NAME),
            ">pattern_%s" % state_ident[state],
            PATTERNS[state] % (
                pattern_cf, pattern_rs, patterns_h[state], pattern_ac,
                patterns_a[state], pattern_cb))

        return folders

    def get_accessibility_positions(self):
        start = self.pos[0]
        rs = list(self._elements[rs_e.Type.target_site][0].pos)
        rs[0] += 3  # HACK!!!
        try:
            ac = self._elements[rs_e.Type.access_constraint][0].pos
        except IndexError:
            raise ValueError("Access constraint should be available!")
        a_b = sorted(self._elements[rs_e.Type.aptamer])[rs_e.State.bound]
        h_0, h_1 = self._elements[rs_e.Type.hairpin]
        l = max(a_b.pos[1] - a_b.pos[0],
                h_0.pos[1] - h_0.pos[0],
                h_1.pos[1] - h_1.pos[0]) + 2
        return (rs[0] - start, ac[0] - start), (rs[1] - rs[0], ac[1] - ac[0]), l

    def get_constraints(self):
        riboswitch_elements = (
            itertools.chain.from_iterable(self._elements))
        return get_constraints(self.pos, riboswitch_elements)

    def get_constraints_functional_site(self):
        # since the functional site can overlap into the aptamer site in this
        # case a part of the aptamer sequence has to be considered in the
        # functional site
        # 'how much' sequence of the aptamer site is going to be considered is
        # decided in get_constraints(..)
        a_0, a_1 = self._elements[rs_e.Type.aptamer]
        if a_0.pos[0] < a_1.pos[0]:
            aptamer_seq = rs_e.SequenceConstraint(a_0.pos, a_0.seq)
        elif a_0.pos[0] > a_1.pos[0]:
            aptamer_seq = rs_e.SequenceConstraint(a_1.pos, a_1.seq)
        else:
            aptamer_seq = None

        if aptamer_seq:
            func_site_elements = (
                itertools.chain(
                    self._elements[rs_e.Type.hairpin],
                    self._elements[rs_e.Type.target_site],
                    self._elements[rs_e.Type.seq_constraint],
                    (aptamer_seq,)))
        else:
            func_site_elements = (
                itertools.chain(
                    self._elements[rs_e.Type.hairpin],
                    self._elements[rs_e.Type.target_site],
                    self._elements[rs_e.Type.seq_constraint]))

        return get_constraints(self.pos_instance, func_site_elements)

    def get_constraints_riboswitch(self):
        riboswitch_elements = (
            itertools.chain(
                self._elements[rs_e.Type.hairpin],
                self._elements[rs_e.Type.target_site],
                self._elements[rs_e.Type.seq_constraint],
                self._elements[rs_e.Type.aptamer]))
        return get_constraints(self.pos_riboswitch, riboswitch_elements)

    def get_raw_info(self):
        try:
            return self.__raw_info
        except AttributeError:
            info = []

            info.append("pos_fs=%s" % str(tuple(self.pos_instance)))

            try:
                cf = self._elements[rs_e.Type.context_front][0]
                info.append(repr(cf))
            except IndexError:
                pass

            target_site = self._elements[rs_e.Type.target_site][0]
            info.append(repr(target_site))

            for hairpin in sorted(self._elements[rs_e.Type.hairpin],
                                  key=lambda h: (h.state, h.pos)):
                info.append(repr(hairpin))

            for aptamer in sorted(self._elements[rs_e.Type.aptamer],
                                  key=lambda a: (a.state, a.pos)):
                info.append(repr(aptamer))

            for seq_constraint in sorted(self._elements[rs_e.Type.seq_constraint],
                                         key=lambda s: s.pos):
                info.append(repr(seq_constraint))

            try:
                acc = self._elements[rs_e.Type.access_constraint][0]
                info.append(repr(acc))
            except IndexError:
                pass

            try:
                cb = self._elements[rs_e.Type.context_back][0]
                info.append(repr(cb))
            except IndexError:
                pass

            self.__raw_info = ';'.join(info)
            return self.__raw_info

    def get_full_info(self):
        info = []
        # start = self.pos_instance[0]
        start = self.pos[0]

        (struct_0, struct_1), seq = self.get_constraints_riboswitch()
        pos = self.pos_riboswitch[0] - start, self.pos_riboswitch[1] - start
        info.append("Riboswitch %s %s %s %s" %
                    (''.join(struct_0), ''.join(struct_1), ''.join(seq), pos))

        (struct_0, struct_1), seq = self.get_constraints()
        pos = self.pos[0] - start, self.pos[1] - start
        info.append("All constraints %s %s %s %s" %
                    (''.join(struct_0), ''.join(struct_1), ''.join(seq), pos))

        info.append("Main riboswitch_elements:")
        for hairpin in self._elements[rs_e.Type.hairpin]:
            pos = hairpin.pos[0] - start, hairpin.pos[1] - start
            info.append("   Hairpin %s %s %s" %
                        (hairpin.struct, pos,
                         rs_e.State.get_str(hairpin.state)))
        for aptamer in self._elements[rs_e.Type.aptamer]:
            pos = aptamer.pos[0] - start, aptamer.pos[1] - start
            info.append("   Aptamer %s %s %s %s" %
                        (aptamer.struct, aptamer.seq, pos,
                         rs_e.State.get_str(aptamer.state)))
        target_site = (
            self._elements[rs_e.Type.target_site][0])
        pos = target_site.pos[0] - start, target_site.pos[1] - start
        info.append("   Restriction site %s %s" % (target_site.seq, pos))
        for seq_constraint in self._elements[rs_e.Type.seq_constraint]:
            pos = seq_constraint.pos[0] - start, seq_constraint.pos[1] - start
            info.append("   Sequence constraint %s %s" %
                        (seq_constraint.seq, pos))

        if (len(self._elements[rs_e.Type.context_front]) or
            len(self._elements[rs_e.Type.context_back])):
            info.append("Additional riboswitch_elements:")
            try:
                context_front = (
                    self._elements[rs_e.Type.context_front][0])
                pos = context_front.pos[0] - start, context_front.pos[1] - start
                info.append("   Kontext front %s %s %s" %
                            (context_front.struct, context_front.seq, pos))
            except IndexError:
                pass
            try:
                context_back = (
                    self._elements[rs_e.Type.context_back][0])
                pos = context_back.pos[0] - start, context_back.pos[1] - start
                info.append("   Kontext back %s %s %s" %
                            (context_back.struct, context_back.seq, pos))
            except IndexError:
                pass

        return '\n'.join(info)

    def validate(self):
        def validate_overlapping_part(rs_elem_a, rs_elem_b):
            # check if the riboswitch_elements are not overlapping
            if rs_elem_a.pos[1] < rs_elem_b.pos[0]:
                return
            if rs_elem_a.pos[0] > rs_elem_b.pos[1]:
                return
            if not rs_elem_a.seq or not rs_elem_b.seq:
                return

            start = (rs_elem_a.pos[0]
                     if rs_elem_a.pos[0] > rs_elem_b.pos[0]
                     else rs_elem_b.pos[0])
            end = (rs_elem_a.pos[1]
                   if rs_elem_a.pos[1] < rs_elem_b.pos[1]
                   else rs_elem_b.pos[1])

            # check if they are compatible
            compatible = True
            for pos in xrange(start, end):
                # compare the sequence constraints base by base
                seq_elem_a = rs_elem_a.seq[pos - rs_elem_a.pos[0]]
                seq_elem_b = rs_elem_b.seq[pos - rs_elem_b.pos[0]]
                compatible = False
                for b_id in xrange(rna.BaseId.count):
                    iupac_id_i = getattr(rna.IUPAC_Id, seq_elem_a)
                    iupac_id_j = getattr(rna.IUPAC_Id, seq_elem_b)
                    if (rna.base_valid_IUPAC(iupac_id_i, b_id) and
                        rna.base_valid_IUPAC(iupac_id_j, b_id)):
                        compatible = True
                        break
                if not compatible:
                    raise ValueError(
                        " '%s' (%i, %i) and rs_elem '%s' (%i, %i) are not "
                        "compatible." %
                        (rs_elem_a, rs_elem_a.pos[0], rs_elem_a.pos[1],
                         rs_elem_b, rs_elem_b.pos[0], rs_elem_b.pos[1]))

        ### validate length of the functional site #############################
        size = self.pos_instance[1] - self.pos_instance[0]
        if size > EXPRESSION_PLATFORM_MAX_LEN:
            raise ValueError(
                "Functional site exceeds maximal length of the functional "
                "site.")

        ### validate the context front #########################################
        # check if builder contains more than one KontextFront
        if len(self._elements[rs_e.Type.context_front]) > 1:
            raise ValueError("Riboswitch contains more than one front kontext.")
        # tests concerning the functional site are done later

        ### validate the target site ######################################
        target_site = (
            self._elements[rs_e.Type.target_site])
        # check if builder contains another TargetSite
        if len(target_site) != 1:
            raise ValueError(
                "Riboswitch contains no or too many TargetSite "
                "elements.")
        target_site = target_site[0]
        # check if the TargetSite intersects with the context front
        if len(self._elements[rs_e.Type.context_front]):
            context_front = (
                self._elements[rs_e.Type.context_front][0])
            if target_site.pos[0] < context_front.pos[1]:
                raise ValueError(
                    "Restriction site intersects with front kontext.")
        # tests concerning the functional site are done later

        ### validate sequence constraints ######################################
        for seq_constraint in self._elements[rs_e.Type.seq_constraint]:
            # check if the seq intersects with the context front
            if len(self._elements[rs_e.Type.context_front]):
                context_front = (
                    self._elements[rs_e.Type.context_front][0])
                if seq_constraint.pos[0] < context_front.pos[1]:
                    raise ValueError(
                        "Riboswitch element sequence intersects with "
                        "front kontext.")

        ### validate the hairpins ##############################################
        hairpins = self._elements[rs_e.Type.hairpin]
        # check if at least one Hairpin exists
        if not len(hairpins):
            raise ValueError("Riboswitch contains no Hairpin.")
        for i, hairpin in enumerate(hairpins):
            # check if the Hairpin intersects with the context front
            if len(self._elements[rs_e.Type.context_front]):
                context_front = (
                    self._elements[rs_e.Type.context_front][0])
                if hairpin.pos[0] < context_front.pos[1]:
                    raise ValueError(
                        "Hairpin intersects with front kontext.")
            # validate overlapping parts with other Hairpin elements
            add = True
            for hairpin_2 in hairpins[i + 1:]:
                if hairpin.state == hairpin_2.state:
                    # check if the Hairpins of the same state are
                    # overlapping
                    start, end = hairpin.pos
                    start_b, end_b = hairpin_2.pos
                    if end <= start_b or end_b <= start:
                        add = True
                    elif start >= end_b or start_b >= end:
                        add = True
                    else:
                        add = False
                        break
                else:
                    # compare the sequences at the overlapping positions for
                    # hairpin elements of different states
                    validate_overlapping_part(hairpin, hairpin_2)
            if not add:
                raise ValueError(
                    "Riboswitch contains Hairpin elements that do not"
                    " match.")
            # validate overlapping parts with the target site
            target_site = (
                self._elements[rs_e.Type.target_site][0])
            validate_overlapping_part(hairpin, target_site)
            # # the target site has to lie within the bound hairpin element
            # if hairpin.state == rs_e.State.bound:
            #     if target_site.pos[1] <= hairpin.pos[0]:
            #         raise ValueError(
            #             "The target site has to lie within the bound "
            #             "hairpin element.")
            # # the target site cannot overlap with the unbound hairpin
            # # element
            # else:
            #     if target_site.pos[1] > hairpin.pos[1]:
            #         raise ValueError(
            #             "The target site cannot overlap with the unbound "
            #             "hairpin element.")
            # validate overlapping parts with sequence constraint elements
            for seq_constraint in self._elements[rs_e.Type.seq_constraint]:
                validate_overlapping_part(hairpin, seq_constraint)

        ### validate the aptamer site ##########################################
        aptamers = self._elements[rs_e.Type.aptamer]
        # check if aptamer site intersects with the functional site
        for aptamer in aptamers:
            if any(hairpin
                   for hairpin
                   in self._elements[rs_e.Type.hairpin]
                   if aptamer.state == hairpin.state and
                      aptamer.pos[0] < hairpin.pos[1]):
                raise ValueError(
                    "Aptamer intersects with functional site.")
            target_site = (
                self._elements[rs_e.Type.target_site][0])
            if aptamer.pos[0] < target_site.pos[1]:
                raise ValueError(
                    "Aptamer intersects with functional site.")
            if any(seq_constraint
                   for seq_constraint
                   in self._elements[rs_e.Type.seq_constraint]
                   if aptamer.pos[0] < seq_constraint.pos[1]):
                raise ValueError(
                    "Aptamer intersects with functional site.")
        # check if builder contains an aptamer element for both states
        if len(aptamers) != 2:
            raise ValueError(
                "Aptamer is missing or the aptamer is defined for only one"
                "bounding state.")
        if aptamers[0].state == aptamers[1].state:
            raise ValueError(
                "Riboswitch contains two aptamers for the same bounding state.")
        # validate the overlapping parts of both aptamer riboswitch_elements
        validate_overlapping_part(aptamers[0], aptamers[1])

        ### validate the context back ##########################################
        context_back = self._elements[rs_e.Type.context_back]
        if len(context_back):
            # check if builder contains more than one KontextBack
            if len(context_back) > 1:
                raise ValueError(
                    "Riboswitch contains more than one back kontext")
            context_back = context_back[0]
            # check if the KontextBack intersects with the aptamer site
            if any(aptamer
                   for aptamer
                   in self._elements[rs_e.Type.aptamer]
                   if context_back.pos[0] < aptamer.pos[1]):
                raise ValueError("Back kontext intersects with aptamer site.")
