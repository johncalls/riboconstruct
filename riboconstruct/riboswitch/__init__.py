from . import riboswitch as rs
from . import element as rs_e

from .. import rna


def iter_riboswitch(config):
    """
    Yields the parsed riboswitch elements that are encoded in *config*
    in a FASTA-like format.

    Example: ::

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
        .............................................
        AAGCUAUACCAAGCAUACAAUCAACUCCAAGCUAGAUCUCUUAAG
        >Context_back|103,148
        .............................................
        AUCUAGCGCUGGUACCAUCCCAUUUAACUGUAAGAAGAAUUGCAC
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
            elif re_type == rs_e.Type.restriction_site:
                seq = rna.IUPACSequence(line.upper())
                if pos[1] - pos[0] != len(seq):
                    raise ValueError("Restriction site element has wrong size.")
                yield rs_e.RestrictionSite(pos, seq)
            elif re_type == rs_e.Type.context_front:
                struct = rna.Structure(line)
                seq = rna.IUPACSequence(next().upper())
                if pos[1] - pos[0] != len(struct) != len(seq):
                    raise ValueError("Context front element has wrong size.")
                yield rs_e.ContextFront(pos, struct, seq)
            elif re_type == rs_e.Type.context_back:
                struct = rna.Structure(line)
                seq = rna.IUPACSequence(next().upper())
                if pos[1] - pos[0] != len(struct) != len(seq):
                    raise ValueError("Context back element has wrong size.")
                yield rs_e.ContextBack(pos, struct, seq)
            elif re_type == rs_e.Type.seq:
                seq = rna.IUPACSequence(line.upper())
                if pos[1] - pos[0] != len(seq):
                    raise ValueError("Sequence constraint element has wrong "
                                     "size.")
                yield rs_e.Sequence(pos, seq)
            elif re_type == rs_e.Type.access_constraint:
                seq = rna.IUPACSequence(line.upper())
                if pos[1] - pos[0] != len(seq):
                    raise ValueError(
                        "Access constraint element has wrong size.")
                yield rs_e.AccessConstraint(pos, seq)
            else:
                raise TypeError("Riboswitch config file contains invalid "
                                "riboswitch element.")
            re_type = None
        else:
            raise ValueError("Config has wrong format.")


def get_riboswitch_from_config_file(config_file):
    """
    Returns a riboswitch instance by parsing the elements from a
    FASTA-like format using :func:`riboconstruct.helper.iter_riboswitch`.
    """
    riboswitch = rs.Riboswitch()
    # read initial settings from config
    with open(config_file) as config:
        map(riboswitch.add, iter_riboswitch(config))
    # check if the rs elements 'fit' together
    riboswitch.validate()
    return riboswitch
