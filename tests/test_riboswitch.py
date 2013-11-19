#!/usr/bin/env python

import unittest

from riboconstruct import rna

from riboconstruct import riboswitch as rs
from riboconstruct.riboswitch import element as rs_e


class TestRiboswitch(unittest.TestCase):
    def setUp(self):
        rna.STRUCT_ELEM_UNSPEC = '*'

        self.riboswitch = rs.Riboswitch()
        self.h_ub_1 = (
            rs_e.Hairpin(rs_e.State.unbound, (6, 11), '(...)'))
        self.h_b = (
            rs_e.Hairpin(rs_e.State.bound, (9, 14), '(...)'))
        self.h_ub_2 = (
            rs_e.Hairpin(rs_e.State.unbound, (11, 16), '(...)'))
        self.cf = rs_e.ContextFront((1, 6), "AUUUU")
        self.cb = rs_e.ContextBack((25, 32), 'GCCAUAA')
        self.rs = rs_e.TargetSite((9, 13), 'AGUU')
        self.a_ub = (
            rs_e.Aptamer(
                rs_e.State.unbound, (17, 25), '.(...)..', 'CGCCAUAA'))
        self.a_b = (
            rs_e.Aptamer(
                rs_e.State.bound, (14, 22), '.(...)..', 'UGUCGCCA'))

    def tearDown(self):
        rna.STRUCT_ELEM_UNSPEC = '.'

    def test_fsm_initialization(self):
        rs_e.Hairpin(rs_e.State.unbound, (0, 5), '(...)')
        rs_e.Hairpin(rs_e.State.unbound, (0, 7), '(...)..')

    def test_add_fsm(self):
        self.riboswitch.add(self.h_ub_1)

    def test_constraint_building_1(self):
        self.riboswitch.add(self.h_ub_1)
        self.riboswitch.add(self.h_b)
        self.riboswitch.add(self.h_ub_2)
        self.riboswitch.add(self.rs)
        self.riboswitch.add(self.a_ub)
        self.riboswitch.add(self.a_b)

        self.assertEqual([6, 25], self.riboswitch.pos_instance)
        self.assertEqual([6, 25], self.riboswitch.pos_riboswitch)

        (struct_ub, struct_b), seq = self.riboswitch.get_constraints()
        self.assertEqual('(...)(...)*.(...)..', ''.join(struct_ub))
        self.assertEqual('***(...).(...)..***', ''.join(struct_b))
        self.assertEqual('NNNAGUUNUGUCGCCAUAA', ''.join(seq))

        (struct_ub, struct_b), seq = self.riboswitch.get_constraints_riboswitch()
        self.assertEqual('(...)(...)*.(...)..', ''.join(struct_ub))
        self.assertEqual('***(...).(...)..***', ''.join(struct_b))
        self.assertEqual('NNNAGUUNUGUCGCCAUAA', ''.join(seq))

    def test_constraint_building_2(self):
        self.riboswitch.add(self.h_ub_1)
        self.riboswitch.add(self.h_b)
        self.riboswitch.add(self.h_ub_2)
        self.riboswitch.add(self.cf)
        self.riboswitch.add(self.cb)
        self.riboswitch.add(self.rs)
        self.riboswitch.add(self.a_ub)
        self.riboswitch.add(self.a_b)

        self.assertEqual([1, 32], self.riboswitch.pos_instance)
        self.assertEqual([6, 25], self.riboswitch.pos_riboswitch)

        (struct_ub, struct_b), seq = self.riboswitch.get_constraints()
        self.assertEqual('*****(...)(...)*.(...)..*******', ''.join(struct_ub))
        self.assertEqual('********(...).(...)..**********', ''.join(struct_b))
        self.assertEqual('AUUUUNNNAGUUNUGUCGCCAUAAGCCAUAA', ''.join(seq))

        (struct_ub, struct_b), seq = (
            self.riboswitch.get_constraints_riboswitch())
        self.assertEqual('(...)(...)*.(...)..', ''.join(struct_ub))
        self.assertEqual('***(...).(...)..***', ''.join(struct_b))
        self.assertEqual('NNNAGUUNUGUCGCCAUAA', ''.join(seq))

    def test_serialization(self):
        self.riboswitch.add(self.h_ub_1)
        self.riboswitch.add(self.h_b)
        self.riboswitch.add(self.h_ub_2)
        self.riboswitch.add(self.rs)
        self.riboswitch.add(self.a_ub)
        self.riboswitch.add(self.a_b)

        rs_str = repr(self.riboswitch)
        riboswitch_from_str = rs.get_riboswitch_from_str(rs_str)
        self.assertEqual(self.riboswitch, riboswitch_from_str)


class TestRiboswitchParser(unittest.TestCase):
    def setUp(self):
        self.configs_example = (
            '>  aPtamer  |17,25|unbound',
            ' .(...)..  ',
            'CGCCAUAA  ',
            ' >hairpin|6, 11 |  UNbound',
            '(...)',
            '>hairpin| 9  ,14|bOund  |',
            '(...)',
            '>hairpin|11,16|unbound',
            ' (...)  ',
            '>context_frONt|1,6|',
            'AUUUU',
            ' > context_back| 25 , 32 | ',
            'GCCAUAA',
            '>target_sitE|9,13|',
            'AGUU',
            '>aptamer|14,22|bound',
            '.(...)..',
            'UGUCGCCA',
            '>access_constraint|17,18',
            'c')

    def test_riboswitch_parser(self):
        for element in rs.iter_riboswitch(self.configs_example):
            if element.type == rs_e.Type.aptamer:
                if element.state == rs_e.State.unbound:
                    self.assertTupleEqual((17, 25), element.pos)
                    self.assertEqual('.(...)..', str(element.struct))
                    self.assertEqual('CGCCAUAA', str(element.seq))
                elif element.state == rs_e.State.bound:
                    self.assertTupleEqual((14, 22), element.pos)
                    self.assertEqual('.(...)..', str(element.struct))
                    self.assertEqual('UGUCGCCA', str(element.seq))
                else:
                    raise ValueError('Wrong state.')
            elif element.type == rs_e.Type.hairpin:
                if element.state == rs_e.State.unbound:
                    pass
                elif element.state == rs_e.State.bound:
                    self.assertTupleEqual((9, 14), element.pos)
                    self.assertEqual('(...)', str(element.struct))
                else:
                    raise ValueError('Wrong state.')
            elif element.type == rs_e.Type.target_site:
                self.assertTupleEqual((9, 13), element.pos)
                self.assertEqual('AGUU', str(element.seq))
            elif element.type == rs_e.Type.context_front:
                self.assertTupleEqual((1, 6), element.pos)
                self.assertEqual('AUUUU', str(element.seq))
            elif element.type == rs_e.Type.context_back:
                self.assertTupleEqual((25, 32), element.pos)
                self.assertEqual('GCCAUAA', str(element.seq))
            elif element.type == rs_e.Type.access_constraint:
                self.assertTupleEqual((17, 18), element.pos)
            else:
                raise ValueError('Wrong mode.')

    def test_aptamer_element_parser_fails(self):
        # missing element
        config = ('>', ' .(...)..  ', 'CGCCAUAA  ')
        with self.assertRaises(ValueError):
            for element in rs.iter_riboswitch(config):
                pass
        # wrong element
        config = ('>multiloop|1,10', ' .(...)..  ', 'CGCCAUAA  ')
        with self.assertRaises(AttributeError):
            for element in rs.iter_riboswitch(config):
                pass
        # each element line (starting with '>') should contain position info at
        # second
        config = ('>aptamer ', ' .(...)..  ', 'CGCCAUAA  ')
        with self.assertRaises(ValueError):
            for element in rs.iter_riboswitch(config):
                pass
        # the position should be simply two numbers separated by a ','
        config = ('> aptamer|(1,5)', ' .(...)..  ', 'CGCCAUAA  ')
        with self.assertRaises(ValueError):
            for element in rs.iter_riboswitch(config):
                pass
        # missing state info
        config = (' >aptamer|1,5|', ' .(...)..  ', 'CGCCAUAA  ')
        with self.assertRaises(ValueError):
            for element in rs.iter_riboswitch(config):
                pass
        # wrong sequence --> L is not a valid IUPAC base
        config = (' >aptamer|1,5|unbound', ' .(...)..  ', 'CGCLAUAA  ')
        with self.assertRaises(TypeError):
            for element in rs.iter_riboswitch(config):
                pass
        # wrong struct
        config = (' >aptamer|1,5|unbound', ' ).(...)..  ', 'CGCCAUAA  ')
        with self.assertRaises(ValueError):
            for element in rs.iter_riboswitch(config):
                pass


if __name__ == "__main__":
    unittest.main()
