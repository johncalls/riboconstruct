#!/usr/bin/env python

import unittest

from tests import test_dependency_graph
from tests import test_helper
from tests import test_inverse_fold
from tests import test_local_refinement
from tests import test_riboswitch
from tests import test_rna
from tests import test_structure
from tests import test_two_target_inverse_fold
from tests import test_two_target_subsolution


if __name__ == "__main__":
    suite_test_dependency_graph = (
        unittest.TestLoader().loadTestsFromModule(test_dependency_graph))
    suite_test_helper = unittest.TestLoader().loadTestsFromModule(test_helper)
    suite_test_inverse_fold = (
        unittest.TestLoader().loadTestsFromModule(test_inverse_fold))
    suite_test_local_refinement = (
        unittest.TestLoader().loadTestsFromModule(test_local_refinement))
    suite_test_riboswitch = (
        unittest.TestLoader().loadTestsFromModule(test_riboswitch))
    suite_test_rna = unittest.TestLoader().loadTestsFromModule(test_rna)
    suite_test_structure = (
        unittest.TestLoader().loadTestsFromModule(test_structure))
    suite_test_two_target_inverse_fold = (
        unittest.TestLoader().loadTestsFromModule(test_two_target_inverse_fold))
    suite_test_two_target_subsolution = (
        unittest.TestLoader().loadTestsFromModule(test_two_target_subsolution))

    alltests = unittest.TestSuite(
        [suite_test_dependency_graph,
         suite_test_inverse_fold,
         suite_test_local_refinement,
         suite_test_helper,
         suite_test_riboswitch,
         suite_test_rna,
         suite_test_structure,
         suite_test_two_target_inverse_fold,
         suite_test_two_target_subsolution])

    unittest.TextTestRunner().run(alltests)
