.. module:: riboconstruct.inverse_folding

riboconstruct.inverse_folding package
=====================================

Inverse folding is the task of identifying a sequence that is likely to fold into a given RNA target structure. :mod:`riboconstruct.inverse_folding` implements the inverse folding for the one and two target case.

The package is organized as follows:

* :class:`riboconstruct.inverse_folding.structure` provides methods to handle RNA structures and their substructures.

* :mod:`riboconstruct.inverse_folding.fillings` provides methods to calculate the free energy fractions between two base pairs according to the nearest neighbour energy model and to assign bases such that a minimal free energy fraction is provided.

* :mod:`riboconstruct.inverse_folding.subsolution` defines methods to traverse the target structure along its substructures and to assign bases to its positions. 

* :mod:`riboconstruct.inverse_folding.two_target_subsolution` extends the methods from the single target case to calculate RNA sequences likely to fold into two target structures.

  * :mod:`riboconstruct.inverse_folding.dependency_graph` is used to mark energy dependent base positions.
  
  * :mod:`riboconstruct.inverse_folding.two_target_local_refinement` is used to improve the identified RNA sequences by local search.

Calculations:

* :func:`riboconstruct.inverse_folding.inverse_fold` calculates the RNA sequence for the one target case.

* The methods in :mod:`riboconstruct.inverse_folding.two_target_inverse_fold` are used to calculate RNA sequences for the two target case.


.. automodule:: riboconstruct.inverse_folding
    :members:
    :undoc-members:
    :show-inheritance:


Submodules
----------

inverse_fold
^^^^^^^^^^^^

.. automodule:: riboconstruct.inverse_folding.inverse_fold
    :members:
    :undoc-members:
    :show-inheritance:

two_target_inverse_fold
^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: riboconstruct.inverse_folding.two_target_inverse_fold
    :members:
    :undoc-members:
    :show-inheritance:

structure
^^^^^^^^^

.. automodule:: riboconstruct.inverse_folding.structure
    :members:
    :undoc-members:
    :show-inheritance:

subsolution
^^^^^^^^^^^

.. automodule:: riboconstruct.inverse_folding.subsolution
    :members:
    :undoc-members:
    :show-inheritance:

dependency_graph
^^^^^^^^^^^^^^^^

.. automodule:: riboconstruct.inverse_folding.dependency_graph
    :members:
    :undoc-members:
    :show-inheritance:

two_target_local_refinement
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: riboconstruct.inverse_folding.two_target_local_refinement
    :members:
    :undoc-members:
    :show-inheritance:
