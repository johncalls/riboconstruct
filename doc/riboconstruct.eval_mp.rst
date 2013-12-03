.. module:: riboconstruct.eval_mp

riboconstruct.eval_mp package
=============================

:mod:`riboconstruct.eval_mp` provides functionalities to generate and evaluate riboswitches on a multiprocessing system in parallel using Python's :mod:`multiprocessing`.

The idea is that each submodule within the package is responsible for one step in the generation and evaluation of riboswitches. The submodules communicate with each other either by the file system (i.e. reading/writing information from/to the local file system) or directly via :class:`multiprocessing.Queue`\ s.

The package is organized as follows:

* :class:`~riboconstruct.eval_mp.creator` generates the evaluation directory structure and either writes the riboswitches which have to be evaluated to the evaluation directory or puts them into the output queue.

* :class:`~riboconstruct.eval_mp.generator`

* :class:`~riboconstruct.eval_mp.improver`

* :class:`~riboconstruct.eval_mp.evaluator`

* :class:`~riboconstruct.eval_mp.selector`


Submodules
----------

creator
^^^^^^^

.. automodule:: riboconstruct.eval_mp.creator
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

evaluator
^^^^^^^^^

.. automodule:: riboconstruct.eval_mp.evaluator
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

generator
^^^^^^^^^

.. automodule:: riboconstruct.eval_mp.generator
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

improver
^^^^^^^^

.. automodule:: riboconstruct.eval_mp.improver
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

selector
^^^^^^^^

.. automodule:: riboconstruct.eval_mp.selector
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. %settings
.. ^^^^^^^^

.. automodule:: riboconstruct.eval_mp.settings
..     :members:
..     :undoc-members:
..     :show-inheritance:
..     :member-order: bysource
