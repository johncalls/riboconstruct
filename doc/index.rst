.. _multiprocessing system: http://en.wikipedia.org/wiki/Multiprocessing
.. _riboswitches: http://en.wikipedia.org/wiki/Riboswitch

===========================================================
:mod:`riboconstruct` --- Generate and evaluate riboswitches
===========================================================

This module provides methods and techniques to generate `riboswitches`_ for a given model and evaluate these with respect to their folding behaviour.

:mod:`riboconstruct` is mainly based on two subpackages:

* :mod:`riboconstruct.riboswitch` provides methods to traverse and evaluate the space of riboswitch instances defined by the underlying model. A riboswitch instance fixes the different parts of the riboswitch model into a concrete setting by defining the size, structure or position of its elements.

  The actual evaluation of a riboswitch instance is done using :mod:`riboconstruct.rna_f`.

* Once a riboswitch instance is fixed while iterating the riboswitch space, only the two structures and *some* of the bases are known. Therefore, a RNA sequence that is likely to fold into the structures has to be identified. The necessary methods are provided in :mod:`riboconstruct.inverse_folding`.

To generate and evaluate riboswitches the following subpackage can be used:

* :mod:`riboconstruct.eval_mp` provides methods to generate and evaluate a riboswitch model within given constraints. The methods are based on Python's :mod:`multiprocessing` and can be used to do the calculations on a `multiprocessing system`_ in parallel.

In addition, :mod:`riboconstruct` offers helper classes and functions e.g. to represent RNA structures or sequences:

* :mod:`riboconstruct.rna`
* :mod:`riboconstruct.helper`


Examples
========

How to use instance spaces and instance iterators of the :mod:`~riboconstruct.riboswitch.generation` module.

.. literalinclude:: ../use_cases/splice_site_iter.py


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

