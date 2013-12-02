.. module:: riboconstruct.riboswitch

riboconstruct.riboswitch package
================================

:mod:`riboconstruct.riboswitch` provides functionalities to represent riboswitches and to iterate all instances of a certain riboswitch within defined boundaries.

* :class:`~riboconstruct.riboswitch.Riboswitch` is used to represent riboswitches. It is basically a collection of riboswitch  :class:`~riboconstruct.riboswitch.element.Element`\ s which define the actual riboswitch by their size, structure and sequence.

  These elements can simpy be added and removed from the riboswitch. Furthermore, methods like :func:`~riboconstruct.riboswitch.Riboswitch.get_constraints` can be used to get the actual structures for the bound and unbound state as well as the sequence of a specific riboswitch instance.

* :class:`~riboconstruct.riboswitch.element` and its submodules provide functionalities to handle the elements of a riboswitch like hairpins, the aptamer or target sites within a riboswitch. These elements define the actual riboswitch.

* :mod:`~riboconstruct.riboswitch.generation` provides functionalities to iterate all possible riboswitch instances within defined boundaries and to generate the information files that needed by ``RNAf`` in order to evaluate the riboswitch instances for their folding behaviour.


Module
~~~~~~

.. automodule:: riboconstruct.riboswitch
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource
    :special-members: __eq__, __repr__, __hash__


Submodules
~~~~~~~~~~

element
^^^^^^^

.. automodule:: riboconstruct.riboswitch.element
   :members:
   :undoc-members:
   :show-inheritance:
   :member-order: bysource
   :private-members:

generation
^^^^^^^^^^

.. automodule:: riboconstruct.riboswitch.generation
   :members:
   :undoc-members:
   :show-inheritance:
   :inherited-members:
   :member-order: bysource
