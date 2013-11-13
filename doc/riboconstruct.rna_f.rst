.. module:: riboconstruct.rna_f

riboconstruct.rna_f package
===========================

The module provides methods to evaluate the switching behaviour of riboswitches by wrapping **RNAf** which is by Robert Kleinkauf.

.. note::
  
  In order to use the python wrapper, the underlying C wrapper has to be compiled first.
  
  1. Go to the ``./riboconstruct/rna_f`` subdirectory

  2. Customize the following paths in the *Makefile*:
  
     * VIENNA: path/to/the/vienna/rna/package
     
     * RNA_F: path/to/rna_d
     
     * INSTALL_DIR: path/to/this/subdirectory
     
  3. Execute ``make`` in the same subdirectory


.. autoclass:: riboconstruct.rna_f.RNAf
    :members:
    :undoc-members:
    :show-inheritance:
