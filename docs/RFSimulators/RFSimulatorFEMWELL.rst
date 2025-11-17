RFSimulatorFEMWELL
==================

.. _FEMWELL: https://femwell.readthedocs.io
.. _skfem: https://scikit-fem.readthedocs.io/en/latest/
.. _MESHWELL: https://simbilod.github.io/meshwell/intro.html


The `RFSimulatorFEMWELL` is a FEM mode solver built around FEMWELL_ and skfem_. With this tool we can calculate the RF mode supported by our structure. We can retrieve the propagation constant, and calculate the characteristic impedance, small signal S parameters and equivalent RLGC transmission line model parameters [:cite:`marks_general_1992`, :cite:`rizzi_microwave_1988`]. 


.. autoclass:: imodulator.RFSimulator.RFSimulatorFEMWELL
   :members:
   :special-members: __init__


References
----------

.. bibliography::
   :filter: docname in docnames
   :style: plain
