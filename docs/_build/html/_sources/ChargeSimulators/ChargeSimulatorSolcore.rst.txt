.. _solcore: http://docs.solcore.solar/en/develop/ 

ChargeSimulatorSolcore
======================

.. _openbandparams: https://github.com/duarte-jfs/openbandparams

The `ChargeSimulatorSolcore` is a 1D poisson-drift-diffusion simulator based on the open-source software Solcore_, made compatible with the openbandparams_ database. 

.. warning::
   There are currently issues with the solcore simulator, namely:

   - The voltages applied to the device that are returned to the user are not necessarily correct. Depending on the voltage values array you give it it may return the values in a different order. Please see https://github.com/qpv-research-group/solcore5/issues/294

.. figure:: imgs/charge_transport_1D_interpolation.png
   :width: 80%
   :align: center

   Illustration of the 1D to 2D data interpolation performed with the `ChargeSimmulatorSolcore`.

.. autoclass:: imodulator.ChargeSimulator.ChargeSimulatorSolcore
   :members:
   :special-members: __init__


References
----------

.. bibliography::
   :filter: docname in docnames
   :style: plain
