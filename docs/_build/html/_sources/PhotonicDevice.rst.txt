Photonic Device
===============

.. _FEMWELL: https://femwell.readthedocs.io

.. _MESHWELL: https://simbilod.github.io/meshwell/intro.html

.. _sk-fem: https://scikit-fem.readthedocs.io/en/latest/

The ``PhotonicDevice`` is at the core of the ``imodulator`` package. It is a data centre for the modulator, meant to store information about the geometry and materials of the device. It is built in such a way to allow seamless integration with the various solvers available in the package, in order to minimize the amount of time spent on parametrization between different simulations. 

.. autoclass:: imodulator.PhotonicDevice.PhotonicDevice
   :members:
   :special-members: __init__


Some remarks about the creation of a ``PhotonicDevice``:

1. It is reccomended to always create a ``PhotonicPolygon`` called 'background'. If it is not given one with such a name, it shall be created one and added to the photopolygons list. This is necessary so as to be able to mesh your background.
2. The ``PhotonicDevice`` will hold all the results from both the mode and the charge transport simulations in the form of interpolators. In case you would like to access such information they shall be stored in the ``PhotonicDevice.mode`` and ``PhotonicDevice.charge`` attributes, respectively.

References
----------

.. bibliography::
   :filter: docname in docnames
   :style: plain