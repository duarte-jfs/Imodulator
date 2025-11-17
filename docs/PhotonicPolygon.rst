Photonic Polygons
=================

.. _shapely: https://shapely.readthedocs.io/en/stable/

A ``PhotonicPolygon`` is a dataclass object that is meant to work as a subdomain for a ``PhotonicDevice``, meaning that it should hold all the information about a single part of your photonic device. In particular it should hold the information about:

- The geometry of the polygon, which need not be a rectangle.
- The material properties of the polygon, including the optical and the RF relative electrical permitivity;
- The various mesh settings for optical, RF, charge transport and electro-optical calculations;
- Additional charge transport and optical simulator arguments.

This is a list of the ``PhotonicPolygon`` classes available in the ``imodulator`` package.	

.. autoclass:: imodulator.PhotonicPolygon.SemiconductorPolygon
   :members:

.. autoclass:: imodulator.PhotonicPolygon.MetalPolygon
   :members:

.. autoclass:: imodulator.PhotonicPolygon.InsulatorPolygon
   :members: