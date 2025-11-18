from __future__ import annotations

from pint import UnitRegistry

from matplotlib import pyplot as plt
from shapely.geometry import (
    Polygon,
    MultiPolygon,
    LineString,
    MultiLineString,
    LinearRing,
    box,
)

from collections import OrderedDict

import numpy as np

from imodulator.PhotonicPolygon import (
    SemiconductorPolygon,
    MetalPolygon,
    InsulatorPolygon,
)


PhotonicPolygon = SemiconductorPolygon | MetalPolygon | InsulatorPolygon
Line = LineString | MultiLineString | LinearRing


class PhotonicDevice:
    """
    Base class for a ``PhotonicDevice``. This class holds all the information needed to do a full simulation of an electro-optical modulator.

    """

    def __init__(self, photo_polygons: list[PhotonicPolygon]):
        """
        Initializes the class.

        Based on the resolutions of the mesh settings given to each :class:`PhotoPolygon`, it will also look for the different interfaces between semiconductor, insulator and metal polygons and generate the junction entities. Furthermore, it will also create the line entities based on the ``PhotoPolygon.calculate_current`` and ``PhotoPolygon.d_buffer_current`` attributes.

        .. warning::
            The algorithm to search for junctions is not active at the moment. It is a feature that will be implemented in the future.

        Args:
            photo_polygons: a ``list`` of :class:`PhotonicPolygon` objects. These objects will be used to create the mesh and the dielectric tensors. Note that the polygons in this list will dictate the hierarchy of the polygons for the meshing algorithm. The first polygon in the list will have the highest priority.


        """
        self.reg = UnitRegistry()

        self.e = 1.602176634e-19 * self.reg.coulomb
        self.e0 = 8.854e-12 * self.reg.farad * self.reg.meter**-1
        self.c = 3e8 * self.reg.meter * self.reg.second**-1  # m s^-1
        self.mu0 = (
            4 * np.pi * 1e-7 * self.reg.henry / self.reg.meter
        )  # vacuum magnetic permeability

        self.photo_polygons = photo_polygons #this being a list makes me lose my mind
        self.line_entities = (OrderedDict())  # This will handle the line integral polygons
        self.polygon_entities = OrderedDict()  # This will handle the bulk polygons
        self.junction_entities = (OrderedDict())  # This will hand the SC-Sc and SC-M junctions
        self.resolutions_rf = dict()        #these are the mesh resolutions for rf
        self.resolutions_optical = dict()   #these are the mesh resolutions for optical
        self.resolutions_charge = dict()    #these are the mesh resolutions for charge
        self.resolutions_eo = dict()        #these are the mesh resolutions for electro-optical
        self.mode = dict() #holds interpolators for mode fields
        self.charge = dict() #holds charge outputs
        
        polygon_names = [polygon.name for polygon in self.photo_polygons]
    
        if "background" not in polygon_names:
            # Find the bounding box of all the polygons if it does not have a background polygon
            for i, polygon in enumerate(self.photo_polygons):
                bounds = polygon.polygon.bounds
                if i == 0:
                    xmin, ymin, xmax, ymax = bounds
                else:
                    if bounds[0] < xmin:
                        xmin = bounds[0]
                    if bounds[1] < ymin:
                        ymin = bounds[1]
                    if bounds[2] > xmax:
                        xmax = bounds[2]
                    if bounds[3] > ymax:
                        ymax = bounds[3]

            self.photo_polygons.append(InsulatorPolygon(
                box(xmin, ymin, xmax, ymax).buffer(0.1*max(xmax-xmin, ymax-ymin), 
                                                   join_style = 'bevel',
                                                   cap_style = 'square'),
                name="background",
                optical_material=1.0,
                rf_eps=1.0,
                eo_mesh_settings={"resolution": 100, "SizeMax": 100, "distance": 0.1},
                charge_mesh_settings={"resolution": 100, "SizeMax": 100, "distance": 0.1},
                rf_mesh_settings={"resolution": 100, "SizeMax": 100, "distance": 0.1},
                optical_mesh_settings={"resolution": 100, "SizeMax": 100, "distance": 0.1},
            ))

        for polygon in self.photo_polygons:

            self.polygon_entities[polygon.name] = polygon.polygon

            #Transfer the mesh settings to the resolutions dictionaries
            self.resolutions_rf[polygon.name] = {
                key: value
                for key, value in polygon.rf_mesh_settings.items()
                if key in ["resolution", "SizeMax", "distance"]
            }
            
            self.resolutions_charge[polygon.name] = {
                key: value
                for key, value in polygon.charge_mesh_settings.items()
                if key in ["resolution", "SizeMax", "distance"]
            }

            self.resolutions_optical[polygon.name] = {
                key: value
                for key, value in polygon.optical_mesh_settings.items()
                if key in ["resolution", "SizeMax", "distance" 
                        #    ,"dx" ,"dy" #in case we want seperate 
                           ]
            }

            self.resolutions_eo[polygon.name] = {
                key: value
                for key, value in polygon.eo_mesh_settings.items()
                if key in ["resolution", "SizeMax", "distance"]
            }

            if polygon.calculate_current:

                line = LineString(
                    polygon.polygon.buffer(
                        polygon.d_buffer_current, join_style="bevel"
                    ).exterior
                )
                self.line_entities[polygon.name + "line_current"] = line

                #Transfer the mesh settings to the resolutions dictionaries
                # This one is only necessary for the RF simulator
                self.resolutions_rf[polygon.name + "line_current"] = {
                    key: value
                    for key, value in polygon.rf_mesh_settings.items()
                    if key in ["resolution", "SizeMax", "distance"]
                }

        # # Search for junctions without double counting
        # for i, polygon1 in enumerate(self.photo_polygons):
        #     for polygon2 in self.photo_polygons[i:]:
        #         if isinstance(polygon1, SemiconductorPolygon) and isinstance(
        #             polygon2, SemiconductorPolygon
        #         ):
        #             if polygon1.polygon.touches(polygon2.polygon):

        #                 # Find in which direction the SC polygon is relative to the metal
        #                 p1 = polygon1.polygon.centroid
        #                 p2 = polygon2.polygon.centroid

        #                 vec = p2 - p1
        #                 dx, dy = vec.x, vec.y

        #                 angle = np.arctan2(dy, dx)
        #                 direction = np.sign(np.sin(angle))

        #                 boundary = polygon1.polygon.intersection(polygon2.polygon)
        #                 junction_up = boundary.buffer(
        #                     direction * polygon1.rf_mesh_settings["distance_junction"],
        #                     single_sided=True,
        #                 )
        #                 junction_down = boundary.buffer(
        #                     -direction * polygon2.rf_mesh_settings["distance_junction"],
        #                     single_sided=True,
        #                 )

        #                 if not junction_up.is_empty:
        #                     self.junction_entities[
        #                         f"junction_up_{polygon1.name}_{polygon2.name}"
        #                     ] = junction_up

        #                     #This part is to transfer the mesh settings to the resolutions dictionaries
        #                     self.resolutions_rf[
        #                         f"junction_up_{polygon1.name}_{polygon2.name}"
        #                     ] = {
        #                         "resolution": min(
        #                             polygon1.rf_mesh_settings["resolution_junction"],
        #                             polygon2.rf_mesh_settings["resolution_junction"],
        #                         ),
        #                         "SizeMax": min(
        #                             polygon1.rf_mesh_settings["resolution"],
        #                             polygon2.rf_mesh_settings["resolution"],
        #                         ),
        #                         "distance": 0,
        #                     }

        #                     self.resolutions_optical[
        #                         f"junction_up_{polygon1.name}_{polygon2.name}"
        #                     ] = {
        #                         "resolution": min(
        #                             polygon1.optical_mesh_settings["resolution_junction"],
        #                             polygon2.optical_mesh_settings["resolution_junction"],
        #                         ),
        #                         "SizeMax": min(
        #                             polygon1.optical_mesh_settings["resolution"],
        #                             polygon2.optical_mesh_settings["resolution"],
        #                         ),
        #                         "distance": 0,
        #                     }

        #                     self.resolutions_eo[
        #                         f"junction_up_{polygon1.name}_{polygon2.name}"
        #                     ] = {
        #                         "resolution": min(
        #                             polygon1.eo_mesh_settings["resolution_junction"],
        #                             polygon2.eo_mesh_settings["resolution_junction"],
        #                         ),
        #                         "SizeMax": min(
        #                             polygon1.eo_mesh_settings["resolution"],
        #                             polygon2.eo_mesh_settings["resolution"],
        #                         ),
        #                         "distance": 0,
        #                     }
                            

        #                 if not junction_down.is_empty:
        #                     self.junction_entities[
        #                         f"junction_down_{polygon1.name}_{polygon2.name}"
        #                     ] = junction_down

        #                     #This part is to transfer the mesh settings to the resolutions dictionaries
        #                     self.resolutions_rf[
        #                         f"junction_down_{polygon1.name}_{polygon2.name}"
        #                     ] = {
        #                         "resolution": min(
        #                             polygon1.rf_mesh_settings["resolution_junction"],
        #                             polygon2.rf_mesh_settings["resolution_junction"],
        #                         ),
        #                         "SizeMax": min(
        #                             polygon1.rf_mesh_settings["resolution"],
        #                             polygon2.rf_mesh_settings["resolution"],
        #                         ),
        #                         "distance": 0,
        #                     }

        #                     self.resolutions_optical[
        #                         f"junction_down_{polygon1.name}_{polygon2.name}"
        #                     ] = {
        #                         "resolution": min(
        #                             polygon1.optical_mesh_settings["resolution_junction"],
        #                             polygon2.optical_mesh_settings["resolution_junction"],
        #                         ),
        #                         "SizeMax": min(
        #                             polygon1.optical_mesh_settings["resolution"],
        #                             polygon2.optical_mesh_settings["resolution"],
        #                         ),
        #                         "distance": 0,
        #                     }

        #                     self.resolutions_eo[
        #                         f"junction_down_{polygon1.name}_{polygon2.name}"
        #                     ] = {
        #                         "resolution": min(
        #                             polygon1.eo_mesh_settings["resolution_junction"],
        #                             polygon2.eo_mesh_settings["resolution_junction"],
        #                         ),
        #                         "SizeMax": min(
        #                             polygon1.eo_mesh_settings["resolution"],
        #                             polygon2.eo_mesh_settings["resolution"],
        #                         ),
        #                         "distance": 0,
        #                     }

        #         if (
        #             isinstance(polygon1, SemiconductorPolygon)
        #             and isinstance(polygon2, MetalPolygon)
        #         ) or (
        #             isinstance(polygon1, MetalPolygon)
        #             and isinstance(polygon2, SemiconductorPolygon)
        #         ):
        #             metal_polygon = (
        #                 polygon1 if isinstance(polygon1, MetalPolygon) else polygon2
        #             )
        #             SC_polygon = (
        #                 polygon1
        #                 if isinstance(polygon1, SemiconductorPolygon)
        #                 else polygon2
        #             )

        #             # Find in which direction the SC polygon is relative to the metal
        #             p1 = metal_polygon.polygon.centroid
        #             p2 = SC_polygon.polygon.centroid

        #             vec = p2 - p1
        #             dx, dy = vec.x, vec.y

        #             angle = np.arctan2(dy, dx)
        #             direction = np.sign(np.sin(angle))

        #             if polygon1.polygon.touches(polygon2.polygon):
        #                 boundary = polygon1.polygon.intersection(polygon2.polygon)
        #                 junction = boundary.buffer(
        #                     direction * polygon1.rf_mesh_settings["distance_junction"],
        #                     single_sided=True,
        #                 )

        #                 if not junction.is_empty:
        #                     self.junction_entities[
        #                         f"junction_{polygon1.name}_{polygon2.name}"
        #                     ] = junction

        #                     self.resolutions_eo[
        #                         f"junction_{polygon1.name}_{polygon2.name}"
        #                     ] = {
        #                         "resolution": min(
        #                             polygon1.eo_mesh_settings["resolution_junction"],
        #                             polygon2.eo_mesh_settings["resolution_junction"],
        #                         ),
        #                         "SizeMax": min(
        #                             polygon1.eo_mesh_settings["resolution"],
        #                             polygon2.eo_mesh_settings["resolution"],
        #                         ),
        #                         "distance": 0,
        #                     }

        #                     self.resolutions_rf[
        #                         f"junction_{polygon1.name}_{polygon2.name}"
        #                     ] = {
        #                         "resolution": min(
        #                             polygon1.rf_mesh_settings["resolution_junction"],
        #                             polygon2.rf_mesh_settings["resolution_junction"],
        #                         ),
        #                         "SizeMax": min(
        #                             polygon1.rf_mesh_settings["resolution"],
        #                             polygon2.rf_mesh_settings["resolution"],
        #                         ),
        #                         "distance": 0,
        #                     }

        #                     self.resolutions_optical[
        #                         f"junction_{polygon1.name}_{polygon2.name}"
        #                     ] = {
        #                         "resolution": min(
        #                             polygon1.optical_mesh_settings["resolution_junction"],
        #                             polygon2.optical_mesh_settings["resolution_junction"],
        #                         ),
        #                         "SizeMax": min(
        #                             polygon1.optical_mesh_settings["resolution"],
        #                             polygon2.optical_mesh_settings["resolution"],
        #                         ),
        #                         "distance": 0,
        #                     }

        self.entities = OrderedDict(
            list(self.line_entities.items())
            + list(self.junction_entities.items())
            + list(self.polygon_entities.items())
        )

        # self.mesh = None

        self.epsilon_optical = None
        self.epsilon_rf = None

        self.charge_transport_sim_flag = (
            False  # THIS MUST BE CHANGED TO TRUE ONCE THE CHARGE TRANSPORT SIM IS DONE
        )
        self.mode_sim_flag = (
            False  # THIS MUST BE CHANGED TO TRUE ONCE THE MODE TRANSPORT SIM IS DONE
        )

    def plot_polygons(
        self,
        color_polygon: str = "black",
        color_line: str = "green",
        color_junctions: str = "blue",
        fill_polygons: bool = False,
        poly_list_color: Optional[Dict[str, str]] = None,
        fig=None,
        ax=None,
    ):
        """
        Plots the polygons of the :class:`PhotonicDevice` object.
    
        Args:
            color_polygon: The color to use for the polygons.
            color_line: The color to use for the lines.
            color_junctions: The color to use for the junctions.
            fill_polygons: If ``True``, the polygons will be filled.
            poly_list_color: A dictionary to limit the plotted polygons and to map polygon names to fill colors.
            {"name of the polygon":color,...}
            fig: The figure to plot on. If ``None``, a new figure is created.
            ax: The axis to plot on. If ``None``, a new axis is created.
    
        Returns:
            fig, ax: The figure and axis objects.
        """
        if fig is None and ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
    
        for name, poly in self.entities.items():
            if isinstance(poly, Polygon):
                ax.plot(
                        *poly.exterior.xy,
                        color=color_polygon if "junction" not in name else color_junctions,
                        )
                
                if fill_polygons:
                    if poly_list_color is not None:
                        if name in poly_list_color:
                            fill_color = (
                                poly_list_color.get(name, np.random.rand(3,))
                                if poly_list_color  
                                else np.random.rand(3,)
                            )
                            ax.fill(
                                *poly.exterior.xy,
                                color=fill_color,
                                alpha=0.5,
                            )
                    else:
                        fill_color = np.random.rand(3,)
                        ax.fill(
                            *poly.exterior.xy,
                            color=fill_color,
                            alpha=0.5,
                        )
            elif isinstance(poly, Line):
                ax.plot(*poly.xy, color=color_line)
                
            



if __name__ == "__main__":
    test_Device = PhotonicDevice(1)
