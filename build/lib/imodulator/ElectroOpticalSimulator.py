"""
The role of the ElectroOpticalSimulator is to:

- Generate delta_epsilon_optical;
- Perform electro_optical calculations;

"""

from __future__ import annotations
from typing import Optional

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from shapely.geometry import (
    Polygon,
    MultiPolygon,
    LineString,
    MultiLineString,
    LinearRing,
)
from shapely.ops import clip_by_rect

from femwell.mesh import mesh_from_OrderedDict

from skfem.io.meshio import from_meshio
from skfem.visuals.matplotlib import draw_mesh2d
from skfem import Basis, ElementTriP1, ElementVector, ElementDG, Functional
from skfem.helpers import inner, cross



import numpy as np
import copy
from pint import Quantity

from imodulator import PhotonicDevice
from imodulator.ElectroOpticalModel import ElectroOpticalModel

from imodulator.PhotonicPolygon import (
    SemiconductorPolygon,
    MetalPolygon,
    InsulatorPolygon,
)

from collections import OrderedDict
import warnings

PhotonicPolygon = SemiconductorPolygon | MetalPolygon | InsulatorPolygon
Line = LineString | MultiLineString | LinearRing


class ElectroOpticalSimulator:

    def __init__(
            self,
            device: PhotonicDevice,
            simulation_window: Polygon | None = None,
    ):
        """
        Initializes the ElectroOpticalSimulator with a PhotonicDevice and an optional simulation window.

        The simulation windows will define the electro optical simulation region. 

        Args:
            device (PhotonicDevice): The photonic device to simulate.
            simulation_window (Polygon, optional): The polygon defining the simulation region. If None, the entire device is used.
        """
        
        self.photodevice = device
        self.reg = self.photodevice.reg

        self.e = 1.602176634e-19 * self.reg.coulomb
        self.e0 = 8.854e-12 * self.reg.farad * self.reg.meter**-1
        self.c = 3e8 * self.reg.meter * self.reg.second**-1  # m s^-1
        self.mu0 = (
            4 * np.pi * 1e-7 * self.reg.henry / self.reg.meter
        )  # vacuum magnetic permeability

        self.simulation_window = simulation_window

        self.eo_photopolygons = copy.deepcopy(self.photodevice.photo_polygons)

        self.line_entities = OrderedDict()
        self.polygon_entities = OrderedDict()
        self.junction_entities = OrderedDict()
        self.resolutions = dict()

        if simulation_window is None:
            for name, poly in self.photodevice.junction_entities.items():
                self.junction_entities[name] = poly

        elif simulation_window is not None:
            #Enforce the simulation plane to be a rectangle
            if not np.isclose(
                simulation_window.minimum_rotated_rectangle.area, simulation_window.area
            ):
                raise ValueError("symmetry plane must be a rectangle")

            # Cut all photopolygons by the simulation plane
            idxs_to_pop = []
            for i, poly in enumerate(self.eo_photopolygons):
                if poly.polygon.intersects(
                    simulation_window
                ) and not simulation_window.contains(poly.polygon):
                    poly_tmp = clip_by_rect(poly.polygon, *simulation_window.bounds)

                    if poly_tmp.is_empty:
                        idxs_to_pop.append(i)
                    else:
                        self.eo_photopolygons[i].polygon = poly_tmp

                elif not poly.polygon.intersects(simulation_window):
                    idxs_to_pop.append(i)

            for index in sorted(idxs_to_pop, reverse=True):
                del self.eo_photopolygons[index]

            # Select all junctions cutted by the symmetry plane
            for name, poly in self.photodevice.junction_entities.items():

                if poly.intersects(simulation_window) and not simulation_window.contains(
                    poly
                ):
                    poly_tmp = clip_by_rect(polygon, *simulation_window.bounds)

                    if not poly_tmp.is_empty:
                        self.junction_entities[name] = poly_tmp

                elif simulation_window.contains(poly):
                    self.junction_entities[name] = poly

        for polygon in self.eo_photopolygons:
            self.polygon_entities[polygon.name] = polygon.polygon
        # We now have all the photopolygons cut by the plane. Let us finally add the boundaries

        surf_bounds = self.polygon_entities["background"].bounds

        left = LineString(
            [(surf_bounds[0], surf_bounds[1]), (surf_bounds[0], surf_bounds[3])]
        )

        bottom = LineString(
            [(surf_bounds[0], surf_bounds[1]), (surf_bounds[2], surf_bounds[1])]
        )

        right = LineString(
            [(surf_bounds[2], surf_bounds[1]), (surf_bounds[2], surf_bounds[3])]
        )

        top = LineString(
            [(surf_bounds[0], surf_bounds[3]), (surf_bounds[2], surf_bounds[3])]
        )

        self.line_entities["left"] = left
        self.line_entities["bottom"] = bottom
        self.line_entities["right"] = right
        self.line_entities["top"] = top

        self.entities = OrderedDict(
            list(self.line_entities.items())
            + list(self.junction_entities.items())
            + list(self.polygon_entities.items())
        )



        # Transfer the resolutions from the photonic device to the optical_simulator
        for name in self.entities.keys():
            if name in self.photodevice.resolutions_optical.keys():
                self.resolutions[name] = self.photodevice.resolutions_eo[name]

    def make_mesh(
        self,
        default_resolution_min: float = 1e-12,
        default_resolution_max: float = 100,
        filename: Optional[str] = None,
        gmsh_algorithm: int = 5,
        global_quad: Optional[bool] = False,
        verbose: bool = False,
        mesh_scaling_factor: float = 1.0,
    ):
        """
        Returns a gmsh msh with physicals corresponding to the self.entities

        """
        self.mesh = from_meshio(
            mesh_from_OrderedDict(
                self.entities,
                self.resolutions,
                default_resolution_max=default_resolution_max,
                default_resolution_min=default_resolution_min,
                filename=filename,
                gmsh_algorithm=gmsh_algorithm,
                global_quad=global_quad,
                verbose=verbose,
                mesh_scaling_factor=mesh_scaling_factor,
            )
        )
        # Choosing each element as ElementTriP1 is crucial as the dofs that assumes are the mesh vertices
        self.basis = Basis(self.mesh, ElementTriP1(), intorder=4)

    def plot_polygons(
        self,
        color_polygon="black",
        color_line="green",
        color_junctions="blue",
        fig=None,
        ax=None,
    ):
        "plots all the polygons and boundaries"
        if fig is None and ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        for name, poly in self.entities.items():
            if isinstance(poly, Polygon):
                ax.plot(
                    *poly.exterior.xy,
                    color=color_polygon if "junction" not in name else color_junctions,
                )
            elif isinstance(poly, Line):
                ax.plot(*poly.xy, color=color_line)

    def plot_mesh(
        self,
        plot_polygons: bool = True,
    ):
        """
        Plots the mesh
        """
        ax = draw_mesh2d(self.mesh)
        ax.set_axis_on()

        fig = ax.get_figure()

        if plot_polygons:
            self.plot_polygons(color_polygon="red", fig=fig, ax=ax)

        return fig, ax
    
    def get_epsilon_optical(
        self,
    ):
        """
        This function will calculate the CHANGE in epsilon optical due to the electro optical effects.

        By default, it will return all the contributors for the change in refractive index and corresponding labels. The output will follow this structure:

        epsilon_optical = {
            ElectroOpticalModelName: {
                'dperms': np.zeros((3, 3, self.mesh.nvertices, N_bias_points, n_effects), dtype=np.complex128),
                'labels': [label1, label2, ...]
            },

            ElectroOpticalModelName2: {
                'dperms': np.zeros((3, 3, self.mesh.nvertices, N_bias_points, n_effects), dtype=np.complex128),
                'labels': [label1, label2, ...]
        }

        """

        if self.mesh is None:
            raise RuntimeError(
                "Cannot attribute an optical epsilon. The mesh of the optical simulator is not yet generated"
            )


        N_bias_points = len(self.photodevice.charge['V'])


        # Account for the case where the user wants tot epsilon optical to be fractioned into the various effects that the electro optical models gives
        epsilon_optical_all = {}

        #Loop over all the polygons and find all the electro optic models present
        for photo_polygon in self.eo_photopolygons:
            # print(photo_polygon.electro_optic_module)
            # print(issubclass(photo_polygon.electro_optic_module,ElectroOpticalModel))
            if hasattr(photo_polygon, "electro_optic_module"):
                if photo_polygon.electro_optic_module is not None:
                    if issubclass(photo_polygon.electro_optic_module, ElectroOpticalModel):
                        if photo_polygon.electro_optic_module.__name__ not in epsilon_optical_all:
                            epsilon_optical_all[photo_polygon.electro_optic_module.__name__] = {
                                'dperms': np.zeros(
                                    (3, 3, self.mesh.nvertices, N_bias_points, photo_polygon.electro_optic_module.n_effects), dtype=np.complex128
                                ),
                                'labels': []
                            }
                        

        # the self.photo_polygons is created so that idx 0 has higher priority over idx 1
        # Here, however, if we loop through the photo_polygons from idx 0 to idx N
        # the hierarchy on the boundaries will be inverted. That is,
        #lower lying polygons in hierarchy will dominate the boundaries. Therefore, we need to
        # loop over the inverted list of photo_polygons

        for photo_polygon in self.eo_photopolygons[::-1]:


            if hasattr(photo_polygon, "electro_optic_module"):
                if photo_polygon.electro_optic_module is not None:
                    if issubclass(photo_polygon.electro_optic_module, ElectroOpticalModel):
                        

                        elements_idxs = self.mesh.subdomains[photo_polygon.name]
                        triangles = self.mesh.t[:, elements_idxs]
                        vertices_idxs = np.unique(triangles.flatten())

                        # Loop through all the voltages
                            
                        EO_model_name = photo_polygon.electro_optic_module.__name__

                        for voltage_idx in range(N_bias_points):

                            x = self.mesh.p[0, vertices_idxs]
                            y = self.mesh.p[1, vertices_idxs]

                            mun = self.photodevice.charge['mun'][voltage_idx](x,y)
                            mup = self.photodevice.charge['mup'][voltage_idx](x,y)
                            Ec = self.photodevice.charge['Ec'][voltage_idx](x,y)
                            Ev = self.photodevice.charge['Ev'][voltage_idx](x,y)
                            Efp = self.photodevice.charge['Efp'][voltage_idx](x,y)
                            Efn = self.photodevice.charge['Efn'][voltage_idx](x,y)
                            P = self.photodevice.charge['P'][voltage_idx](x,y)
                            N = self.photodevice.charge['N'][voltage_idx](x,y)
                            Efield = self.photodevice.charge['Efield'][voltage_idx](x,y)

                            data = photo_polygon.electro_optic_module(
                                mup=mup,
                                mun=mun,
                                Ec=Ec,
                                Ev=Ev,
                                Efp=Efp,
                                Efn=Efn,
                                P=P,
                                N=N,
                                Efield=Efield,
                                reg=self.photodevice.reg,
                                **photo_polygon.electro_optic_module_kwargs,
                            ).get_dperm(fractions = True)

                            labels = data[0]
                            perms = data[1:]



                            for perm_idx, perm in enumerate(perms):
                                epsilon_optical_all[EO_model_name]['dperms'][:, :, vertices_idxs, voltage_idx, perm_idx] = perm
                            
                            epsilon_optical_all[EO_model_name]['labels'] = labels

                                
            
        self.epsilon_optical = epsilon_optical_all

    def _calculate_eo_integrand(
            self,
            voltage_idx: int,
            rot_x: float = 0,
            rot_y: float = 0,
            rot_z: float = 0,
            base_epsilon_voltage_idx: int = 0,
            optical_mode: str = 'TE'
    ):
        
        """
        This function is used to calculate the integrands for the change in propagation constant due to change in epsilon optical.

        Mathematically it folloes integrand = Eopt*dn*Eopt^H, where dn is the change in epsilon optical due to the electro optical effects.

        The output from this follows the same structure as the epsilon_optical dictionary, that is:

        integrands = {
            ElectroOpticalModelName: {
                'integrand': np.zeros((self.mesh.nvertices, N_effects), dtype=np.complex128),
                'labels': [label1, label2, ...]
            },  

            ElectroOpticalModelName2: {
                'integrand': np.zeros((self.mesh.nvertices, N_effects), dtype=np.complex128),
                'labels': [label1, label2, ...]
            }
        }  

        """

        def Rz(alpha):
            mat=np.array([[np.cos(alpha), -np.sin(alpha), 0],
                        [np.sin(alpha), np.cos(alpha), 0],
                        [0,0,1]])
            return mat

        def Ry(beta):
            mat=np.array([[np.cos(beta), 0, np.sin(beta)],
                        [0           , 1, 0],
                        [-np.sin(beta),0,np.cos(beta)]])
            return mat

        def Rx(gamma):
            mat=np.array([[1           , 0            , 0],
                        [0           , np.cos(gamma), -np.sin(gamma)],
                        [0           , np.sin(gamma),np.cos(gamma)]])
            return mat

        ALPHA=rot_z
        BETA=rot_y
        GAMMA=rot_x

        rot= Rx(ALPHA) @ Ry(BETA) @ Rz(GAMMA)

        #Generate the optical field
        x = self.mesh.p[0, :]
        y = self.mesh.p[1, :]

        self.Eopt = np.asarray([
            self.photodevice.mode[optical_mode]['Ex']((x, y)),
            self.photodevice.mode[optical_mode]['Ey']((x, y)),
            self.photodevice.mode[optical_mode]['Ez']((x, y))
        ])

        self.Hopt = np.asarray([
            self.photodevice.mode[optical_mode]['Hx']((x, y)),
            self.photodevice.mode[optical_mode]['Hy']((x, y)),
            self.photodevice.mode[optical_mode]['Hz']((x, y))
        ])
        
            
        integrands = {}

        for electro_optic_model_name, data in self.epsilon_optical.items():
            dperms = data['dperms']
            labels = data['labels']

            
            integrands_tmp = np.zeros((self.mesh.nvertices, len(labels)), dtype=np.complex128)
            for dperm_idx in range(dperms.shape[-1]): #The last index is the neffects index
                base_epsilon = dperms[:, :, :, base_epsilon_voltage_idx, dperm_idx] # shape (3,3,self.mesh.nvertices)

                rot_epsilon = (np.einsum('ir,rkl,kj->ijl', rot, dperms[:,:,:,voltage_idx, dperm_idx], rot.T) - 
                            np.einsum('ir,rkl,kj->ijl', rot, base_epsilon, rot.T))

                integrands_tmp[:,dperm_idx] = np.einsum('ik,ijk,jk->k', self.Eopt.conjugate(), rot_epsilon, self.Eopt)

            integrands[electro_optic_model_name] = {
                'integrand': integrands_tmp,
                'labels': labels
            }

        return integrands

    def calculate_EO_response(
        self,
        voltage_idx: int,
        rot_x: float = 0,
        rot_y: float = 0,
        rot_z: float = 0,
        base_epsilon_voltage_idx: int = 0,
        optical_mode: str = 'TE',
        wavelength: float = 1.55
    ):
        """
        This function calculates the electro optical response of the device.
        It calculates the change in propagation constant due to the electro optical effects.
        The output will follow the same structure as the epsilon_optical dictionary, that is:

        results = {
            ElectroOpticalModelName: {
                'results': np.zeros((N_effects), dtype=np.complex128),
                'labels': [label1, label2, ...]
            },

            ElectroOpticalModelName2: {
                'results': np.zeros((N_effects), dtype=np.complex128),
                'labels': [label1, label2, ...]
            }
        }

        Args:
            voltage_idx (int): The index of the voltage to use for the electro optical effects.
            rot_x (float): The rotation around the x-axis in radians.
            rot_y (float): The rotation around the y-axis in radians.
            rot_z (float): The rotation around the z-axis in radians.
            base_epsilon_voltage_idx (int): The index of the voltage to use as the base epsilon optical.
            optical_mode (str): The optical mode to use for the electro optical effects. Default is 'TE'.
            wavelength (float): The wavelength of the optical mode in micrometers. Default is 1.55 micrometers.
        Returns:
            results (dict): A dictionary containing the results of the electro optical effects.
            The keys are the names of the electro optical models and the values are dictionaries with the results and labels.
            
        """
        
        @Functional(dtype=np.complex128)
        def integral_form(w):
            return w.integrand
        

        integrands = self._calculate_eo_integrand(
            voltage_idx, rot_x, rot_y, rot_z, base_epsilon_voltage_idx, optical_mode)
        
        I_integrand = (np.cross(self.Eopt.T.conjugate(), self.Hopt.T) + np.cross(self.Eopt.T, self.Hopt.T.conjugate()))[:,2]
        power = integral_form.assemble(self.basis, integrand = I_integrand)

        #Loop over all the integrands and remove any nan values that may be present

        for electro_optic_model_name, data in integrands.items():
            integrand = data['integrand']
            integrand = np.nan_to_num(integrand, nan=0)
        
        results = {}

        for electro_optic_model_name, data in integrands.items():

            results_tmp = np.zeros(data['integrand'].shape[-1], dtype=np.complex128)

            for idx in range(data['integrand'].shape[-1]): #Loop over the effects
                
                integ = data['integrand'][:, idx]

                # print('INSIDE THE SIM', integ)
                integral_result = integral_form.assemble(self.basis, integrand = integ/np.abs(power)) 
                integral_result *= 2 * np.pi * 3e8 / (wavelength*1e-6) * 8.85e-12

                results_tmp[idx] = integral_result

            if electro_optic_model_name not in results:
                results[electro_optic_model_name] = {}

            results[electro_optic_model_name]['results'] = results_tmp
            results[electro_optic_model_name]['labels'] = data['labels']
            
        return results