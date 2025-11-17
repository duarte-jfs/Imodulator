from __future__ import annotations

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from shapely.geometry import (
    Polygon,
    MultiPolygon,
    LineString,
    MultiLineString,
    LinearRing,
)
from shapely.ops import clip_by_rect, linemerge, unary_union

from skfem.io.meshio import from_meshio
from skfem.visuals.matplotlib import draw_mesh2d
from skfem import Basis, ElementTriP1, ElementVector, ElementDG, Functional
from skfem.helpers import inner

from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import (
    Modes,
    Mode,
    calculate_scalar_product,
    compute_modes,
)

import numpy as np
import copy
from pint import Quantity

from imodulator import PhotonicDevice
from collections import OrderedDict
import warnings

from imodulator.PhotonicPolygon import (
    SemiconductorPolygon,
    MetalPolygon,
    InsulatorPolygon,
)

PhotonicPolygon = SemiconductorPolygon | MetalPolygon | InsulatorPolygon
Line = LineString | MultiLineString | LinearRing


class RFSimulatorFEMWELL:
    """
    Base class for RF simulation of a PhotonicDevice. The role of the RFsimulator is to:

    .. warning::
        The algorithm to search for junctions is not active at the moment. It is a feature that will be implemented in the future.

    * generate the RF mesh;
    * apply RF simmetry planes;
    * Generate field visualizations of the RF modes;
    * Generate mesh visuals;
    * Solve for the RF modes;
    * Return small-signal S parameters;
    * Return RLGC parameters based on waveguide theory;

    """

    def __init__(
        self,
        device: PhotonicDevice,
        simulation_window: Polygon | None = None,
    ):
        """
        Initializes the simulator with the photonic device and the simulation window.

        Args:
            device: The :class:`PhotonicDevice` object to simulate.
            simulation_window: The simulation window to use for the RF simulation. If None, the entire device is simulated.

        .. note::
            The simulation window MUST be a rectangle. If not, the simulation will fail. Also beware that the definition of the simulation window is how you can make use of symmetry planes through the definition of metal boundaries or not. See `compute_modes` for more information on how to use the symmetry planes.
        
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

        self.rf_photopolygons = copy.deepcopy(self.photodevice.photo_polygons)
        self.line_entities = OrderedDict()
        self.polygon_entities = OrderedDict()
        self.junction_entities = OrderedDict()
        self.resolutions = dict()

        # Add the lines for line integrals of currents
        for polygon in self.rf_photopolygons:

            if polygon.calculate_current:

                line = LineString(
                    polygon.polygon.buffer(
                        polygon.d_buffer_current, join_style="bevel"
                    ).exterior
                )
                self.line_entities[polygon.name + "line_current"] = line

        # THIS NEEDS A REFACTOR!!!
        if simulation_window is None:
            # Select all junctions cutted by the simulation window
            for name, poly in self.photodevice.junction_entities.items():
                self.junction_entities[name] = poly

        elif simulation_window is not None:
            if not np.isclose(
                simulation_window.minimum_rotated_rectangle.area, simulation_window.area
            ):
                raise ValueError("simulation window must be a rectangle")

            # Cut all photopolygons by the simulation window
            idxs_to_pop = []
            for i, poly in enumerate(self.rf_photopolygons):
                if poly.polygon.intersects(
                    simulation_window
                ) and not simulation_window.contains(poly.polygon):
                    poly_tmp = clip_by_rect(poly.polygon, *simulation_window.bounds)

                    if poly_tmp.is_empty:
                        idxs_to_pop.append(i)
                    else:
                        self.rf_photopolygons[i].polygon = poly_tmp

                elif not poly.polygon.intersects(simulation_window):
                    idxs_to_pop.append(i)

            for index in sorted(idxs_to_pop, reverse=True):
                del self.rf_photopolygons[index]

            # Cut all the line entities by the simulation window
            keys_to_pop = []
            for key, poly in self.line_entities.items():
                if poly.intersects(simulation_window) and not simulation_window.contains(
                    poly
                ):
                    poly_tmp = clip_by_rect(poly, *simulation_window.bounds)

                    if type(poly_tmp) == MultiLineString:
                        self.line_entities[key] = linemerge(poly_tmp)

                    elif poly_tmp.is_empty:
                        keys_to_pop.append(key)
                    else:
                        self.line_entities[key] = poly_tmp
                elif not poly.intersects(simulation_window):
                    keys_to_pop.append(key)

            for key in keys_to_pop:
                self.line_entities.pop(key)

            # Select all junctions cutted by the simulation window
            for name, poly in self.photodevice.junction_entities.items():

                if poly.intersects(simulation_window) and not simulation_window.contains(
                    poly
                ):
                    poly_tmp = clip_by_rect(polygon, *simulation_window.bounds)

                    if not poly_tmp.is_empty:
                        self.junction_entities[name] = poly_tmp

                elif simulation_window.contains(poly):
                    self.junction_entities[name] = poly

        for polygon in self.rf_photopolygons:
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

        # Transfer the resolutions from the photonic device to the rf_simulator
        for name in self.entities.keys():
            if name in self.photodevice.resolutions_rf.keys():
                self.resolutions[name] = self.photodevice.resolutions_rf[name]

    def make_mesh(
        self,
        default_resolution_min: float = 1e-12,
        default_resolution_max: float = 100,
        filename: str = None,
        gmsh_algorithm: int = 5,
        global_quad: bool = False,
        verbose: bool = False,
        mesh_scaling_factor: float = 1.0,
    ):
        """
        Returns a `gmsh <https://gmsh.info/>`_ with the geometric information of the photonic device. The mesh generation is currently handled by `FEMWELL`_, but it will later be deprecated in favour of `MESHWELL`_.

        Args:
            default_resolution_min: The minimum resolution to use for the mesh generation if no resolution is specified for a given entity.
            default_resolution_max: The maximum resolution to use for the mesh generation if no resolution is specified for a given entity.
            filename: The filename to save the mesh to. If `None`, the mesh is not saved.
            gmsh_algorithm: The algorithm to use for the mesh generation.
            global_quad: Whether to use global quadrature for the mesh generation.
            verbose: Whether to print verbose output during the mesh generation.
            mesh_scaling_factor: The scaling factor to apply to the mesh.

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

    def get_epsilon_rf(self, 
                       frequency: float,
                       use_charge_transport_data: bool = False,
                       voltage_idx: int = 0):
        """
        This function will return the :math:`\epsilon_{RF}` tensor with the signature ``self.epsilon_rf[vertice_idx, voltage_idx]`` where ``vertice_idx`` is the index of the vertice in the mesh and `voltage_idx` is the index of the voltage in the bias points. In case no bias dependent data is available (i.e. no charge transport data is available) the ``voltage_idx`` must be 0 as the second axis of the array will have size 1.

        Args:
            frequency: The frequency at which to compute the permittivity tensor. The frequency must be in GHz.
            use_charge_transport_data: Whether to use the charge transport data to compute the permittivity tensor. Doing so will yield a :math:`\sigma(x,y,V)`. Make sure your mesh is appropriate for it.
            voltage_idx: The index of the voltage to use for the permittivity tensor.

        """
        omega = 2 * np.pi * frequency * self.reg.GHz

        self.epsilon_rf = np.zeros((self.mesh.nvertices, 1), dtype=np.complex128)

        # the self.photo_polygons is created so that idx 0 has higher priority over idx 1
        # Here, however, if we loop through the photo_polygons from idx 0 to idx N
        # the hierarchy on the boundaries will be inverted. That is,
        #lower lying polygons in hierarchy will dominate the boundaries. Therefore, we need to
        # loop over the inverted list of photo_polygons
        for photo_polygon in self.rf_photopolygons[::-1]:
            # Sweep through all the junctions and see which fall within the semiconductor
            extra_vertices = np.asarray([], dtype=int)
            for name, junction_poly in self.junction_entities.items():
                if photo_polygon.polygon.contains(junction_poly):
                    elements_idxs = self.mesh.subdomains[name]
                    triangles = self.mesh.t[:, elements_idxs]
                    vertices_idxs = np.unique(triangles.flatten())
                    extra_vertices = np.concatenate((extra_vertices, vertices_idxs))

            elements_idxs = self.mesh.subdomains[photo_polygon.name]
            triangles = self.mesh.t[:, elements_idxs]
            vertices_idxs = np.unique(triangles.flatten())
            vertices_idxs = np.concatenate((vertices_idxs, extra_vertices))
            
            self.epsilon_rf[vertices_idxs, 0] = photo_polygon.rf_eps(omega)

        if use_charge_transport_data == True:
            # Expand the epsilon_rf array to the voltage values
            # for photo_poly in self.rf_photopolygons:
            #     if isinstance(photo_poly, SemiconductorPolygon):
            #         N_bias_points = len(photo_poly.Ec)
            #         break

            N_bias_points = len(self.photodevice.charge['V'])

            if self.epsilon_rf is None:
                self.epsilon_rf = np.zeros((self.mesh.nvertices, N_bias_points), dtype=np.complex128)
            
            elif self.epsilon_rf.shape[1] == 1:  
                # This is to account for the case where it has been first initialized 
                #with the charge transport data so we extend the array to another dimension
                self.epsilon_rf = (
                    np.tile(self.epsilon_rf[:, 0], N_bias_points)
                    .reshape(N_bias_points, -1)
                    .T
                )

            for voltage_idx in range(N_bias_points):
                # Sweep through all the semiconductors
                # the self.photo_polygons is created so that idx 0 has higher priority over idx 1
                # Here, however, if we loop through the photo_polygons from idx 0 to idx N
                # the hierarchy on the boundaries will be inverted. That is,
                #lower lying polygons in hierarchy will dominate the boundaries. Therefore, we need to
                # loop over the inverted list of photo_polygons
                for photo_polygon in self.rf_photopolygons[::-1]:
                    # Sweep through all the junctions and see which fall within the semiconductor
                    extra_vertices = np.asarray([], dtype=int)
                    for name, junction_poly in self.junction_entities.items():
                        if photo_polygon.polygon.contains(junction_poly):
                            elements_idxs = self.mesh.subdomains[name]
                            triangles = self.mesh.t[:, elements_idxs]
                            vertices_idxs = np.unique(triangles.flatten())
                            extra_vertices = np.concatenate((extra_vertices, vertices_idxs))

                    elements_idxs = self.mesh.subdomains[photo_polygon.name]
                    triangles = self.mesh.t[:, elements_idxs]
                    vertices_idxs = np.unique(triangles.flatten())
                    vertices_idxs = np.concatenate((vertices_idxs, extra_vertices))

                    if (
                        isinstance(photo_polygon, SemiconductorPolygon) and
                        photo_polygon.has_charge_transport_data
                    ):
                        x = self.mesh.p[0, vertices_idxs]
                        y = self.mesh.p[1, vertices_idxs]

                        self.epsilon_rf[vertices_idxs, voltage_idx] = self.epsilon_rf[
                            vertices_idxs, voltage_idx
                        ].real + (
                            -1j
                            * (
                                self.e
                                * (self.photodevice.charge['N'][voltage_idx](x,y) * 
                                   self.photodevice.charge['mun'][voltage_idx](x,y) + 
                                   self.photodevice.charge['P'][voltage_idx](x,y) * 
                                   self.photodevice.charge['mup'][voltage_idx](x,y))
                                / omega
                                / self.e0
                            )
                            .to(self.reg.dimensionless)
                            .magnitude
                        )
                    else:
                        self.epsilon_rf[vertices_idxs, voltage_idx] = photo_polygon.rf_eps(omega)

    def compute_modes(
        self,
        frequency: float = 10,
        voltage_idx: int = 0,
        num_modes: int = 1,
        order: int = 1,
        metallic_boundaries: list | str | bool = False,
        n_guess: float = 4.0,
        return_modes: bool = False,
        use_charge_transport_data: bool = False,
    ) -> Modes:
        """
        Compute the electromagnetic RF modes of a photonic device at a given frequency. The modes are computed via `FEMWELL`_.

        Args:
            frequency: The frequency at which to compute the modes. The frequency must be in GHz.
            voltage_idx: The index of the voltage to use for the permittivity tensor.
            num_modes: The number of modes to compute.
            order: Order of the basis functions to use in the EM solver.
            metallic_boundaries: The boundaries to treat as metallic. If `False`, no boundaries are treated as metallic. If `True`, all boundaries are treated as metallic. If a list of strings, the boundaries with the given names are treated as metallic. At the moment, the simulation window is treated as a square, therefore, the metallic boundaries can be ``left`, ``right``, ``top`` and ``bottom`` boundaries.
            n_guess: Initial guess for the effective index.
            return_modes: Whether to return the computed modes. 
            use_charge_transport_data: Whether to use the charge transport data to compute the permittivity tensor. Doing so will yield a :math:`\sigma(x,y,V)`. Make sure your mesh is appropriate for it.

        Returns:
            The computed modes if `return_modes` is `True`.
       
        """

        self.get_epsilon_rf(frequency, use_charge_transport_data=use_charge_transport_data, voltage_idx=voltage_idx)

        modes = compute_modes(
            self.basis,
            self.epsilon_rf[:, voltage_idx],
            (self.c / (frequency * self.reg.GHz)).to(self.reg.micrometer).magnitude,
            mu_r=1,
            num_modes=num_modes,
            order=order,
            metallic_boundaries=metallic_boundaries,
            n_guess=n_guess,
        )
        

        self.modes = modes.sorted(lambda mode: mode.n_eff.real)

        if return_modes:
            return self.modes

    def plot_eps_rf(
            self,
            frequency: float = 10,
            voltage_idx: int = 0,
            use_charge_transport_data: bool = False,
            log_scale_im: bool = True,
            log_scale_re: bool = True,
            cmap = "jet",
    ):

        """
        Plots the real and imaginary parts of the RF permittivity tensor.
        
        Args:
            frequency: The frequency at which to compute the permittivity tensor. The frequency must be in GHz.
            voltage_idx: The index of the voltage to use for the permittivity tensor.
            use_charge_transport_data: Whether to use the charge transport data to compute the permittivity tensor. Doing so will yield a :math:`\sigma(x,y,V)`.
            log_scale_im: Whether to plot the imaginary part of the permittivity tensor in logarithmic scale.
            log_scale_re: Whether to plot the real part of the permittivity tensor in logarithmic scale.
            cmap: The colormap to use for plotting.
        
        returns:
            fig: The matplotlib figure containing the plots.
            ax1: The axis for the imaginary part of the permittivity tensor.
            ax2: The axis for the real part of the permittivity tensor.
        
        """
        self.get_epsilon_rf(frequency=frequency, use_charge_transport_data=use_charge_transport_data)

        fig = plt.figure(figsize=(10, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        #Plot the polygons
        self.plot_polygons(
            ax=ax1,
        )
        self.plot_polygons(
            ax=ax2,
        )

        #Plot the imaginary part of the permittivity
        data_to_plot_im = self.epsilon_rf[:, voltage_idx].imag
        data_to_plot_re = self.epsilon_rf[:, voltage_idx].real

        if log_scale_im:
            data_to_plot_im = np.log10(-data_to_plot_im)
        if log_scale_re:    
            data_to_plot_re = np.log10(data_to_plot_re)

        self.basis.plot(
            data_to_plot_im,
            cmap=cmap,
            colorbar=True,
            ax=ax1,
        )

        self.basis.plot(
            data_to_plot_re,
            cmap=cmap,
            colorbar=True,
            ax=ax2,
        )

        if use_charge_transport_data:
            ax1.set_title(f"Im(-$\epsilon_{{RF}}$) at {frequency} GHz, V={self.photodevice.charge['V'][voltage_idx]} V")
            ax2.set_title(f"Re($\epsilon_{{RF}}$) at {frequency} GHz, V={self.photodevice.charge['V'][voltage_idx]} V")
        else:
            ax1.set_title(f"Im(-$\epsilon_{{RF}}$) at {frequency} GHz")
            ax2.set_title(f"Re($\epsilon_{{RF}}$) at {frequency} GHz")

        ax1.set_xlabel("x (um)")
        ax1.set_ylabel("y (um)")

        ax2.set_xlabel("x (um)")
        ax2.set_ylabel("y (um)")

        return fig, ax1, ax2

    def plot_mode(
        
        self,
        mode: Mode,
        Nx: int = 50,
        Ny: int = 50,
        xmin: float = 0,
        xmax: float = 10,
        ymin: float = -10,
        ymax: float = +10,
        figsize: tuple[float, float] = (10, 4),
        wspace: float = 0.5,
        color_polygons: str = "white",
        color_integral_lines: str = "red",
        cmap: str = "jet",
        color_vectors: str = "black",
    ):
        """
        Plots the electric and magnetic fields of a mode.

        Args:
            mode: The mode object containing the electric and magnetic field data to be plotted.
            Nx: Number of grid points along the x-axis.
            Ny: Number of grid points along the y-axis.
            xmin: Minimum x-coordinate for the plot. 
            xmax: Maximum x-coordinate for the plot. 
            ymin: Minimum y-coordinate for the plot. 
            ymax: Maximum y-coordinate for the plot. 
            figsize: Figure size as (width, height). 
            wspace: The amount of width reserved for blank space between subplots.
            color_polygons: Color of the polygons in the plot.
            color_integral_lines: Color of the integral lines in the plot. 
            cmap: Colormap used for plotting the magnitude of the fields. 
            color_vectors: Color of the vector field streamlines. 
        Returns:
            None: This method does not return any value. It generates a plot.

        .. note::

            * This method uses the `mode` object to interpolate the electric and magnetic fields onto a rectangular grid.
            * The fields are split into transverse and longitudinal components for plotting.
            * The resulting plot includes both the magnitude of the fields and their respective vector field streamlines.
            * Polygons and integral lines can be overlaid on the plot for additional context.
            * The projection into a rectangular grid is quite heavy, so it is reccomended to use this method only for small regions to avoid large rectangular grids.

        """
        grid_x, grid_y = np.meshgrid(
            np.linspace(xmin, xmax, Nx), np.linspace(ymin, ymax, Ny)
        )

        # Transform the data onto a rectangular regular grid
        grid_data_E = np.zeros((Ny, Nx, 3), dtype=complex)
        grid_data_H = np.zeros((Ny, Nx, 3), dtype=complex)

        basis = mode.basis
        basis_fix = basis.with_element(ElementVector(ElementDG(ElementTriP1())))

        (et, et_basis), (ez, ez_basis) = basis.split(mode.E)
        (et_x, et_x_basis), (et_y, et_y_basis) = basis_fix.split(
            basis_fix.project(et_basis.interpolate(et))
        )

        coordinates = np.array([grid_x.flatten(), grid_y.flatten()])

        grid_data = np.array(
            (
                et_x_basis.interpolator(et_x)(coordinates),
                et_y_basis.interpolator(et_y)(coordinates),
                ez_basis.interpolator(ez)(coordinates),
            )
        ).T

        grid_data_E = grid_data.reshape((*grid_x.shape, -1))

        (et, et_basis), (ez, ez_basis) = basis.split(mode.H)
        (et_x, et_x_basis), (et_y, et_y_basis) = basis_fix.split(
            basis_fix.project(et_basis.interpolate(et))
        )

        coordinates = np.array([grid_x.flatten(), grid_y.flatten()])

        grid_data = np.array(
            (
                et_x_basis.interpolator(et_x)(coordinates),
                et_y_basis.interpolator(et_y)(coordinates),
                ez_basis.interpolator(ez)(coordinates),
            )
        ).T

        grid_data_H = grid_data.reshape((*grid_x.shape, -1))

        # Plot the data

        fig = plt.figure(figsize=figsize)
        gs = GridSpec(1, 2, wspace=wspace)

        ax_E = fig.add_subplot(gs[0, 0])
        ax_H = fig.add_subplot(gs[0, 1])

        for ax, data, label in zip(
            [ax_E, ax_H], [grid_data_E, grid_data_H], [r"$|E(x,y)|$", r"$|H(x,y)|"]
        ):

            ax.imshow(
                np.sqrt(np.sum(np.abs(data)**2, axis=2)),
                origin="lower",
                extent=[xmin, xmax, ymin, ymax],
                cmap=cmap,
                interpolation="bicubic",
                aspect="auto",
            )

            ax.streamplot(
                grid_x,
                grid_y,
                data.real[:, :, 0],
                data.real[:, :, 1],
                color=color_vectors,
                linewidth=0.5,
            )

            self.plot_polygons(
                fig=fig,
                ax=ax,
                color_polygon=color_polygons,
                color_line=color_integral_lines,
            )

            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)

            ax.set_xlabel("x (um)")
            ax.set_ylabel("y (um)")

            ax.set_title(label)

    def plot_polygons(
        self,
        color_polygon="black",
        color_line="green",
        color_junctions="blue",
        fig=None,
        ax=None,
    ):
        """
        Plots the polygon and line entities in the simulation window.

        Args:
            color_polygon: The color of the polygons in the plot.
            color_line: The color of the line entities in the plot.
            color_junctions: The color of the junction entities in the plot.
            fig: The matplotlib figure to plot on. If `None`, a new figure is created.
            ax: The matplotlib axis to plot on. If `None`, a new axis is created.
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
            elif isinstance(poly, Line):
                ax.plot(*poly.xy, color=color_line)

    def plot_mesh(
        self,
        plot_polygons: bool = True,
    ):
        """
        Plots the mesh of the photonic device.

        Args:
            plot_polygons: Whether to plot the polygons on the mesh. If ``True`` it will call the :meth:`imodulator.RFSimulator.plot_polygons` method.

        Returns:
            fig: The matplotlib figure.
            ax: The matplotlib axis.
        """
        ax = draw_mesh2d(self.mesh)
        ax.set_axis_on()

        fig = ax.get_figure()

        if plot_polygons:
            self.plot_polygons(color_polygon="red", fig=fig, ax=ax)

        return fig, ax

    def get_currents(self, mode: Mode):
        r"""
        This function calculates the currents in every polygon marked as ``calculate_current = True`` via:

        .. math::


            i_0 = \int_{\partial \Omega} \mathbf{n} \times \mathbf{H} \cdot d\mathbf{s}

        Note that every field returned by `FEMWELL`_ is **power normalized by orthogonality relations**. Therefore, if you're using a symmetry plane, you have to take into account that this code will return a value taking into account that HALF the field has power 1. Please adjust the values to your specific case/symmetry. The reccomended approach is to calculate once with the full structure and then find the proper scaling factors.

        Args:
            mode: The mode object containing the electric and magnetic field data to be used for the calculation.

        Returns:
            p0 (complex): The power normalization factor. It will be **very** close to 1, but it is not exactly 1 due to the numerical integration.
            currents (dict): A dictionary containing the calculated currents for each polygon.
            impedances (dict): A dictionary containing the calculated impedances for each polygon.
        """

        @Functional(dtype=np.complex128)
        def current_form(w):
            """
            What this does is it takes the normal vector to the boundary and rotates it 90deg
            Then takes the inner product with the magnetic field in it's complex form
            """
            return inner(np.array([w.n[1], -w.n[0]]), w.H)

        p0 = calculate_scalar_product(
            mode.basis, np.conjugate(mode.E), mode.basis, np.conjugate(mode.H)
        )  # The conjugate is due to femwell internal definition as conj(E) x H

        (ht, ht_basis), (hz, hz_basis) = mode.basis.split(mode.H)

        currents = {}
        impedances = {}

        for photopoly in self.rf_photopolygons:
            if photopoly.calculate_current:
                facet_basis = ht_basis.boundary(
                    facets=self.mesh.boundaries[photopoly.name + "line_current"]
                )
                i0 = current_form.assemble(facet_basis, H=facet_basis.interpolate(ht))
                currents[photopoly.name] = i0
                impedances[photopoly.name] = p0 / np.abs(i0) ** 2

        return p0, currents, impedances

    def get_RLGC(self, mode: Mode, i0: complex = 1.0, v0: complex = 1.0):
        r"""
        This function calculates the RLGC parameters calculated from waveguide circuit theory for a TEM field. 

        The :math:`i_0` and :math:`v_0` values provided should be the values of current and voltage across the signal and ground tracks of the equivalent transmission line. This is not always straightforward, particularly if you're dealing with multiconductor transmission lines. It is reccomended to use this formalism only when you have a good equivalent transmission line of the present mode.

        Please note that every field returned by `FEMWELL`_ is power normalized by orthogonality relations. Therefore, if you're using a symmetry plane, you have to take into account that this code will return a value taking into account that HALF the field has power 1. Please adjust the values to your specific case/symmetry. The reccomended approach is to calculate once with the full structure and then find the proper scaling factors.

        The RLGC parameters are calculated according to :cite:t:`marks_general_1992`:

        .. math::

            C = \frac{1}{|v_0|^2} \left[ \int_S \Re\{\epsilon_{RF}\} |E_t|^2dS - \int_S \Re\{\mu_{RF}\}|H_z|^2 dS \right] \\
            L = \frac{1}{|i_0|^2} \left[ \int_S \Re\{\mu_{RF}\} |H_t|^2dS - \int_S \Re\{\epsilon_{RF}\}|E_z|^2 dS \right] \\
            G = \frac{\omega}{|v_0|^2} \left[ \int_S \Im\{\epsilon_{RF}\} |E_t|^2dS + \int_S \Im\{\mu_{RF}\}|H_z|^2 dS \right] \\
            R = \frac{\omega}{|i_0|^2} \left[ \int_S \Im\{\mu_{RF}\} |H_t|^2dS + \int_S \Im\{\epsilon_{RF}\}|E_z|^2 dS \right] 
        
        Args:
            mode: The mode object containing the electric and magnetic field data to be used for the calculation.
            i0: The current value to use for the calculation.
            v0: The voltage value to use for the calculation.

        Returns:
            R (float): The resistance per unit length in Ohms per micrometer.
            L (float): The inductance per unit length in Henrys per micrometer.
            G (float): The conductance per unit length in Siemens per micrometer.
            C (float): The capacitance per unit length in Farads per micrometer.
        """

        reg = self.photodevice.reg
        frequency = 3e8 / mode.wavelength / 1e3 * reg.GHz
        omega = 2 * np.pi * frequency
        e0 = self.photodevice.e0
        mu0 = self.photodevice.mu0
        c = self.photodevice.c

        @Functional(dtype=np.complex64)
        def C_form(w):

            return (
                1
                / np.abs(w.v0) ** 2
                * (
                    np.real(w.epsilon) * inner(w.et, np.conj(w.et))
                    - np.real(w.mu) * inner(w.hz, np.conj(w.hz))
                )
            )

        @Functional(dtype=np.complex64)
        def L_form(w):

            return (
                1
                / np.abs(w.i0) ** 2
                * (
                    np.real(w.mu) * inner(w.ht, np.conj(w.ht))
                    - np.real(w.epsilon) * inner(w.ez, np.conj(w.ez))
                )
            )

        @Functional(dtype=np.complex64)
        def G_form(w):
            # The minus sign account for the fact that in the paper they define eps = eps_r-1j*eps_i
            # whereas with python we have eps = eps_r+1j*eps_i.
            return (
                -w.omega
                / np.abs(w.v0) ** 2
                * (
                    np.imag(w.epsilon) * inner(w.et, np.conj(w.et))
                    + np.imag(w.mu) * inner(w.hz, np.conj(w.hz))
                )
            )

        @Functional(dtype=np.complex64)
        def R_form(w):
            # The minus sign account for the fact that in the paper they define eps = eps_r-1j*eps_i
            # whereas with python we have eps = eps_r+1j*eps_i.
            return (
                -w.omega
                / np.abs(w.i0) ** 2
                * (
                    np.imag(w.epsilon) * inner(w.ez, np.conj(w.ez))
                    + np.imag(w.mu) * inner(w.ht, np.conj(w.ht))
                )
            )

        (ht, ht_basis), (hz, hz_basis) = mode.basis.split(mode.H)
        (et, et_basis), (ez, ez_basis) = mode.basis.split(mode.E)

        epsilon = mode.epsilon_r

        basis_t = et_basis
        basis_z = ez_basis
        basis_eps = mode.basis_epsilon_r

        # Be careful with the units!!
        # In this case just make sure you adjust the length unit to micrometers
        # everything else can stay as is.

        C = C_form.assemble(
            mode.basis,
            epsilon=basis_eps.interpolate(
                epsilon * e0.to(reg.farad / reg.micrometer).magnitude
            ),
            mu=mu0.to(reg.henry / reg.micrometer).magnitude,
            i0=i0,
            omega=omega.to(reg.second**-1).magnitude,
            v0=v0,
            et=basis_t.interpolate(et),
            ez=basis_z.interpolate(ez),
            ht=basis_t.interpolate(ht),
            hz=basis_z.interpolate(hz),
        )

        L = L_form.assemble(
            mode.basis,
            epsilon=basis_eps.interpolate(
                epsilon * e0.to(reg.farad / reg.micrometer).magnitude
            ),
            mu=mu0.to(reg.henry / reg.micrometer).magnitude,
            i0=i0,
            omega=omega.to(reg.second**-1).magnitude,
            v0=v0,
            et=basis_t.interpolate(et),
            ez=basis_z.interpolate(ez),
            ht=basis_t.interpolate(ht),
            hz=basis_z.interpolate(hz),
        )

        G = G_form.assemble(
            mode.basis,
            epsilon=basis_eps.interpolate(
                epsilon * e0.to(reg.farad / reg.micrometer).magnitude
            ),
            mu=mu0.to(reg.henry / reg.micrometer).magnitude,
            i0=i0,
            omega=omega.to(reg.second**-1).magnitude,
            v0=v0,
            et=basis_t.interpolate(et),
            ez=basis_z.interpolate(ez),
            ht=basis_t.interpolate(ht),
            hz=basis_z.interpolate(hz),
        )

        R = R_form.assemble(
            mode.basis,
            epsilon=basis_eps.interpolate(
                epsilon * e0.to(reg.farad / reg.micrometer).magnitude
            ),
            mu=mu0.to(reg.henry / reg.micrometer).magnitude,
            i0=i0,
            omega=omega.to(reg.second**-1).magnitude,
            v0=v0,
            et=basis_t.interpolate(et),
            ez=basis_z.interpolate(ez),
            ht=basis_t.interpolate(ht),
            hz=basis_z.interpolate(hz),
        )

        # The units returned are:
        # mOhm/um
        # picoHenry/um
        # picoSiemens/um
        # femtofarad/um

        return (
            R.real / 1e-3 * self.reg.milliohm / self.reg.micrometer,
            L.real / 1e-12 * self.reg.picohenry / self.reg.micrometer,
            G.real / 1e-15 * self.reg.femtosiemens / self.reg.micrometer,
            C.real / 1e-15 * self.reg.femtofarad / self.reg.micrometer,
        )

    def get_power_loss_per_polygon(self, mode: Mode, frequency: float):
        r"""
        This function returns the power loss per polygon as described in :cite:t:`pozar_microwave_2012` eq. 1.92. The result is given in W/cm:

        .. math::

            P_l = \int_V \frac{\sigma}{2}|E|^2dv + \frac{\omega}{2}\int_V(\Im\{\epsilon\}|E|^2 + \Im\{\mu\}|H|^2)dv


        .. warning:: This part is not developed yet. 
        
        Args:
            mode: The mode object containing the electric and magnetic field data to be used for the calculation.
            frequency: The frequency at which to compute the power loss. The frequency must be in GHz.
        Returns:
            power_lost_all (dict): A dictionary containing the power lost per polygon
        """

        warnings.warn(
            "The loss per polygon functionality is under development!! Results are not correct."
        )

        from skfem.helpers import dot

        @Functional
        def factor(w):
            return w["epsilon"] * (
                dot(np.conj(w["E"][0]), w["E"][0]) + np.conj(w["E"][1]) * w["E"][1]
            )

        power_lost_all = {}
        for photo_polygon in self.rf_photopolygons:
            # Retrieve all vertice idxs that belong to the polygon
            elements = self.mesh.subdomains[photo_polygon.name]

            basis = mode.basis.with_elements(elements)
            basis_epsilon_r = mode.basis_epsilon_r.with_elements(elements)
            power_lost = (
                (2 * np.pi * frequency * 1e9)
                / 2
                * 8.85e-12
                * 1e-6  # To get F/um
                * factor.assemble(
                    basis,
                    E=basis.interpolate(mode.E),
                    epsilon=basis_epsilon_r.interpolate(mode.epsilon_r),
                )
            ).imag * 1e4  # To get W/cm

            power_lost_all[photo_polygon.name] = (
                power_lost * self.reg.watt / self.reg.centimeter
            )

        return power_lost_all

    def get_S(
        self,
        gamma: Quantity,
        Z: Quantity,
        ZL: float = 50.0,
        ZS: float = 50.0,
        L: float = 1e3,
    ):
        r""" 
        We calculate the S parameters based on transmission line theory :cite:t:`rizzi_microwave_1988`, starting from the ABCD matrix and going to S parameters:

        .. math::

            S_{11} = \frac{AZ_L + B - CZ_L*Z_S - DZ_S}{AZ_L + B + CZ_LZ_S + DZ_S} \\
            S_{22} = \frac{2\sqrt{Z_S Z_L}}{AZ_L + B + CZ_LZ_S + DZ_S}

        Args:
            gamma: The propagation constant.
            Z: The characteristic impedance.
            ZL: The load impedance.
            ZS: The source impedance.
            L: The length of the transmission line.

        Returns:
            S11: The S11 parameter.
            S12: The S12 parameter
        """

        ZL = ZL * self.reg.ohm
        ZS = ZS * self.reg.ohm
        L = L * self.reg.micrometer

        self.S = np.zeros((2, 2), dtype=complex)

        # This is to adhere to the same gamma convention of the RF theory where gamma = alpha + 1j*beta
        # whereas here we deal with k = beta - 1j*alpha
        k = np.abs(gamma.imag) + 1j * np.abs(gamma.real)

        A = np.cosh(k * L)
        B = Z * np.sinh(k * L)
        C = np.sinh(k * L) / Z
        D = np.cosh(k * L)

        S11 = (A * ZL + B - C * ZL * ZS - D * ZS) / (A * ZL + B + C * ZL * ZS + D * ZS)
        S12 = (2 * np.sqrt(ZS * ZL)) / (A * ZL + B + C * ZL * ZS + D * ZS)

        return (
            S11.to(self.reg.dimensionless).magnitude,
            S12.to(self.reg.dimensionless).magnitude,
        )
    
