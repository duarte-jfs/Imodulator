"""
The role of the OpticalSimulatorFEMWELL is to:

- generate the optical mesh;
- Apply a simulation window;
- Generate field visualizations of the optical modes;
- Generate mesh visuals;
- Solve for the optical modes;

"""

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
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator

from skfem.io.meshio import from_meshio
from skfem.visuals.matplotlib import draw_mesh2d
from skfem import Basis, ElementTriP1, ElementVector, ElementDG, Functional
from skfem.helpers import inner

from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import (
    Modes,
    Mode,
    compute_modes,
)

import numpy as np
import copy

from imodulator import PhotonicDevice
from imodulator.ElectroOpticalModel import ElectroOpticalModel
import imodulator.Config as Config
lumapi=Config.config_instance.get_lumapi()

from collections import OrderedDict
import warnings

from imodulator.PhotonicPolygon import (
    SemiconductorPolygon,
    MetalPolygon,
    InsulatorPolygon,
)

PhotonicPolygon = SemiconductorPolygon | MetalPolygon | InsulatorPolygon
Line = LineString | MultiLineString | LinearRing


class OpticalSimulatorMODE:
    #The api import can be cleaner dunno how
    
    # The default paths for windows
    
    """
    Base class for Mode simulation of a PhotonicDevice based on Lumerical MODE
    """

    def __init__(
        self,
        device: PhotonicDevice,
        simulation_window: Polygon | None = None,
        include_metals: bool = False,
        save_sim: bool = False,
    ):

        '''
       Initialize the OpticalSimulatorMODE class for MODE simulation.

       Creates geometry and material entities from the photonic device.
       
       With every function call the state of the simulation after the 
       function applied will be saved to the file quicksave.lms.
       quicksave.lms helps with chenking the state of the simulation and 
       prevents multiple of mode simulations to reduce the system memory. 
       
       Way of working:

       1. __init__(); runs 
       2. create_geometry(): creates polygons
       3. create_materials(): creates the material in lumerical material library for each geometry object and assigns them accordingly. The rest will be run by the user to control the simulation parameters
       4. create_fde(): sets the required simulation parameters & creates the simulation region" 
       5. compute_modes(): runs the simulation
       6. select_modes() #select TE and TM mode index or it will be automatically chosen if the top two modes are the fundamental TE and TM mode
       7. plot_modes() #plots the electric field magnitudes for TE and TM modes
       
        
        Args:
            device: The photonic device to simulate
            simulation_window: The simulation window bounds. If None, uses device bounds
            include_metals: Whether to include metal polygons in simulation. Defaults to True
            save_sim: Whether to save simulation files. Defaults to False
            
        .. todo::  
            Permitivity tensor to input n,k values? Not sure if necesary
        

        '''
        
        self.photonicdevice = device
        self.lumapi=lumapi
        self.simulation_window = simulation_window

        self.optical_photopolygons = copy.deepcopy(self.photonicdevice.photo_polygons)

        # remove metals if not include_metals
        if not include_metals:
            self.optical_photopolygons = [
                poly
                for poly in self.optical_photopolygons
                if not isinstance(poly, MetalPolygon)
            ]
            
        self.polygon_entities = OrderedDict() #?
        self.junction_entities = OrderedDict() #?
        self.resolutions = dict()
        self.modefields = dict() 
        self.save_sim = save_sim 


        for polygon in self.optical_photopolygons:
            self.polygon_entities[polygon.name] = polygon.polygon
        # We now have all the photopolygons cut by the plane. Let us finally add the boundaries

        self.entities = OrderedDict(
            list(self.junction_entities.items())
            + list(self.polygon_entities.items())
        )
        # Transfer the resolutions from the photonic device to the optical_simulator
        for name in self.entities.keys():
            if name in self.photonicdevice.resolutions_optical.keys():
                self.resolutions[name] = self.photonicdevice.resolutions_optical[name]
        #Lumerical MODE API starts here      
        self.mode = self.lumapi.MODE()
        self._create_geometry()
        self._create_materials()
        self.mode.save("quicksave"+".lms")

    def _create_geometry(self): #implement a priority order
        """
        Create geometric primitives in Lumerical MODE from polygon entities.
        
        Creates polygon objects in MODE for all non-background polygon entities,
        setting their vertices in micron units and saving the simulation file.
        """
        # create the primitives (change into polygons)
        print("List of polygons in simulation:")
        for key in self.polygon_entities.keys():
        #     if "metal" in key:
                #     continue
            if key == "background":
                continue
            self.mode.addpoly()
            self.mode.set("name", key)
            self.mode.select(key)
            # Extract the coordinates from the entitity
            x, y = self.polygon_entities[key].exterior.coords.xy
            # Combine the coordinates into a single array
            coords = np.column_stack((x, y))* 1e-6 #in microns
            # Set the vertices
            self.mode.set("vertices", coords)
        
        self.mode.save("quicksave"+".lms")
        
        
        
    def _add_nonmetal_mat_mode(self):
        """
        Add non-metallic materials to Lumerical MODE simulation.
        
        Creates (n,k) materials in MODE for semiconductor and insulator polygons,
        setting their refractive index and imaginary refractive index based on
        the optical material properties.
        """
        mode_mats = {}
        for polygon in self.optical_photopolygons:
            if isinstance(polygon, (SemiconductorPolygon,InsulatorPolygon)):
            
                lum_materialname = polygon.name + "_m"
                if self.mode.materialexists(lum_materialname):
                    self.mode.deletematerial(lum_materialname)
                mode_mats[polygon.name] = self.mode.addmaterial("(n,k) Material")
                self.mode.setmaterial(mode_mats[polygon.name], "name", lum_materialname)
                
                n0 = polygon.optical_material.real
                nimag = polygon.optical_material.imag
               
                self.mode.setmaterial(lum_materialname, "Refractive Index", np.sqrt(polygon.optical_material.real))
                self.mode.setmaterial(lum_materialname, "Imaginary Refractive Index", np.sqrt(polygon.optical_material.imag))  
    
    def _add_metal_mat_mode(self) :
        ##
        print("placehold")

        
    def _create_materials(self):#Overwrite meh order for priority
        """
        Create and assign materials to geometric objects in MODE.
        
        Adds non-metallic materials, assigns random contrasting colors to each material,
        and associates materials with their corresponding geometric objects. Background
        materials are set with reduced alpha for transparency.
        """
        ###Semiconductor
        self._add_nonmetal_mat_mode()
        # Number of colors you want
        N = len(self.optical_photopolygons)

        # Generate N random RGB colors
        colors = np.random.randint(0, 256, size=(N, 3))

        # Ensure contrasting colors
        for i in range(1, N):
            while np.linalg.norm(colors[i] - colors[i-1]) < 50:  # Adjust the threshold as needed
                colors[i] = np.random.randint(0, 256, size=3)

        # Convert the colors to a list
        cl = [tuple(color) for color in colors]
        
        for i, polygon in enumerate(self.optical_photopolygons):
            lum_matname=polygon.name+"_m"
            if polygon.name =='background': #implement in fde
                alpha_mat=0.1
            else:
                alpha_mat=1
            self.mode.setmaterial(lum_matname,"color",np.asarray([[cl[i][0] / 255], [cl[i][1] / 255], [cl[i][1] / 255], [alpha_mat]]))
            
            try:
                self.mode.setnamed(polygon.name,"material",lum_matname)
                print(lum_matname)
            except lumapi.LumApiError: #There might be more material defined then geometry, except error and continue in the loop
                continue
        self.mode.save("quicksave"+".lms")
        
    def mesh(self):
        
        for key in self.polygon_entities.keys():
        #     if "metal" in key:
                #     continue
            if key == "background":
                continue
            self.mode.addmesh()
            meshname="mesh_"+key
            self.mode.set("name", meshname)
            self.mode.select(meshname)
            self.mode.set("based on a structure",1)
            self.mode.set("structure",key)
            self.mode.set("buffer",0)
            # self.mode.set("buffer",self.resolutions[key]["distance"])
            self.mode.set("dx",self.resolutions[key]["resolution"]*1e-6)
            self.mode.set("dy",self.resolutions[key]["resolution"]*1e-6)
            
            
    def create_fde(self,
                   num_modes: int = 4,
                   dx=20e-9,
                   dy=10e-9,
                   background_index=1,
                   bc={"x min bc":"PML",
                       "x max bc":"PML",
                       "y min bc":"PML",
                       "y max bc":"PML",
                       },
                   bounds=None,
    ):
        """
        Create and configure the Finite Difference Eigenmode (FDE) solver in Lumerical MODE.
        
        Sets up a 2D Z-normal FDE solver with specified mesh parameters, boundary conditions,
        and simulation region. The solver will find eigenmode solutions for the waveguide
        cross-section defined by the device geometry.
        
        Args:
            num_modes (int, optional): Number of trial modes to calculate. The solver will
                attempt to find this many eigenmodes. Defaults to 4.
                
            dx (float, optional): Maximum mesh step size in the x-direction in meters.
                Smaller values give more accurate results but increase computation time.
                Defaults to 20e-9 (20 nm).
                
            dy (float, optional): Maximum mesh step size in the y-direction in meters.
                Smaller values give more accurate results but increase computation time.
                Defaults to 10e-9 (10 nm).
                
            background_index (float, optional): Background refractive index used as the
                initial guess for the mode search. Should be close to the expected
                effective index of the fundamental mode. Defaults to 1.
                
            bc (dict, optional): Boundary conditions for each edge of the simulation region.
                Keys should be "x min bc", "x max bc", "y min bc", "y max bc".
                Values can be "PML" (Perfectly Matched Layer), "Metal", or "PMC" 
                (Perfect Magnetic Conductor). Defaults to PML on all boundaries.
                
            bounds (list or tuple, optional): Custom simulation region bounds as 
                [x_min, y_min, x_max, y_max] in meters. If None, automatically sets
                bounds based on the background polygon with 1 μm padding on all sides.
                Defaults to None.
        
        Returns:
            None: Configures the FDE solver in the Lumerical MODE session.
            
        Notes:
            - The solver is configured for 2D Z-normal geometry (cross-section analysis)
            - Mesh is defined by maximum step size ("maximum mesh step" method)
            - Mode search uses "near n" method with the background_index as target
            - The "use max index" option is enabled for better convergence
            - Simulation state is automatically saved as "quicksave.lms"
            
        Example:
            >>> # Basic usage with defaults
            >>> Mode.create_fde()
            
            >>> # Custom mesh and more modes
            >>> Mode.create_fde(num_modes=8, dx=10e-9, dy=5e-9, background_index=3.2)
            
            >>> # Custom boundary conditions
            >>> bc_custom = {"x min bc": "Metal", "x max bc": "Metal", 
            ...              "y min bc": "PML", "y max bc": "PML"}
            >>> Mode.create_fde(bc=bc_custom)
            
            >>> # Custom simulation bounds (in meters)
            >>> bounds = [-10e-6, -5e-6, 10e-6, 5e-6]  # 20x10 μm region
            >>> Mode.create_fde(bounds=bounds)
        """
        
        scaled_bounds = tuple(element * 1e-6 for element in self.polygon_entities["background"].bounds)
        
        self.mode.addfde()
        self.mode.select("FDE")
       
        self.mode.set("solver type", "2D Z normal")
        #expend the simulation region
        if bounds==None:
            self.mode.set("x min", scaled_bounds[0]-1e-6)
            self.mode.set("y min", scaled_bounds[1]-1e-6)
            self.mode.set("x max", scaled_bounds[2]+1e-6)
            self.mode.set("y max", scaled_bounds[3]+1e-6)  
        elif bounds != None:
            self.mode.set("x min", bounds[0])
            self.mode.set("y min", bounds[1])
            self.mode.set("x max", bounds[2])
            self.mode.set("y max", bounds[3])

        #Work on the mesh
        # mesh settings
        self.mode.select("FDE")
        self.mode.set("define x mesh by", "maximum mesh step")
        self.mode.set("define y mesh by", "maximum mesh step")
        self.mode.set("dx", dx)
        self.mode.set("dy", dy)

        #Mode search
        self.mode.select("FDE")
        self.mode.set("search", "near n")
        self.mode.set("use max index", True)
        self.mode.set("number of trial modes",num_modes)
        self.mode.set("index", background_index) #background index
        
        for bc_key in bc.keys():
            self.mode.set(bc_key,bc[bc_key])
        
        self.mode.save("quicksave"+".lms")
        
            
        
    def compute_modes(
        self,
        ):  
        """
        Compute optical modes using the FDE solver.
        
        Args:
            wavelength (float): Operating wavelength in microns. Defaults to 1.55
            voltage_idx (int): Voltage index for electro-optic simulations. Defaults to 0
            num_modes (int): Number of modes to compute. Defaults to 4
            order (int): Mode order. Defaults to 1
            return_modes (bool): Whether to return mode data. Defaults to False
            auto_select (bool): Whether to auto-select TE/TM modes. Defaults to True
            
        Runs the MODE findmodes() solver and optionally saves the simulation file with the device class name used.
        """
        self.mode.findmodes()
        self.mode.save("quicksave"+".lms")
        
        #Save the file if asked
        if self.save_sim:
            self.mode.save(self.photonicdevice.__class__.__name__+".lms")
   
       
    def select_modes(
        self,
        auto_select: bool = True,  # works only if the top two mods are the fundamental modes
        TE_TM_idx: tuple[int,int] = (1,2),             
                     ) -> tuple[int,int]:
        """
        Select and extract TE and TM mode field data from simulation results.
        
        Args:
            auto_select (bool): If True, automatically selects TE/TM based on polarization fraction.
                               If False, uses provided indices. Defaults to True
            TE_TM_idx (tuple[int,int]): Mode indices for TE and TM modes when auto_select=False.
                                       Defaults to (1,2)
                                       
        Returns:
            tuple[int,int]: The selected TE and TM mode indices
            
        Extracts Ex, Ey, Ez, Hx, Hy, Hz field components and coordinate arrays,
        storing them in self.modefields dictionary structure.
        """
        if auto_select: 
            mode1TEfrac=self.mode.getresult("FDE::data::mode" + str(1), "TE polarization fraction")
            mode2TEfrac=self.mode.getresult("FDE::data::mode" + str(2), "TE polarization fraction")
            if mode2TEfrac>mode1TEfrac: # make sure the order is TE then TM
                TE_idx=2
                TM_idx=1
            else:
                TE_idx=1
                TM_idx=2 

        if not auto_select: 
            TE_idx=TE_TM_idx[0]
            TE_idx=TE_TM_idx[1]

        
        print("TE and TM mode indexes;\n")
            
        
        TE_Efield=self.mode.getresult("FDE::data::mode" + str(TE_idx), "E")
        TE_Hfield=self.mode.getresult("FDE::data::mode" + str(TE_idx), "H")
        TM_Efield=self.mode.getresult("FDE::data::mode" + str(TM_idx), "E")
        TM_Hfield=self.mode.getresult("FDE::data::mode" + str(TM_idx), "H")
        
        
        self.modefields = {
        "TE": {
            "Ex": TE_Efield["E"][:, :, 0, 0, 0],
            "Ey": TE_Efield["E"][:, :, 0, 0, 1],
            "Ez": TE_Efield["E"][:, :, 0, 0, 2],
            "Hx": TE_Hfield["H"][:, :, 0, 0, 0],
            "Hy": TE_Hfield["H"][:, :, 0, 0, 1],
            "Hz": TE_Hfield["H"][:, :, 0, 0, 2],
            "x": TE_Efield["x"],
            "y": TE_Efield["y"]
        },
        "TM": {
            "Ex": TM_Efield["E"][:, :, 0, 0, 0],
            "Ey": TM_Efield["E"][:, :, 0, 0, 1],
            "Ez": TM_Efield["E"][:, :, 0, 0, 2],
            "Hx": TM_Hfield["H"][:, :, 0, 0, 0],
            "Hy": TM_Hfield["H"][:, :, 0, 0, 1],
            "Hz": TM_Hfield["H"][:, :, 0, 0, 2],
            "x": TM_Efield["x"],
            "y": TM_Efield["y"]
        },
        }

        # xx, yy = np.meshgrid(Efield_modes["x"][:, 0], Efield_modes["y"][:, 0])
        # print(xx.shape, yy.shape)
        self.mode.save("quicksave"+".lms")
        return TE_TM_idx

    def transfer_results_to_device(
            self,
    ):
        """
        This function will transfer the results from the simulation into the self.device.mode dictionary.

        Returns:
            photonicdevice.mode
        """
        interpolatordict={}
        for pol in ["TE", "TM"]:
            fields = self.modefields[pol]
            interp_fields = {}
            x = np.unique(fields["x"])
            y = np.unique(fields["y"])
            for comp in ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz"]:

                interp_fields[comp] = RegularGridInterpolator(
                    (x, y),
                    fields[comp],
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )
            # interp_fields["x"] = fields["x"]
            # interp_fields["y"] = fields["y"]
            interpolatordict[pol] = interp_fields
        self.photonicdevice.mode=interpolatordict
           
    def _plot_polygons_on_axis(self, ax, color_polygons="white"):
        """
        Helper function to plot polygon outlines on a given matplotlib axis.
        
        Args:
            ax: Matplotlib axes object to plot on
            color_polygons (str): Color for polygon outlines. Defaults to "white"
            
        Converts polygon coordinates to micron scale and plots outlines with
        specified color and transparency.
        """
        for name, poly in self.polygon_entities.items():
            if isinstance(poly, Polygon):
                # Convert coordinates from meters to match the field plot scale
                x_coords, y_coords = poly.exterior.xy
                x_coords = np.array(x_coords) * 1e-6  # Convert to microns to match field data
                y_coords = np.array(y_coords) * 1e-6  # Convert to microns to match field data
                ax.plot(x_coords, y_coords, color=color_polygons, linewidth=1, alpha=1)
            
    def plot_mode(
        self,
        figsize: tuple[float, float] = (10,5),
        color_polygons: str = "white",
        show_polygons: bool = True,
        aspect = None,
        ):
        """
        Plot TE and TM mode intensity distributions.
        
        Args:
            figsize: Figure size (width, height). Defaults to (10,5)
            color_polygons: Color for polygon overlays. Defaults to "white"
            show_polygons: Whether to overlay device geometry. Defaults to True
            aspect: The aspect ratio of the Axes ('equal', 'auto', float, or None). 
                This parameter is particularly relevant for images since it determines 
                whether data pixels are square. Defaults to None.
        
        Returns:
            tuple: (figure, (ax1, ax2)) - matplotlib figure and axes objects
        
            
        Creates side-by-side plots of :math:`|E|^2` intensity for TE and TM modes with
        optional device geometry overlay and colorbars.
        """        
        # Extract TE mode fields
        TE_Ex = self.modefields["TE"]["Ex"].transpose()
        TE_Ey = self.modefields["TE"]["Ey"].transpose()
        TE_Ez = self.modefields["TE"]["Ez"].transpose()
        TE_x = self.modefields["TE"]["x"]
        TE_y = self.modefields["TE"]["y"]
        
        # Extract TM mode fields
        TM_Ex = self.modefields["TM"]["Ex"].transpose()
        TM_Ey = self.modefields["TM"]["Ey"].transpose()
        TM_Ez = self.modefields["TM"]["Ez"].transpose()
        TM_x = self.modefields["TM"]["x"]
        TM_y = self.modefields["TM"]["y"]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # Plot TE mode 
        TE_intensity = np.abs(TE_Ex) ** 2 + np.abs(TE_Ey) ** 2 + np.abs(TE_Ez) ** 2
        im1 = ax1.imshow(
            TE_intensity,
            cmap="jet",
            extent=[
                np.min(TE_x[:, 0]),
                np.max(TE_x[:, 0]),
                np.min(TE_y[:, 0]),
                np.max(TE_y[:, 0]),
            ],
            origin="lower",
            aspect = aspect
        )
        
        # Add polygons to TE plot
        if show_polygons:
            self._plot_polygons_on_axis(ax1, color_polygons)
        
        fig.colorbar(im1, ax=ax1, shrink=0.8)
        
        # Plot TM mode
        TM_intensity = np.abs(TM_Ex) ** 2 + np.abs(TM_Ey) ** 2 + np.abs(TM_Ez) ** 2
        im2 = ax2.imshow(
            TM_intensity,
            cmap="jet",
            extent=[
                np.min(TM_x[:, 0]),
                np.max(TM_x[:, 0]),
                np.min(TM_y[:, 0]),
                np.max(TM_y[:, 0]),
            ],
            origin="lower",
            vmin=0,
            aspect = aspect
        )
        
        # Add polygons to TM plot
        if show_polygons:
            self._plot_polygons_on_axis(ax2, color_polygons)
        
        fig.colorbar(im2, ax=ax2, shrink=0.8)
        
        # Set labels and titles
        ax1.set_ylabel("y (m)")
        ax1.set_xlabel("x (m)")
        ax1.set_title(r"$|E|^2$ | TE")
        
        ax2.set_ylabel("y (m)")
        ax2.set_xlabel("x (m)")
        ax2.set_title(r"$|E|^2$ | TM")
        
        fig.tight_layout()
        
        return fig, (ax1, ax2)   


# class OpticalSimulatorFEMWELL:
#     """
#     Base class for optical simulation of a PhotonicDevice using `FEMWELL`_. The role of the optical simulator is to:


#     * generate the optical mesh;
#     * apply optical symmetry planes;
#     * Generate field visualizations of the optical modes;
#     * Generate mesh visuals;
#     * Solve for the optical modes;
#     * Transfer the optical modes to a new mesh (for EO calculations)

#     """

#     def __init__(
#         self,
#         device: PhotonicDevice,
#         simulation_window: Polygon | None = None,
#         include_metals: bool = True,
#     ):

#         '''
#         Initializes the simulator with the photonic device, the simulation window and the include_metals flag.

#         Args:
#             device: The :class:`PhotonicDevice` object to simulate.
#             simulation_window: The simulation window to use for the optical simulation. If None, the entire device is simulated.
#             include_metals: if `True`, :class:`MetalPolygon`s will be included in the simulation. If `False`, they will be ignored and their optical properties assigned to :math:`\epsilon_{opt} = 1`.

#         '''
#         self.photonicdevice = device
#         self.reg = self.photonicdevice.reg

#         self.e = 1.602176634e-19 * self.reg.coulomb
#         self.e0 = 8.854e-12 * self.reg.farad * self.reg.meter**-1
#         self.c = 3e8 * self.reg.meter * self.reg.second**-1  # m s^-1
#         self.mu0 = (
#             4 * np.pi * 1e-7 * self.reg.henry / self.reg.meter
#         )  # vacuum magnetic permeability

#         self.simulation_window = simulation_window

#         self.optical_photopolygons = copy.deepcopy(self.photonicdevice.photo_polygons)

#         # remove metals if not include_metals
#         if not include_metals:
#             self.optical_photopolygons = [
#                 poly
#                 for poly in self.optical_photopolygons
#                 if not isinstance(poly, MetalPolygon)
#             ]

#         self.line_entities = OrderedDict()
#         self.polygon_entities = OrderedDict()
#         self.junction_entities = OrderedDict()
#         self.resolutions = dict()


#         # THIS NEEDS A REFACTOR!!!

#         if simulation_window is None:
#             for name, poly in self.photonicdevice.junction_entities.items():
#                 self.junction_entities[name] = poly

#         elif simulation_window is not None:
#             #Enforce the simulation plane to be a rectangle
#             if not np.isclose(
#                 simulation_window.minimum_rotated_rectangle.area, simulation_window.area
#             ):
#                 raise ValueError("symmetry plane must be a rectangle")

#             # Cut all photopolygons by the simulation plane
#             idxs_to_pop = []
#             for i, poly in enumerate(self.optical_photopolygons):
#                 if poly.polygon.intersects(
#                     simulation_window
#                 ) and not simulation_window.contains(poly.polygon):
#                     poly_tmp = clip_by_rect(poly.polygon, *simulation_window.bounds)

#                     if poly_tmp.is_empty:
#                         idxs_to_pop.append(i)
#                     else:
#                         self.optical_photopolygons[i].polygon = poly_tmp

#                 elif not poly.polygon.intersects(simulation_window):
#                     idxs_to_pop.append(i)

#             for index in sorted(idxs_to_pop, reverse=True):
#                 del self.optical_photopolygons[index]

#             # Select all junctions cutted by the symmetry plane
#             for name, poly in self.photonicdevice.junction_entities.items():

#                 if poly.intersects(simulation_window) and not simulation_window.contains(
#                     poly
#                 ):
#                     poly_tmp = clip_by_rect(polygon, *simulation_window.bounds)

#                     if not poly_tmp.is_empty:
#                         self.junction_entities[name] = poly_tmp

#                 elif simulation_window.contains(poly):
#                     self.junction_entities[name] = poly

#         for polygon in self.optical_photopolygons:
#             self.polygon_entities[polygon.name] = polygon.polygon
#         # We now have all the photopolygons cut by the plane. Let us finally add the boundaries

#         surf_bounds = self.polygon_entities["background"].bounds

#         left = LineString(
#             [(surf_bounds[0], surf_bounds[1]), (surf_bounds[0], surf_bounds[3])]
#         )

#         bottom = LineString(
#             [(surf_bounds[0], surf_bounds[1]), (surf_bounds[2], surf_bounds[1])]
#         )

#         right = LineString(
#             [(surf_bounds[2], surf_bounds[1]), (surf_bounds[2], surf_bounds[3])]
#         )

#         top = LineString(
#             [(surf_bounds[0], surf_bounds[3]), (surf_bounds[2], surf_bounds[3])]
#         )

#         self.line_entities["left"] = left
#         self.line_entities["bottom"] = bottom
#         self.line_entities["right"] = right
#         self.line_entities["top"] = top

#         self.entities = OrderedDict(
#             list(self.line_entities.items())
#             + list(self.junction_entities.items())
#             + list(self.polygon_entities.items())
#         )



#         # Transfer the resolutions from the photonic device to the optical_simulator
#         for name in self.entities.keys():
#             if name in self.photonicdevice.resolutions_optical.keys():
#                 self.resolutions[name] = self.photonicdevice.resolutions_optical[name]

#     def get_epsilon_optical(
#         self,
#         use_charge_transport_data: bool = False,
#     ):
#         """
#         This function will return the :math:`\epsilon_{opt}` tensor with the signature ``self.epsilon_optical[vertice_idx, voltage_idx]`` where ``vertice_idx`` is the index of the vertice in the mesh and `voltage_idx` is the index of the voltage in the bias points. In case no bias dependent data is available (i.e. no charge transport data is available) the ``voltage_idx`` must be 0 as the second axis of the array will have size 1.

#         Args:
#             use_charge_transport_data: Whether to use the charge transport data to compute the permittivity tensor. Doing so will yield a :math:`\epsilon_{opt}(x,y,V)`. Make sure your mesh is appropriate for it. In order to use this option you must define a :class:`ElectroOpticalModel` for each of the polygons you want to compute for. Those polygons that do not have an :class:`ElectroOpticalModel` will have their permittivity set to ``optical_material**2``. Default is False.
#             voltage_idx: The index of the voltage to use for the permittivity tensor.

#         """

#         if self.mesh is None:
#             raise RuntimeError(
#                 "Cannot attribute an optical epsilon. The mesh of the optical simulator is not yet generated"
#             )

#         if use_charge_transport_data == False:
#             self.epsilon_optical = np.zeros(
#                 (3, 3, self.mesh.nvertices, 1), dtype=np.complex128
#             )

#             self.epsilon_optical[0, 0, :, :] = 1
#             self.epsilon_optical[1, 1, :, :] = 1
#             self.epsilon_optical[2, 2, :, :] = 1

#         else:
#             for photo_poly in self.optical_photopolygons:
#                 if isinstance(photo_poly, SemiconductorPolygon):
#                     N_bias_points = len(photo_poly.Ec)
#                     break
#             self.epsilon_optical = np.zeros(
#                 (3, 3, self.mesh.nvertices, N_bias_points), dtype=np.complex128
#             )

#             self.epsilon_optical[0, 0, :, :] = 1
#             self.epsilon_optical[1, 1, :, :] = 1
#             self.epsilon_optical[2, 2, :, :] = 1

        
#         # the self.photo_polygons is created so that idx 0 has higher priority over idx 1
#         # Here, however, if we loop through the photo_polygons from idx 0 to idx N
#         # the hierarchy on the boundaries will be inverted. That is,
#         #lower lying polygons in hierarchy will dominate the boundaries. Therefore, we need to
#         # loop over the inverted list of photo_polygons

#         for photo_polygon in self.optical_photopolygons[::-1]:
#             if (
#                 hasattr(photo_polygon, "electro_optic_module")
#                 and photo_polygon.electro_optic_module is None
#             ) or (
#                 hasattr(photo_polygon, "electro_optic_module")
#                 and use_charge_transport_data is False
#             ):
#                 # Retrieve all vertice idxs that belong to the polygon
#                 elements_idxs = self.mesh.subdomains[photo_polygon.name]
#                 triangles = self.mesh.t[:, elements_idxs]
#                 vertices_idxs = np.unique(triangles.flatten())

#                 # if it has the optical_material as a number then it will put that number on. If not then it will pass it
#                 if isinstance(photo_polygon.optical_material, float | int | complex):
#                     for i in range(3):
#                         self.epsilon_optical[i, i, vertices_idxs] = (
#                             photo_polygon.optical_material
#                         )

#                 elif isinstance(photo_polygon.optical_material, str):
#                     warnings.warn(
#                         f"The optical_material of {photo_polygon.name} is a str. Make sure you insert the material in your mode solver.",
#                         stacklevel=2,
#                     )

#             elif (
#                 hasattr(photo_polygon, "electro_optic_module")
#                 and photo_polygon.electro_optic_module is ElectroOpticalModel
#                 and use_charge_transport_data is True
#             ):

#                 elements_idxs = self.mesh.subdomains[photo_polygon.name]
#                 triangles = self.mesh.t[:, elements_idxs]
#                 vertices_idxs = np.unique(triangles.flatten())

#                 # Loop through all the voltages

#                 for voltage_idx in range(N_bias_points):
#                     # Here we generate the mun, mup, Ec, Ev, Efp, Efn, P, N, Efield for each polygon
#                     x = self.mesh.p[0, vertices_idxs]
#                     y = self.mesh.p[1, vertices_idxs]

#                     mun = photo_polygon.mun[voltage_idx](x,y)
#                     mup = photo_polygon.mup[voltage_idx](x,y)
#                     Ec = photo_polygon.Ec[voltage_idx](x,y)
#                     Ev = photo_polygon.Ev[voltage_idx](x,y)
#                     Efp = photo_polygon.Efp[voltage_idx](x,y)
#                     Efn = photo_polygon.Efn[voltage_idx](x,y)
#                     P = photo_polygon.P[voltage_idx](x,y)
#                     N = photo_polygon.N[voltage_idx](x,y)
#                     Efield = photo_polygon.Efield[voltage_idx](x,y)
                    
#                     self.epsilon_optical[:, :, vertices_idxs, voltage_idx] = (
#                         photo_polygon.electro_optic_module(
#                             mup=mup,
#                             mun=mun,
#                             Ec=Ec,
#                             Ev=Ev,
#                             Efp=Efp,
#                             Efn=Efn,
#                             P=P,
#                             N=N,
#                             Efield=Efield,
#                             reg=self.photonicdevice.reg,
#                             **photo_polygon.electro_optic_module_kwargs,
#                         ).get_perm()
#                     )




#     def make_mesh(
#         self,
#         default_resolution_min: float = 1e-12,
#         default_resolution_max: float = 100,
#         filename: str = None,
#         gmsh_algorithm: int = 5,
#         global_quad: bool = False,
#         verbose: bool = False,
#         mesh_scaling_factor: float = 1.0,
#     ):
#         """
#         Returns a `gmsh <https://gmsh.info/>`_ with the geometric information of the photonic device. The mesh generation is currently handled by `FEMWELL`_, but it will later be deprecated in favour of `MESHWELL`_.

#         Args:
#             default_resolution_min: The minimum resolution to use for the mesh generation if no resolution is specified for a given entity.
#             default_resolution_max: The maximum resolution to use for the mesh generation if no resolution is specified for a given entity.
#             filename: The filename to save the mesh to. If `None`, the mesh is not saved.
#             gmsh_algorithm: The algorithm to use for the mesh generation.
#             global_quad: Whether to use global quadrature for the mesh generation.
#             verbose: Whether to print verbose output during the mesh generation.
#             mesh_scaling_factor: The scaling factor to apply to the mesh.

#         """
#         self.mesh = from_meshio(
#             mesh_from_OrderedDict(
#                 self.entities,
#                 self.resolutions,
#                 default_resolution_max=default_resolution_max,
#                 default_resolution_min=default_resolution_min,
#                 filename=filename,
#                 gmsh_algorithm=gmsh_algorithm,
#                 global_quad=global_quad,
#                 verbose=verbose,
#                 mesh_scaling_factor=mesh_scaling_factor,
#             )
#         )
#         # Choosing each element as ElementTriP1 is crucial as the dofs that assumes are the mesh vertices
#         self.basis = Basis(self.mesh, ElementTriP1(), intorder=4)

#     def compute_modes(
#         self,
#         wavelength: float = 1.55,
#         voltage_idx: int = 0,
#         num_modes: int = 1,
#         order: int = 1,
#         metallic_boundaries: list | str | bool = False,
#         n_guess: float = 4.0,
#         return_modes: bool = False,
#     ) -> Modes:
#         """
#         Compute the electromagnetic RF modes of a photonic device at a given wavelength. The modes are computed via `FEMWELL`_.

#         Args:
#             wavelength: The wavelength at which to compute the modes. The wavelength must be in um.
#             voltage_idx: The index of the voltage to use for the permittivity tensor.
#             num_modes: The number of modes to compute.
#             order: Order of the basis functions to use in the EM solver.
#             metallic_boundaries: The boundaries to treat as metallic. If `False`, no boundaries are treated as metallic. If `True`, all boundaries are treated as metallic. If a list of strings, the boundaries with the given names are treated as metallic. At the moment, the simulation window is treated as a square, therefore, the metallic boundaries can be ``left``, ``right``, ``top`` and ``bottom`` boundaries.
#             n_guess: Initial guess for the effective index.
#             return_modes: Whether to return the computed modes. 
#             use_charge_transport_data: Whether to use the charge transport data to compute the permittivity tensor. Doing so will yield a :math:`\sigma(x,y,V)`. Make sure your mesh is appropriate for it.

#         Returns:
#             The computed modes if `return_modes` is `True`.
#         """

#         self.get_epsilon_optical()

#         modes = compute_modes(
#             self.basis,
#             self.epsilon_optical[0,0,:,voltage_idx],
#             wavelength,
#             mu_r=1,
#             num_modes=num_modes,
#             order=order,
#             metallic_boundaries=metallic_boundaries,
#             n_guess=n_guess,
#         )

#         self.modes = modes.sorted(lambda mode: mode.n_eff.real)[::-1]

#         if return_modes:
#             return self.modes

#     def plot_mode(
#         self,
#         mode: Mode,
#         figsize: tuple[float, float] = (10,3),
#         color_polygons: str = "black"
#     ):
#         """
#         Plots the electric and magnetic fields components of a given mode. The plotting is handled by `sk-fem`_.

#         Args:
#             mode: The mode object containing the electric and magnetic field data to be plotted.
            
#             figsize: Figure size as (width, height). 
            
#             color_polygons: Color of the polygons in the plot. 

#         Returns
#         -------
#         None
#             This method does not return any value. It generates a plot.


#         """
        
#         fig = plt.figure(figsize = figsize)
#         gs = GridSpec(1,3)
#         ax1 = fig.add_subplot(gs[0,0])
#         ax2 = fig.add_subplot(gs[0,1])
#         ax3 = fig.add_subplot(gs[0,2])

#         self.plot_polygons(fig = fig, ax = ax1, color_polygon=color_polygons)
#         self.plot_polygons(fig = fig, ax = ax2, color_polygon=color_polygons)
#         self.plot_polygons(fig = fig, ax = ax3, color_polygon=color_polygons)


#         mode.plot_component('E', 'x', colorbar = True, ax = ax1)
#         mode.plot_component('E', 'y', colorbar = True, ax = ax2)
#         mode.plot_component('E', 'z', colorbar = True, ax = ax3)

#         fig.tight_layout()
    
#     def plot_polygons(
#         self,
#         color_polygon="black",
#         color_line="green",
#         color_junctions="blue",
#         fig=None,
#         ax=None,
#     ):
#         """
#         Plots the polygons of the :class:`PhotonicDevice` object.

#         Args:
#             color_polygon: The color to use for the polygons.
#             color_line: The color to use for the lines.
#             color_junctions: The color to use for the junctions.
#             fig: The figure to plot on. If ``None``, a new figure is created.
#             ax: The axis to plot on. If ``None``, a new axis is created.

#         Returns:
#             fig, ax: The figure and axis objects.
#         """
#         if fig is None and ax is None:
#             fig = plt.figure()
#             ax = fig.add_subplot(111)

#         for name, poly in self.entities.items():
#             if isinstance(poly, Polygon):
#                 ax.plot(
#                     *poly.exterior.xy,
#                     color=color_polygon if "junction" not in name else color_junctions,
#                 )
#             elif isinstance(poly, Line):
#                 ax.plot(*poly.xy, color=color_line)

#     def plot_mesh(
#         self,
#         plot_polygons: bool = True,
#     ):
#         """
#         Plots the mesh on a new figure

#         Args:
#             plot_polygons: Whether to plot the polygons over the mesh or not.

#         Returns:
#             fig, ax: The figure and axis objects.
#         """
#         ax = draw_mesh2d(self.mesh)
#         ax.set_axis_on()

#         fig = ax.get_figure()

#         if plot_polygons:
#             self.plot_polygons(color_polygon="red", fig=fig, ax=ax)

#         return fig, ax
    

#     def transfer_mode_to_new_mesh(
#             self,
#             mode_idx: int,
#             new_mesh: np.ndarray,
#     ):
#         """
#         Transfers the optical mode to a new mesh. The created objects will be used by the :class:`ElectroOpticalSimulator` to compute the electro-optical response of the device. It creates objects of the shape ``(3, new_mesh.p.shape[1])`` where the first axis is the component of the field and the second axis is the vertice index.

#         Args:
#             mode_idx: The index of the mode to transfer.
#             new_mesh: The new mesh to transfer the mode to. This object is expected to be the ``Mesh.p`` object from `sk-fem`_.

#         Returns:
#             Eopt, Hopt: The electric and magnetic field components of the mode in the new mesh.
#         """
#         mode = self.modes[mode_idx]

#         old_mesh = self.mesh

#         Eopt = np.zeros((3,new_mesh.p.shape[1]), dtype = np.complex128)
#         Hopt = np.zeros((3,new_mesh.p.shape[1]), dtype = np.complex128)

#         idx_in = np.where((new_mesh.p[0] > old_mesh.p[0].min()) & 
#                         (new_mesh.p[0] < old_mesh.p[0].max()) & 
#                         (new_mesh.p[1] > old_mesh.p[1].min()) &
#                         (new_mesh.p[1] < old_mesh.p[1].max()))[0]
        
#         basis = mode.basis
#         basis_fix = basis.with_element(ElementVector(ElementTriP1()))

#         (et, et_basis), (ez, ez_basis) = mode.basis.split(mode.E)
#         (ex, ex_basis), (ey, ey_basis) = basis_fix.split(basis_fix.project(et_basis.interpolate(et)))

#         (ht, ht_basis), (hz, hz_basis) = mode.basis.split(mode.H)
#         (hx, hx_basis), (hy, hy_basis) = basis_fix.split(basis_fix.project(ht_basis.interpolate(ht)))

#         Eopt[0, idx_in] = ex_basis.interpolator(ex)(new_mesh.p[:,idx_in])
#         Eopt[1, idx_in] = ey_basis.interpolator(ey)(new_mesh.p[:,idx_in])
#         Eopt[2, idx_in] = ez_basis.interpolator(ez)(new_mesh.p[:,idx_in])

#         Hopt[0, idx_in] = hx_basis.interpolator(hx)(new_mesh.p[:,idx_in])
#         Hopt[1, idx_in] = hy_basis.interpolator(hy)(new_mesh.p[:,idx_in])
#         Hopt[2, idx_in] = hz_basis.interpolator(hz)(new_mesh.p[:,idx_in])


#         return Eopt, Hopt