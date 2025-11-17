from __future__ import annotations

import os
import copy
from collections import OrderedDict

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.cm as cm

from shapely.geometry import (
    Polygon,
    LineString,
    MultiLineString,
    LinearRing,
)
import shapely

from scipy.interpolate import RegularGridInterpolator

from imodulator import PhotonicDevice
from imodulator.PhotonicPolygon import (
    SemiconductorPolygon,
    MetalPolygon,
    InsulatorPolygon,
)

####### SOLCORE imports ##########
from solcore.material_system.material_system import BaseMaterial
from solcore.material_data.mobility import (
    calculate_InAlAs,
    calculate_InGaAs,
    calculate_InGaAsP,
    calculate_InGaP,
    calculate_AlGaAs,
    mobility_low_field,
)
import json
import inspect

import solcore
from solcore import config
from solcore.parameter_system import ParameterSystem
from solcore.poisson_drift_diffusion.DeviceStructure import DefaultProperties
from configparser import ConfigParser

from solcore.solar_cell import Junction, SolarCell, Layer
from solcore.state import State
from solcore.solar_cell_solver import solar_cell_solver

import gmsh

PhotonicPolygon = SemiconductorPolygon | MetalPolygon | InsulatorPolygon
Line = LineString | MultiLineString | LinearRing

# from imodulator.ElectroOpticalModel import InGaAsPElectroOpticalModel
##Configured imports
from imodulator.Config import config_instance
# Get access to imported modules
nn = config_instance.get_nextnanopy()


#References
# https://www.nextnano.com/documentation/tools/nextnanopy/index.html

def get_normalized_vector(line: LineString):
    start = np.array(line.coords[0])
    end = np.array(line.coords[-1])
    vector = end - start
    norm = np.linalg.norm(vector)
    if norm == 0:
        return np.zeros_like(vector)
    return vector / norm

class ChargeSimulatorNN:
    #The api import can be cleaner dunno how
    
    # The default paths for windows
    
    """
        Initialize the ChargeSimulatorNN.

        Args:
            device (PhotonicDevice): The photonic device to simulate.
            simulation_line (LineString): The line along which to perform 1D simulation.
            inputfile_name (str, optional): Name for the nextnano input file. Defaults to "quicksave".
            output_directory (str, optional): Directory for simulation output. Defaults to config value.
            temperature (float, optional): Simulation temperature in Kelvin. Defaults to 300.0.
            bias_start_stop_step (list, optional): Voltage sweep [start, stop, step]. Defaults to [0,1,1].
    
        Way of working:

        1. The PhotonicDevice should be passed
        2. The simualtion line should be defined as shapely.geometry.LineString. 
        
        .. warning::
            Note that this linestring will define the field direction that will later be used for the electro optic calculations, meaning that the order of the points also matters!!
        
        3. The boundaries of the simulation line will act as contacts of 10nm where the starting point will be contact1...
        and the end of the line will be contact2 
        
        4. The defined voltage will be applied thru contact1
        
        5. Check your simulation line in the geometry via

            >>> self.plot_with_simulation_line()
        
        6. In order to run the NNInfile: #the commands are from the nextnanopy

            >>> self.NNinputf.execute(show_log=False,convergenceCheck=True,convergence_check_mode="continue")
        
        7. In order to load already completed results

            >>> self.load_output_data(folderpath=r"AbsolutePath") #right click on the folder with the self.inputfile_name and copy path

        8. In case you will use EO and RF module, move results to photonic device as interpolators to be used by other modules
    """

    def __init__(
        self,
        device: PhotonicDevice, 
        simulation_line: LineString,
        inputfile_name: str ="quicksave", 
        output_directory:str=nn.config.config['nextnano++']['outputdirectory'],
        temperature: float = 300.0,  # Add temperature parameter
        bias_start_stop_step: list = [0,1,1], #contact1 is the bias electrode decide - or + accordingly
        # save_sim: bool = False,
    ):
        
        """
        Initialize the ChargeSimulatorNN.

        Args:
            device: PhotonicDevice instance containing the device geometry and materials.
            simulation_line: LineString defining the simulation line along which to perform 1D simulation.
            inputfile_name: Name for the nextnano input file. Defaults to "quicksave".
            output_directory: Directory for simulation output. Defaults to config value.
            temperature: Simulation temperature in Kelvin. Defaults to 300.0.
            bias_start_stop_step: Voltage sweep [start, stop, step]. Defaults to [0,1,1].
            
        """
             
        self.temperature = temperature
        self.inputfile_name = inputfile_name 
        self.output_directory = output_directory
        self.photonicdevice = device
        self.bias_start_stop_step=bias_start_stop_step

        self.optical_photopolygons = copy.deepcopy(self.photonicdevice.photo_polygons)

        self.polygon_entities = OrderedDict()

        for polygon in self.optical_photopolygons:
            self.polygon_entities[polygon.name] = polygon.polygon
        
        self._select_line(simulation_line=simulation_line)
        self._create_in_file()

        self.sim_vector_norm = get_normalized_vector(simulation_line)
        self.sim_vector_norm = np.array([self.sim_vector_norm[0], self.sim_vector_norm[1], 0])

        # self.input_file = nn.InputFile(inputfile_name+".in")
        # self.input_file.config = nn.config

    def _select_line(self,
                simulation_line: LineString, #to select the region
                ): 
        """
        Find intersections between the simulation line and device polygons.

        Args:
            simulation_line (LineString): The line along which to perform 1D simulation.

        Returns:
            None. Stores intersection segments in self.line_segments.
        """
        line_segments = OrderedDict()

        names_to_print = []
        for polygon_name, polygon_geom in self.polygon_entities.items():
            # Find intersection between simulation line and polygon
            
            intersection = simulation_line.intersection(polygon_geom)
            if (polygon_name == "substrate") or (polygon_name == "background"):
                continue
            # Handle different intersection types
            if intersection.is_empty:
                continue

            # If intersection is a single LineString
            if intersection.geom_type == 'LineString':
                line_segments[f'{polygon_name}'] = intersection

                ## Loop over the photopolygons of the PhotonicDevice to find the proper PhotonicPolygon that will allow for charge transport data to be loaded on:
                for poly in self.photonicdevice.photo_polygons:
                    if poly.name == polygon_name:
                        poly.has_charge_transport_data = True

                names_to_print.append(polygon_name)

        print('Charge transport will take place with:')
        print(*names_to_print, sep='\n')

        # Store the line segments for later use
        self.line_segments = line_segments
        self.simulation_line = simulation_line

    def _create_in_file(self):
        
        """
        Create and write the nextnano input file from PhotonicDevice data.

        Returns:
            None. Writes file to disk and stores input file object.
        """
        
        output_path = os.path.join(self.output_directory, f"{self.inputfile_name}.in") 
        # Ensure output directory exists
        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)    
            
        # Build the complete file content
        content_sections = [
            self._create_global_section(),
            self._create_grid_section(),
            self._create_structure_section(),
            self._create_impurities_section(),
            self._create_classical_section(),
            self._create_poisson_section(),
            self._create_currents_section(),
            self._create_contacts_section(),
            self._create_run_section()
        ]
            # Join all sections
        self.complete_content = "\n".join(content_sections)
            # Write to file
        with open(output_path, 'w') as f:
            f.write(self.complete_content)
            
        print(f"Input file created: {output_path}")
        self.NNinputf=nn.InputFile(output_path)
        self.NNinputf.config=nn.config#makes sure you use the config
        
    def _create_global_section(self):
        """Create the global section of the nextnano input file"""
        output = f"""
        global{{
        $sweepcheck = 1

        simulate1D{{}}

        temperature = {self.temperature} # Kelvin
        substrate{{ name = "InP" }}
        crystal_zb{{
            x_hkl = [1, 0, 0]
            y_hkl = [0, 1, 0]
        }}
        }}
        
        """
        return output

    def _create_grid_section(self):
        """Create the grid section based on line segments"""
        line_definitions=[]
        cummulative_pos=0
        total_items=len(self.line_segments.items())
        for i, (segment_name, line_segment) in enumerate(self.line_segments.items()):
            spacing = self.photonicdevice.resolutions_charge[segment_name]["resolution"]*1e3 #convert from um to nm
            if i == 0: #first contact + initial position 
                #start of the contact
                line_definitions.append(f"\tline{{pos = {round(-10+line_segment.xy[1][0]*10**3,2)} spacing = 2}}")
                # layer1
                line_definitions.append(f"\t\t\tline{{pos = {round(line_segment.xy[1][0]*10**3,2)} spacing = {spacing}}}")
                # cummulative_pos=+line_segment.length*10**3
                cummulative_pos=line_segment.xy[1][1]*10**3
       
            # elif i == 1: 
            #     line_definitions.append(f"\t\tline{{pos = {round(cummulative_pos,2)} spacing = {spacing}}}")
       
            elif 1 <= i < total_items-1:  #non contact lines
                line_definitions.append(f"\t\t\tline{{pos = {round(cummulative_pos,2)} spacing = {spacing}}}")
                cummulative_pos+=line_segment.length*10**3
            
            elif i == total_items-1: #last contact
                line_definitions.append(f"\t\t\tline{{pos = {round(cummulative_pos,2)} spacing = {spacing}}}")
                cummulative_pos+=line_segment.length*10**3
                line_definitions.append(f"\t\t\tline{{pos = {round(cummulative_pos,2)} spacing = 2}}")
                cummulative_pos+=10
                line_definitions.append(f"\t\t\tline{{pos = {round(cummulative_pos,2)} spacing = 2}}")
            

        joined_line_definitions="\n".join(line_definitions)
        output = f"""
        grid{{
            xgrid{{
        {joined_line_definitions}   
            }}
                  
        }}
        """

        return output#"\n".join(output)
    
    def _create_structure_section(self):
        """Create the structure section based on PhotonicDevice polygons"""
        total_items=len(self.line_segments.items())
        region_definitions=[]
        cummulative_pos=0
        self.contact_thickness=10
        for i, (segment_name, line_segment) in enumerate(self.line_segments.items()):
            for polygon_idx, polygon in enumerate(self.optical_photopolygons):
                if segment_name==polygon.name:
                    segment_charge_kwargs=polygon.charge_transport_simulator_kwargs
                    if segment_charge_kwargs["material_definition"] == None:
                        material_def = "Ga(x)In(1-x)As(y)P(1-y)"
                    else:
                        material_def = segment_charge_kwargs["material_definition"]
                    if i == 0: #first contact + initial position 
                        cummulative_pos=-self.contact_thickness+line_segment.xy[1][0]*10**3
                        line = f"line{{x = [{cummulative_pos:.2f},{cummulative_pos+self.contact_thickness:.2f}] }}"

                        region_definitions.append(f"""
            region{{
                {line}
                {f"contact{{name = contact1}}"}
            }}
            """)

                    #non contact lines
                    line = f"line{{x = [{line_segment.xy[1][0]*10**3:.2f},{line_segment.xy[1][1]*10**3:.2f}]}}"
                    region_definitions.append(f"""
            region{{
                {line}
                quaternary_constant{{
                    name = "{material_def}"    
                    alloy_x = {segment_charge_kwargs["alloy_x"]:.2f}
                    alloy_y = {segment_charge_kwargs["alloy_y"]:.2f}
                }}
                doping{{
                    constant{{
                        name= "{segment_charge_kwargs["doping_type"]}-type"
                        conc= {segment_charge_kwargs["doping_conc"]:.2e}    
                    }}
                }}        
            }}
            """)

                    if i == total_items - 1: #last contact
                        line = f"line{{x = [{line_segment.xy[1][1]*10**3:.2f},{self.contact_thickness+line_segment.xy[1][1]*10**3:.2f}] }}"

                        region_definitions.append(f"""
            region{{
                {line}
                {f"contact{{name = contact2}}"}
            }}
            """)
                    
                    break
                else:
                    continue

        joined_region_definitions="".join(region_definitions)
        output = f"""
        structure{{
            output_region_index{{ }}
            output_material_index{{ }}
            output_user_index{{ }}
            output_contact_index{{ }}
            output_alloy_composition{{ }}
            output_impurities{{ }}
            
            region{{
                everywhere{{}}
                binary{{name = 'Air'}}
            }}
            
        {joined_region_definitions}
                  
        }}
        """

        return output#"\n".join(output)
    
    def _create_impurities_section(self):
        """Create the impurities section"""
        return f"""
        impurities{{
            donor {{ name = "n-type" energy = -1000 degeneracy = 2 }}
            acceptor {{ name = "p-type" energy = -1000 degeneracy = 4 }}
        }}"""
    
    def _create_classical_section(self):
        """Create the classical section"""
        return """
        classical{
            Gamma{}
            HH{}
            # LH{}
            # SO{}
            # X{}
            # L{}

            output_bandedges{ averaged = yes}
            output_carrier_densities{}
        }"""
    
    def _create_poisson_section(self):
        """Create the poisson section"""
        return """
        poisson{
            charge_neutral{}
            output_electric_field{}
        }"""
    
    def _create_currents_section(self):
        """Create the currents section"""
        return """
        currents{
            output_mobilities{}
            recombination_model{} #required by the runner
        }"""
    
    def _create_contacts_section(self):
        """Create the contacts section with voltage sweep parameters"""
        return f"""
        contacts{{
            ohmic{{ name = "contact1" bias = [{self.bias_start_stop_step[0]}, {self.bias_start_stop_step[1]}] steps = {self.bias_start_stop_step[2]}}}
            ohmic{{ name = "contact2" bias = 0.0 }}
        }}"""
    
    def _create_run_section(self):
        """Create the run section"""
        return """
        run{
            current_poisson{
                iterations = 10000
                output_log = yes
            }
        }"""
        
    def load_output_data(self,folderpath=None):
        """
        Load and store simulation output data from nextnano results.

        Args:
            folderpath (str, optional): Path to results folder. Uses default if None.

        Reads bias-dependent data from output files indicated in the InFile
        and stores as 2D arrays indexed by bias point and position along the simulation line.

        Creates instance variables:
            grid: Spatial grid positions [nm]
            Ec: Conduction band edge energy [eV] 
            Ev: Valence band edge energy [eV]
            Efn: Electron quasi-Fermi level [eV]
            Efp: Hole quasi-Fermi level [eV] 
            density_electron: Electron density [cm⁻³]
            density_hole: Hole density [cm⁻³]
            electric_field: Electric field [V/cm]
            mobility_electron: Electron mobility [cm²/Vs]
            mobility_hole: Hole mobility [cm²/Vs]

        All arrays have shape (n_bias_points, n_grid_points).
        """
        if folderpath == None:
            nndata=nn.DataFolder(self.NNinputf.folder_output)
        else:
            nndata=nn.DataFolder(folderpath)
        #file locations to be processed
        f_iv = [f for f in nndata.files if 'IV_characteristics.dat' in f][0]
        self.V= pd.read_csv(f_iv,delim_whitespace=True).iloc[:,0]
        f_grid = [f for f in nndata.files if 'grid_x.dat' in f][0]
        self.grid = pd.read_csv(f_grid,delim_whitespace=True)["Position[nm]"].values.tolist()

        self.Ec = np.zeros(shape=(len(self.V),len(self.grid)))
        self.Ev = np.zeros(shape=(len(self.V),len(self.grid)))
        self.Efn = np.zeros(shape=(len(self.V),len(self.grid)))
        self.Efp = np.zeros(shape=(len(self.V),len(self.grid)))
        self.N = np.zeros(shape=(len(self.V),len(self.grid)))
        self.P = np.zeros(shape=(len(self.V),len(self.grid)))
        self.Efield = np.zeros(shape=(len(self.V),len(self.grid)))
        self.mun = np.zeros(shape=(len(self.V),len(self.grid)))
        self.mup = np.zeros(shape=(len(self.V),len(self.grid)))
        #Loops over the bias000x folders
        for i, v in enumerate(self.V):
            #first get the file locations for each data needed
            nnfiles=nndata.folders[i].files
            # Read each .dat file into DataFrames
            self.Ec[i] = pd.read_csv([f for f in nnfiles if 'bandedges.dat' in f][0], delim_whitespace=True)["Gamma_[eV]"]
            self.Ev[i] = pd.read_csv([f for f in nnfiles if 'bandedges.dat' in f][0], delim_whitespace=True)["HH_[eV]"]
            self.Efn[i] = pd.read_csv([f for f in nnfiles if 'bandedges.dat' in f][0], delim_whitespace=True)["electron_Fermi_level_[eV]"]
            self.Efp[i] = pd.read_csv([f for f in nnfiles if 'bandedges.dat' in f][0], delim_whitespace=True)["hole_Fermi_level_[eV]"]
            self.N[i] = pd.read_csv([f for f in nnfiles if 'density_electron.dat' in f][0], delim_whitespace=True).iloc[:,1]
            self.P[i] = pd.read_csv([f for f in nnfiles if 'density_hole.dat' in f][0], delim_whitespace=True).iloc[:,1]
            self.Efield[i][1::] = pd.read_csv([f for f in nnfiles if 'electric_field.dat' in f][0], delim_whitespace=True).iloc[:,1]
            self.Efield[i][0] = pd.read_csv([f for f in nnfiles if 'electric_field.dat' in f][0], delim_whitespace=True).iloc[0,1]
            self.mun[i][1::] = pd.read_csv([f for f in nnfiles if 'mobility_electron.dat' in f][0], delim_whitespace=True).iloc[:,1]
            self.mun[i][0] = pd.read_csv([f for f in nnfiles if 'mobility_electron.dat' in f][0], delim_whitespace=True).iloc[0,1]
            self.mup[i][1::] = pd.read_csv([f for f in nnfiles if 'mobility_hole.dat' in f][0], delim_whitespace=True).iloc[:,1]
            self.mup[i][0] = pd.read_csv([f for f in nnfiles if 'mobility_hole.dat' in f][0], delim_whitespace=True).iloc[0,1]
    
    def plot_results(self,V_idx=None, cmap = 'tab10'):
        """
        Plot simulation results in a 2x1 subplot layout.

        Args:
            V_idx (list, optional): Indices of voltages to plot. Defaults to first and last.
            colors (list, optional): List of colors for plotting. Defaults to default color cycle.

        Returns:
            tuple: Figure and axes objects
        """
        if V_idx is None:
            V_idx = [0, len(self.V)-1]

        cmap = plt.get_cmap(cmap)
        norm = plt.Normalize(vmin=min(V_idx), vmax=max(V_idx))
        colors = [cmap(norm(v)) for v in V_idx]
        
        if V_idx == None:
            V_idx = [0,len(self.V)-1]
            
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 8),sharex=True)
        # ax2r = ax2.twinx()
        for i, v in enumerate(V_idx):
            ax1.plot(self.grid, self.Ec[v], "-", color=colors[i], label=r"$E_c(eV) @ V={{{:.1f}}} V)$".format(self.V[v]))
            ax1.plot(self.grid, self.Ev[v], "-", color=colors[i])
            # Plot quasi-Fermi levels
            ax1.plot(self.grid, self.Efn[v], "-.", color=colors[i], linewidth=0.5)
            ax1.plot(self.grid, self.Efp[v], "-.", color=colors[i], linewidth=0.5)
        # Configure first subplot
                        
            # ax2 = ax1.twinx()
            ax2.plot(self.grid,  self.N[v],"-", color=colors[i], label=r"e conc. @ V={{{:.1f}}} V)$".format(self.V[v]))
            ax2.plot(self.grid, -self.P[v],"-.", color=colors[i], label=r"-h conc. @ V={{{:.1f}}} V)$".format(self.V[v]))
            
            ax3.plot(self.grid,self.Efield[v],label=r"EField@ V={{{:.1f}}} V)$".format(v), color = colors[i])
            ax3.set_ylim(-300,100)
            
        ax1.set_ylabel('Energy (eV)')
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='best')
        
        ax2.set_ylabel(r"Carrier conc. ($cm^{-3}$)")
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc='best')
        
        ax3.set_ylabel(r"Electric field (kV/cm)")
        ax3.grid(True, alpha=0.3)
        ax3.legend(loc='best')
        
    def transfer_results_to_device(self,
                        dx=0.05,
                        xmin=None,
                        xmax=None):
        
        """
        Interpolate 1D simulation data onto a new 2D mesh.

        .. note::
            This method for the moment it only interpolates the 1D data onto the horizontal dimension.

        Args:
            dx (float, optional): Step size for new mesh in microns. Defaults to 0.05.
            xmin (float): Minimum x value for mesh (required).
            xmax (float): Maximum x value for mesh (required).

        Returns:
            None. Stores interpolators in self.photonicdevice.charge.
        """
        
        if xmin is None or xmax is None:
            raise ValueError("Both xmin and xmax must be provided as numeric values. (e.g. waveguide boundaries)")
        
        reg = self.photonicdevice.reg
        # First part is to make data into 2d and fit the wg
        x = np.arange(xmin, xmax, dx)
        y = np.array(self.grid) * 1e-3  # Convert list to numpy array first

        xx, yy = np.meshgrid(x, y)

        # Initialize 2D arrays for each variable
        Ec_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        Ev_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        Efn_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        Efp_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        N_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        P_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        Efield_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        mun_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        mup_2d = np.zeros(shape=(len(self.V), len(y), len(x)))

        # Store coordinate grids
        self.x_2d = x
        self.y_2d = y
        self.xx_2d = xx
        self.yy_2d = yy

        # For each voltage, replicate 1D data across x-axis
        for i, v in enumerate(self.V):
            # Take 1D data (shape: n_y_points) and replicate across x-axis
            # Using broadcasting: 1D array becomes column, then broadcast to all x positions
            Ec_2d[i] = np.broadcast_to(self.Ec[i][:, np.newaxis], (len(y), len(x)))
            Ev_2d[i] = np.broadcast_to(self.Ev[i][:, np.newaxis], (len(y), len(x)))
            Efn_2d[i] = np.broadcast_to(self.Efn[i][:, np.newaxis], (len(y), len(x)))
            Efp_2d[i] = np.broadcast_to(self.Efp[i][:, np.newaxis], (len(y), len(x)))
            N_2d[i] = np.broadcast_to(self.N[i][:, np.newaxis], (len(y), len(x)))
            P_2d[i] = np.broadcast_to(self.P[i][:, np.newaxis], (len(y), len(x)))
            Efield_2d[i] = np.broadcast_to(self.Efield[i][:, np.newaxis], (len(y), len(x)))
            mun_2d[i] = np.broadcast_to(self.mun[i][:, np.newaxis], (len(y), len(x)))
            mup_2d[i] = np.broadcast_to(self.mup[i][:, np.newaxis], (len(y), len(x)))

        #Transform the Efield into a 3d vector field of shape (Ny, Nx, 3)
        Efield_2d = Efield_2d[..., np.newaxis]*self.sim_vector_norm
        
        #this part needs to poop out the interpolators 
        #if the interpolator is called the out of bound points should return the boundary values
            # Initialize interpolator dictionaries
        Ec_int = []
        Ev_int = []
        Efn_int = []
        Efp_int = []
        N_int = []
        P_int = []
        Efield_int = []
        mun_int = []
        mup_int = []
        
        for i ,v in enumerate(self.V):
            Ec_int.append(
                lambda x,y, arr=Ec_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.eV
            )

            Ev_int.append(
                lambda x,y, arr=Ev_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.eV
            )

            Efn_int.append(
                lambda x,y, arr=Efn_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.eV
            )

            Efp_int.append(
                lambda x,y, arr=Efp_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.eV
            )

            N_int.append(
                lambda x,y, arr=N_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.cm**-3
            )

            P_int.append(
                lambda x,y, arr=P_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.cm**-3
            )

            Efield_int.append(
                lambda x,y, arr=Efield_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.kV / reg.cm
            )

            mun_int.append(
                lambda x,y, arr=mun_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.cm**2 / reg.V / reg.s
            )

            mup_int.append(
                lambda x,y, arr=mup_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.cm**2 / reg.V / reg.s
            )
            
        self.photonicdevice.charge = {
            "Ec": Ec_int,
            "Ev": Ev_int,
            "Efn": Efn_int,
            "Efp": Efp_int,
            "N": N_int,
            "P": P_int,
            "Efield": Efield_int,
            "mun": mun_int,
            "mup": mup_int,
            "V": self.V,
        }
          
    def plot_with_simulation_line(
        self,
        color_polygon="black",
        color_line="green", 
        color_junctions="blue",
        color_simulation_line="red",
        fill_polygons=False,
        linewidth_simulation=2,
        fig=None,
        ax=None,
    ):
        """
        Plot the device polygons with the simulation line overlay.
        
        Args:
            color_polygon: Color for device polygons
            color_line: Color for current calculation lines  
            color_junctions: Color for junction regions
            color_simulation_line: Color for the simulation line
            fill_polygons: Whether to fill the polygons
            linewidth_simulation: Line width for simulation line
            fig: Existing figure object (optional)
            ax: Existing axis object (optional)
            
        Returns:
            fig, ax: matplotlib figure and axis objects
        """
        if fig is None and ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        for name, poly in self.polygon_entities.items():
            if isinstance(poly, Polygon):
                ax.plot(
                    *poly.exterior.xy,
                    color=color_polygon if "junction" not in name else color_junctions,
                )
                if fill_polygons:
                    ax.fill(
                        *poly.exterior.xy,
                        color=np.random.rand(3,),
                        alpha=0.5,
                    )
                    
            elif isinstance(poly, Line):
                ax.plot(*poly.xy, color=color_line)

        
        # Plot the simulation line if it exists
        if hasattr(self, 'simulation_line') and self.simulation_line is not None:
            ax.plot(
                *self.simulation_line.xy, 
                color=color_simulation_line, 
                linewidth=linewidth_simulation,
                label="Simulation Line"
            )
            ax.legend()
            
            if len(self.line_segments) > 0:
                ax.legend()
        
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_title("PhotonicDevice with Simulation Line")
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        return fig, ax
    



class CustomMaterial_OBP(BaseMaterial):

    def __init__(self, name, opb_mat, T = 300, **kwargs):
        BaseMaterial.__init__(
            self,
            T=300, 
            material_string = name,
            **kwargs
        )

        ########## Universal constants ##########
        kb = solcore.constants.kb  # Boltzmann constant in J/K
        hbar = solcore.constants.hbar  # Planck's constant in J·s
        h = hbar * 2 * np.pi  # Convert hbar to Planck's constant
        m0 = solcore.constants.electron_mass  # Electron rest mass in kg
        q = solcore.constants.q  # Elementary charge in C

        self.name = name

        mat_name = name
        PARAMETER_PATH = os.path.join(config.user_folder, "custom_parameters.txt")
        if "custom" not in config.parameters():
            if not os.path.isfile(PARAMETER_PATH):
                open(PARAMETER_PATH, "a").close()
            config["Parameters", "custom"] = PARAMETER_PATH

        CUSTOM_PATH = os.path.join(config.user_folder, "custom_materials")

        # create a folder in the custom materials folders
        folder = os.path.join(CUSTOM_PATH, mat_name + "-Material")
        if not os.path.exists(folder) and folder != "":
            os.makedirs(folder)

        # Save some dummy data in the folder
        n = np.asarray([np.linspace(300, 2000, 10), np.ones(10)]).T
        k = np.asarray([np.linspace(300, 2000, 10), np.zeros(10)]).T

        np.savetxt(os.path.join(folder, "n.txt"), n)
        np.savetxt(os.path.join(folder, "k.txt"), k)

        config["Materials", mat_name] = folder

        self.obp_mat = opb_mat
        self.Na = kwargs.get('Na',1)
        self.Nd = kwargs.get('Nd', 1)
        self.T = T

        # solcore_property: openbandparams_property
        # Map the properties from openbandparams to solcore
        property_map = {
            'lattice_constant': 'a',
            'Eg_Gamma': 'Eg_Gamma',
            'Eg_X': 'Eg_X',
            'Eg_L': 'Eg_L',
            'spin_orbit_splitting': 'Delta_SO',
            'interband_matrix_element': 'Ep',
            'gamma1': 'luttinger1',
            'gamma2': 'luttinger2',
            'gamma3': 'luttinger3',
            'relative_permittivity': 'dielectric',
            'eff_mass_split_off': 'meff_SO',
            'valence_band_offset': 'VBO',
            'a_c': 'a_c',
            'a_v': 'a_v',
            'c11': 'c11',
            'c12': 'c12',
            'c44': 'c44',
            'b': 'b',
            'd': 'd',
            'F': 'F',
            'eff_mass_electron_Gamma': 'meff_e_Gamma',
            'eff_mass_electron_long_L': 'meff_e_L_long',
            'eff_mass_electron_trans_L': 'meff_e_L_trans',
            'eff_mass_electron_long_X': 'meff_e_X_long',
            'eff_mass_electron_trans_X': 'meff_e_X_trans',
            'eff_mass_DOS_L': 'meff_e_L_DOS',
            'eff_mass_DOS_X': 'meff_e_X_DOS',
            'electron_affinity': 'electron_affinity'
        }

        # Loop over attributes and assign them
        for attr, func_name in property_map.items():
            try:
                if func_name == 'meff_SO':  # special case
                    value = opb_mat.meff_SO(literature_value=True, T=T)
                else:
                    func = getattr(opb_mat, func_name)
                    value = func(T=T)

                if attr in [
                    'Eg_Gamma', 
                    'Eg_X', 
                    'Eg_L', 
                    'spin_orbit_splitting', 
                    'interband_matrix_element',
                    'valence_band_offset',
                    'a_c',
                    'a_v',
                    'b',
                    'd',
                    'electron_affinity',
                ]:
                    setattr(self, attr, value*q) #Energies are stored in eV in obp and J in solcore

                elif attr in [
                    'c11',
                    'c12',
                    'c44',
                ]:
                    setattr(self, attr, value*1e9)  # Convert GPa to Pa for solcore

                elif attr in [
                    'lattice_constant'
                ]:
                    setattr(self, attr, value*1e-10)  # Convert Angstrom to m for solcore
                else:
                    setattr(self, attr, value)
            except AttributeError:
                # print(f"The material {name} does not have a defined '{func_name}'")
                pass


        self.band_gap = min(self.Eg_Gamma, self.Eg_X, self.Eg_L)
        self.lowest_band = ['Gamma', 'X', 'L'][[self.Eg_Gamma, self.Eg_X, self.Eg_L].index(self.band_gap)]
        self.m0 = self.eff_mass_split_off*(self.gamma1 - self.interband_matrix_element*self.spin_orbit_splitting/(3*self.band_gap*(self.band_gap+self.spin_orbit_splitting)))
        self.eff_mass_hh_z = opb_mat.meff_hh_100(T = T)
        self.eff_mass_hh_110 = opb_mat.meff_hh_110(T = T)
        self.eff_mass_hh_111 = opb_mat.meff_hh_111(T = T)
        self.eff_mass_lh_z = opb_mat.meff_lh_100(T = T)
        self.eff_mass_lh_110 = opb_mat.meff_lh_110(T = T)
        self.eff_mass_lh_111 = opb_mat.meff_lh_111(T = T)
        self.eff_mass_electron = opb_mat.meff_e_Gamma(T = T)
        self.permittivity = opb_mat.dielectric(T = T)*8.854187817e-12

        self.electron_mobility, self.hole_mobility = self._calculate_mobility()

        default_properties = DefaultProperties

        #Missing parameters
        #electron_minority_lifetime
        #hole_minority_lifetime
        #electron_auger_recombination
        #hole_auger_recombination
        #radiative_recombination

        # Calculate the ni, Nc and Nv

        me = self.eff_mass_electron * m0
        self.Nc = 2*(2*np.pi*me*kb*T/h**2)**(3/2)  # Effective density of states in the conduction band in m^-3

        mh = (self.eff_mass_hh_z**(3/2)+self.eff_mass_lh_z**(3/2))**(2/3) * m0
        self.Nv = 2*(2*np.pi*mh*kb*T/h**2)**(3/2) # Effective density of states in the valence band in m^-3

        self.ni = np.sqrt(self.Nc) * np.sqrt(self.Nv) * np.exp(-self.band_gap/((2*kb*T)))  # Intrinsic carrier concentration in m^-3

        #InP values
        self.electron_minority_lifetime = kwargs.get('electron_minority_lifetime', 1e-9)
        self.hole_minority_lifetime = kwargs.get('hole_minority_lifetime', 1e-9)
        self.electron_auger_recombination = kwargs.get('electron_auger_recombination', 9e-31*1e-12)  # Convert from cm^6/s to m^6/s. Taken from ioffe
        self.hole_auger_recombination = kwargs.get('hole_auger_recombination', 9e-31*1e-12)  # Convert from cm^6/s to m^6/s. Taken from ioffe
        self.radiative_recombination = kwargs.get('radiative_recombination', 1.2e-10*1e-6)  # Convert from cm^3/s to m^3/s. Taken from ioffe


        new_params = [
            'band_gap',
            'lowest_band',
            'm0',
            'eff_mass_hh_z',
            'eff_mass_hh_110',
            'eff_mass_hh_111',
            'eff_mass_lh_z',
            'eff_mass_lh_110',
            'eff_mass_lh_111',
            'eff_mass_electron',
            'permittivity',
            'electron_mobility',
            'hole_mobility'
        ]

        units_dict = {
            'lattice_constant': 'Angstrom',
            'Eg_Gamma': 'eV',
            'Eg_X': 'eV',
            'Eg_L': 'eV',
            'spin_orbit_splitting': 'eV',
            'interband_matrix_element': 'eV',
            'gamma1': '',
            'gamma2': '',
            'gamma3': '',
            'relative_permittivity': '',
            'eff_mass_split_off': '',
            'valence_band_offset': 'eV',
            'a_c': 'eV',
            'a_v': 'eV',
            'c11': 'GPa',
            'c12': 'GPa',
            'c44': 'GPa',
            'b': 'eV',
            'd': 'eV',
            'F': '',
            'eff_mass_electron_Gamma': '',
            'eff_mass_electron_long_L': '',
            'eff_mass_electron_trans_L': '',
            'eff_mass_electron_long_X': '',
            'eff_mass_electron_trans_X': '',
            'eff_mass_DOS_L': '',
            'eff_mass_DOS_X': '',
            'electron_affinity': 'eV',
            'band_gap': 'eV',
            'lowest_band': '',
            'm0': '',
            'eff_mass_hh_z': '',
            'eff_mass_hh_110': '',
            'eff_mass_hh_111': '',
            'eff_mass_lh_z': '',
            'eff_mass_lh_110': '',
            'eff_mass_lh_111': '',
            'eff_mass_electron': '',
            'permittivity': 'F/m',
            'electron_mobility': 'm^2/(V·s)',
            'hole_mobility': 'm^2/(V·s)'
        }

        all_params = list(property_map.keys()) + new_params

        parameter_source = 'params_file_tmp.txt'
        with open(parameter_source, 'w') as f:
            f.write(f"[{mat_name}]\n")
            for param in all_params:
                try:
                    value = object.__getattribute__(self, param)
                except AttributeError:
                    value = None

                if param in [
                    'Eg_Gamma', 
                    'Eg_X', 
                    'Eg_L', 
                    'spin_orbit_splitting', 
                    'interband_matrix_element',
                    'valence_band_offset',
                    'a_c',
                    'a_v',
                    'b',
                    'd',
                    'electron_affinity',
                ]:
                    value = value/q #while energies are stored in J in solcore, we need to write this file with energies in eV

                elif param in [
                    'c11',
                    'c12',
                    'c44',
                ]:
                    value = value*1e-9 # While these values are stored in Pa in solcore, we need to write this file with GPa

                elif param in [
                    'lattice_constant'
                ]:
                    value = value/1e-10  # while these values are stored in m in solcore, we need to write this file with Angstrom

                if value is not None:
                    if isinstance(value, (int, float)):
                        f.write(f"{param}={value:.2f} {units_dict[param]}\n")
                        
                    elif isinstance(value, str):
                        f.write(f"{param}='{value} {units_dict[param]}'\n")

        # append the parameters for the new material
        params = ConfigParser()
        params.optionxform = str
        if parameter_source is not None:
            params.read([PARAMETER_PATH, parameter_source])
            with open(PARAMETER_PATH, "w") as fp:
                params.write(fp)
        # print(PARAMETER_PATH)
        # print(CUSTOM_PATH)

        ParameterSystem().read()

        if os.path.exists(parameter_source):
            os.remove(parameter_source)
                    

    def _calculate_mobility(self):

        N = float(self.Nd)  # Carrier concentration in cm^-3
        P = float(self.Na) # Hole concentration in cm^-3
        T = self.T  # Temperature in K

        elements = self.obp_mat.elements
        fractions = [self.obp_mat.element_fraction(element) for element in elements]

        binary = False
        ternary = False
        quaternary = False

        if len(elements) == 2:
            binary = True
        elif len(elements) == 3:
            ternary = True
        elif len(elements) == 4:
            quaternary = True

        #From the elements of the openbandparams material, we will see which alloy it corresponds to on the solcore side

        if quaternary:
            if set(elements) == set(['In', 'Ga', 'As', 'P']):
                d_e = calculate_InGaAsP(
                    x = fractions[elements.index('In')],
                    y = fractions[elements.index('P')],
                    i = 1,
                    T = 300
                )

                d_h = calculate_InGaAsP(
                    x = fractions[elements.index('In')],
                    y = fractions[elements.index('P')],
                    i = 2,
                    T = 300
                )

            else:
                raise ValueError(f"Quaternary alloy {elements} not recognised. It is not allowed at the moment")
            
        elif ternary:
            if set(elements) == set(['In', 'Ga', 'As']):
                d_e = calculate_InGaAs(
                    x = fractions[elements.index('In')],
                    i = 1,
                )

                d_h = calculate_InGaAs(
                    x = fractions[elements.index('In')],
                    i = 2,
                )

            elif set(elements) == set(['Ga', 'In', 'P']):
                d_e = calculate_InGaP(
                    x = fractions[elements.index('In')],
                    i = 1,
                    T = 300
                )

                d_h = calculate_InGaP(
                    x = fractions[elements.index('In')],
                    i = 2,
                    T = 300
                )
            elif set(elements) == set(['Al', 'Ga', 'As']):
                d_e = calculate_AlGaAs(
                    x = fractions[elements.index('Al')],
                    i = 1,
                    T = 300
                )

                d_h = calculate_AlGaAs(
                    x = fractions[elements.index('Al')],
                    i = 2,
                    T = 300
                )
            elif set(elements) == set(['Al', 'In', 'As']):
                d_e = calculate_InAlAs(
                    x = fractions[elements.index('Al')],
                    i = 1,
                    T = 300
                )

                d_h = calculate_InAlAs(
                    x = fractions[elements.index('Al')],
                    i = 2,
                    T = 300
                )
            else:
                raise ValueError(f"Ternary alloy {elements} not recognised. It is not allowed at the moment")
            
        elif binary:
            this_dir = os.path.dirname(inspect.getfile(solcore.material_data))
            parameters = os.path.join(this_dir, "mobility_parameters.json")
            f = open(parameters, mode="r")
            data = json.load(f)

            if set(elements) == set(['In', 'As']):
                d_e = data["InAs"][1]
                d_h = data["InAs"][2]
            elif set(elements) == set(['Al', 'As']):
                d_e = data["AlAs"][1]
                d_h = data["AlAs"][2]
            elif set(elements) == set(['Ga', 'As']):
                d_e = data["GaAs"][1]
                d_h = data["GaAs"][2]
            elif set(elements) == set(['In', 'P']):
                d_e = data["InP"][1]
                d_h = data["InP"][2]
            elif set(elements) == set(['Ga', 'P']):
                d_e = data["GaP"][1]
                d_h = data["GaP"][2]
            

        #electrons
        muMin_e = d_e["muMin"]
        muMax_e = d_e["muMax"]
        Nref_e = d_e["Nref"]
        l_e = d_e["l"]
        t1_e = d_e["t1"]
        t2_e = d_e["t2"]

        m_e = mobility_low_field(N / 1e6, muMin_e, muMax_e, Nref_e, l_e, t1_e, t2_e, T) / 10000  # To convert it from cm2 to m2

        muMin_h = d_h["muMin"]
        muMax_h = d_h["muMax"]
        Nref_h = d_h["Nref"]
        l_h = d_h["l"]
        t1_h = d_h["t1"]
        t2_h = d_h["t2"]

        m_h = mobility_low_field(P / 1e6, muMin_h, muMax_h, Nref_h, l_h, t1_h, t2_h, T) / 10000  # To convert it from cm2 to m2

        return m_e, m_h
    
    def print_properties(self):
        params = [
        'lattice_constant'
        ,'Eg_Gamma'
        ,'Eg_X'
        ,'Eg_L'
        ,'spin_orbit_splitting'
        ,'interband_matrix_element'
        ,'gamma1'
        ,'gamma2'
        ,'gamma3'
        ,'relative_permittivity'
        ,'eff_mass_split_off'
        ,'valence_band_offset'
        ,'a_c'
        ,'a_v'
        ,'c11'
        ,'c12'
        ,'c44'
        ,'b'
        ,'d'
        ,'F'
        ,'eff_mass_electron_Gamma'
        ,'eff_mass_electron_long_L'
        ,'eff_mass_electron_trans_L'
        ,'eff_mass_electron_long_X'
        ,'eff_mass_electron_trans_X'
        ,'eff_mass_DOS_L'
        ,'eff_mass_DOS_X'
        ,'electron_affinity'
        ,'band_gap'
        ,'lowest_band'
        ,'m0'
        ,'eff_mass_hh_z'
        ,'eff_mass_hh_110'
        ,'eff_mass_hh_111'
        ,'eff_mass_lh_z'
        ,'eff_mass_lh_110'
        ,'eff_mass_lh_111'
        ,'eff_mass_electron'
        ,'permittivity'
        ,'electron_mobility'
        ,'hole_mobility'
        ,'electron_minority_lifetime'
        ]

        print(f'********************{self.name}*******************')
        for prop in params:
            try:
                value = getattr(self, prop)
                print(f'{prop:25}: {value:15.2E}')
            except:
                print(f'{prop:15} not found in {self.name}')



class ChargeSimulatorSolcore:
    '''

        Each SemiconductorPolygon can have the following charge_transport_kwargs:
            - sol_obp_material: a openbandparams material object,
            - sol_Na: concentration of acceptors in cm^-3,
            - sol_Nd: concentration of donors in cm^-3,
            - sol_electron_minority_lifetime: electron minority lifetime in seconds,
            - sol_hole_minority_lifetime: hole minority lifetime in seconds,

        .. warning::
            This solver has a current issue with the ordering of the voltages when it outputs the results. Please visit https://github.com/qpv-research-group/solcore5/issues/294 for more details.
    '''

    def __init__(
        self,
        device: PhotonicDevice, 
        simulation_line: LineString,
        temperature: float = 300.0,  # Add temperature parameter
        bias_start_stop_step: list = [0,1,1], #contact1 is the bias electrode decide - or + accordingly
    ):
        """
        Initialize the ChargeSimulatorNN.

        Args:
            device: PhotonicDevice instance containing the device geometry and materials.
            simulation_line: LineString defining the simulation line along which to perform 1D simulation.
            temperature: Simulation temperature in Kelvin. Defaults to 300.0.
            bias_start_stop_step: Voltage sweep [start, stop, steps]. Defaults to [0,1,1].
            
        """
        self.temperature = temperature
        self.photonicdevice = device
        self.bias_start_stop_step=bias_start_stop_step

        self.optical_photopolygons = copy.deepcopy(self.photonicdevice.photo_polygons)

        self.polygon_entities = OrderedDict()

        for polygon in self.optical_photopolygons:
            self.polygon_entities[polygon.name] = polygon.polygon
        
        self.sim_vector_norm = get_normalized_vector(simulation_line)
        self.sim_vector_norm = np.array([self.sim_vector_norm[0], self.sim_vector_norm[1], 0])

        self._select_line(simulation_line=simulation_line)
        self._create_materials()
        self._create_mesh()
        self._create_junction()
        
    def _create_junction(self):
        """
            Create the solcore.Junction object that will be passed to the solcore PoissonDriftDiffusion simulator.

            It will update the attribute self.junction
        """
        junction_layers = []
        for name, line in self.line_segments.items():

            x0, x1 = line.xy[0]
            y0, y1 = line.xy[1]
            D = y1 - y0  # Distance between the two points in the y direction
            layer = Layer(D*1e-6, material=self.solcore_materials[name])
            junction_layers.append(layer)

        # Create the junction
        self.junction = Junction(
            junction_layers,
            kind = 'sesame_PDD'
        )

        #Add the mesh to the junction
        self.junction.mesh = self.mesh*1e-6 # Convert from um to m

        ## For some reason the solver always crashes if the mesh does not start at 0. Therefore, we will set the initial point of the mesh to 0.
        # self.junction.mesh -= np.min(self.junction.mesh)
        self.junction.mesh -= min(self.junction.mesh, key=lambda x: abs(x - self.simulation_line.coords[0][1]*1e-6))

    def _create_materials(self):
        """
        Create solcore materials for each line segment in the device

        This method extracts the material properties from the PhotonicPolygon's charge_transport_simulator_kwargs
        and creates a CustomMaterial_OBP for each segment. The materials are stored in the
        self.solcore_materials dictionary, keyed by the segment name.
        """

        solcore_materials = OrderedDict()
        for name, segment in self.line_segments.items():
            photo_poly = next((photo_poly for photo_poly in self.optical_photopolygons if photo_poly.name == name), None)

            obp_mat = photo_poly.charge_transport_simulator_kwargs.get('sol_obp_material', None)
            if obp_mat is None:
                raise ValueError(f"PhotonicPolygon {photo_poly.name} does not have a valid 'sol_obp_material' in charge_transport_kwargs.")
            sol_Na = photo_poly.charge_transport_simulator_kwargs.get('sol_Na', 1)*1e6 # # Convert from cm^-3 to m^-3
            sol_Nd = photo_poly.charge_transport_simulator_kwargs.get('sol_Nd', 1)*1e6 # # Convert from cm^-3 to m^-3
            sol_electron_minority_lifetime = photo_poly.charge_transport_simulator_kwargs.get('sol_electron_minority_lifetime', 1e-9)
            sol_hole_minority_lifetime = photo_poly.charge_transport_simulator_kwargs.get('sol_hole_minority_lifetime', 1e-9)

            #No need to default the minority lifetimes as they are already defaulted in the CustomMaterial_OBP class

            # Create the material
            if sol_electron_minority_lifetime is not None and sol_hole_minority_lifetime is not None:
                material = CustomMaterial_OBP(
                    name=name,
                    opb_mat=obp_mat,
                    Na=sol_Na,
                    Nd=sol_Nd,
                    electron_minority_lifetime=sol_electron_minority_lifetime,
                    hole_minority_lifetime=sol_hole_minority_lifetime,
                )
            elif sol_electron_minority_lifetime is not None and sol_hole_minority_lifetime is None:
                material = CustomMaterial_OBP(
                    name=name,
                    opb_mat=obp_mat,
                    Na=sol_Na,
                    Nd=sol_Nd,
                    electron_minority_lifetime=sol_electron_minority_lifetime
                )
            elif sol_electron_minority_lifetime is None and sol_hole_minority_lifetime is not None:
                material = CustomMaterial_OBP(
                    name=name,
                    opb_mat=obp_mat,
                    Na=sol_Na,
                    Nd=sol_Nd,
                    hole_minority_lifetime=sol_hole_minority_lifetime
                )
            else:
                material = CustomMaterial_OBP(
                    name=name,
                    opb_mat=obp_mat,
                    Na=sol_Na,
                    Nd=sol_Nd,
                )

            # Add the material to the dictionary
            solcore_materials[name] = material

        # Store the materials in the simulator
        self.solcore_materials = solcore_materials

    def _create_mesh(self):
        """
        Here we create a mesh for the device based on the line segments. 

        This meshing is handled by gmsh.

        The mesh is stored in self.mesh as a numpy array of x-coordinates.

        """
        
        gmsh.initialize()

        gmsh.model.add('mesh1d')

        point_tags = []
        #Add all the points that comprise our line to the GMSH model
        for name, segment in self.line_segments.items():
            x0, x1 = segment.xy[0]
            y0, y1 = segment.xy[1]

            tag = gmsh.model.geo.addPoint(y0, 0, 0, self.photonicdevice.resolutions_charge[name]['resolution'])
            
            point_tags.append(tag)

        #Add the final point
        x0, x1 = segment.xy[0]
        y0, y1 = segment.xy[1]

        tag = gmsh.model.geo.addPoint(y1, 0, 0, self.photonicdevice.resolutions_charge[name]['resolution'])

        point_tags.append(tag)

        #Add the lines
        for i in range(len(point_tags) - 1):
            gmsh.model.geo.addLine(point_tags[i], point_tags[i+1])

        gmsh.model.geo.synchronize()

        gmsh.model.mesh.generate(1)
        gmsh.model.mesh.removeDuplicateNodes()

        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        node_coords = node_coords.reshape(-1, 3)[:, 0]  # only x-coordinates

        node_coords = sorted(node_coords)
        gmsh.finalize()

        self.mesh = np.asarray(node_coords)

            

    def _select_line(self,
                simulation_line: LineString, #to select the region
                ): 
        """
        Find intersections between the simulation line and device polygons.

        Args:
            simulation_line (LineString): The line along which to perform 1D simulation.

        Returns:
            None. Stores intersection segments in self.line_segments.
        """
        line_segments = OrderedDict()

        names_to_print = []
        for polygon_name, polygon_geom in self.polygon_entities.items():
            # Find intersection between simulation line and polygon
            
            intersection = simulation_line.intersection(polygon_geom)
            if (polygon_name == "substrate") or (polygon_name == "background"):
                continue
            # Handle different intersection types
            if intersection.is_empty:
                continue

            # If intersection is a single LineString
            if intersection.geom_type == 'LineString':
                line_segments[f'{polygon_name}'] = intersection

                ## Loop over the photopolygons of the PhotonicDevice to find the proper PhotonicPolygon that will allow for charge transport data to be loaded on:
                for poly in self.photonicdevice.photo_polygons:
                    if poly.name == polygon_name:
                        poly.has_charge_transport_data = True

                names_to_print.append(polygon_name)

        print('Charge transport will take place with:')
        print(*names_to_print, sep='\n')

        #The order of the line segments in line_Segments is dependent on the order of the self.polygon_entities. However, in this case, we must ensure that the line segments are stored in such a way that the polygons to which they belong appear as if we are to walk along the beggining of the simulation line to the end.

        #To do this we will simple walk the line: we start at point 0 of the line. Then we find the point halfway through point 0 and 1, and find to which polygon it belongs to. Then we move on to point 1 and find the one halfway between 1 and 2 and look for the corresponding polygon. We continue until we're done.
        
        start_point = np.asarray(simulation_line.xy).T[0]
        end_point = np.asarray(simulation_line.xy).T[1]

        new_line_segments = OrderedDict()

        current_point = start_point #This will be updated as we walk the line

        #The starting point must appear only once, so finding the first line segment is easy.
        #For the second point we will loop over all the line segments except those already found, and find the one that contains the current point
        k=0
        while not np.isclose(current_point, end_point).all():
            for name, line in line_segments.items():

                if name in new_line_segments:
                    continue
                # print(name, np.asarray(line.xy))
                points = np.asarray(line.xy).T

                p0 = points[0]
                p1 = points[1]

                if np.allclose(p0, current_point) or np.allclose(p1, current_point):
                    new_line_segments[name] = line
                    current_point = p1 if np.allclose(p0, current_point) else p0
                    break

            k += 1
            if k > 10_000: #This is just a timeout safety measure
                raise ValueError('Something went wrong, could not find all line segments')

        # Store the line segments for later use
        self.line_segments = new_line_segments
        self.simulation_line = simulation_line

    def solve_PDD(
            self,
            verbose = False,
            tol = 1e-6,
            max_iter = 100,
            htp = 1,
            smooth_output = False,
            savgol_window = 15,
    ):
        """
        Solve the Poisson Drift Diffusion equations for the device.

        .. warning::
            The current implementation of the sesame pdd in solcore has a small bug that can cause the outputs to be switched around. Because the solver must first solve the 0V and then it goes to some voltage, the solver does some reordering of the voltages to make sure that the solver converges. However, when storing the voltages corresponding to the outputs in the junction object, it is not saving the correct voltages yet, and for that reason, the outputs can be switched around. For more details go to: https://github.com/qpv-research-group/solcore5/issues/294

        This function will update self object with the following attributes:
            - self.V: The potential distribution in the device.
            - self.Ec: The conduction band edge energy distribution in the device.
            - self.Ev: The valence band edge energy distribution in the device.
            - self.Efn: The Fermi level energy distribution in the device.
            - self.Efp: The Fermi level energy distribution in the device.
            - self.N: The electron concentration distribution in the device.
            - self.P: The hole concentration distribution in the device.
            - self.Efield: The electric field distribution in the device.
            - self.mun: The electron mobility distribution in the device.
            - self.mup: The hole mobility distribution in the device.

        Args:
            verbose (bool): Whether to print detailed output.
            tol (float): Tolerance for convergence.
            max_iter (int): Maximum number of iterations.
            htp (int): Number of homotopic Newton loops to perform

            .. Note::
               For more details on these arguments please consult solcore.sesame_drift_diffusion.solve_pdd 

            smooth_output (bool): If True, it will apply a savitzky-golay filter to each of the attributes specified above. The reasoning behind this is to numerical artifacts with big spikes in the charge distributions and electrical fields. For more details please check https://github.com/qpv-research-group/solcore5/discussions/293
            savgol_window (int): The window size for the savitzky-golay filter. Default is 15. It should be an odd number.

        Returns:
            None

        """

        options = State()
        options.verbose = verbose
        options.T = self.temperature

        Vstart, Vstop, Vsteps = self.bias_start_stop_step
        options.voltages = np.linspace(Vstart, Vstop, Vsteps)
        options.internal_voltages = options.voltages.copy()

        #PDD solver options
        options.sesame_verbose = verbose
        options.sesame_tol = tol
        options.sesame_max_iterations = max_iter
        options.sesame_htp = htp
        options.sesame_periodic = False

        solar_cell_sesame = SolarCell([self.junction])
        solar_cell_solver(solar_cell_sesame, 'iv', options)

        if smooth_output:
            from scipy.signal import savgol_filter
            
            for key in self.junction.pdd_output.keys():
                if key in ['Ec', 'Ev', 'Efe', 'Efh', 'potential', 'n', 'p']:
                    for idx_voltage in range(len(self.junction.pdd_output[key])):
                        self.junction.pdd_output[key][idx_voltage] = savgol_filter(
                            self.junction.pdd_output[key][idx_voltage], savgol_window, 1
                        )

        # Extract the results from the solar cell solver
        #Process the output to retrieve the mobilities
        junction = self.junction
        un = np.zeros_like(junction.mesh)
        up = np.zeros_like(junction.mesh)

        x0 = junction.mesh[0]
        for i, layer in enumerate(junction):
            width = layer.width
            up_layer = layer.material.hole_mobility
            un_layer = layer.material.electron_mobility

            if i == len(junction) - 1:
                idx = np.where((junction.mesh >= x0) & (junction.mesh <= x0 + width))
            else:
                idx = np.where((junction.mesh >= x0) & (junction.mesh < x0 + width))
            up[idx] = up_layer
            un[idx] = un_layer
            x0 += width

        un = np.broadcast_to(un[np.newaxis,:], (len(junction.pdd_output['potential']), len(junction.mesh)))
        up = np.broadcast_to(up[np.newaxis, :], (len(junction.pdd_output['potential']), len(junction.mesh)))

        # Retrieve the electric field
        E_field = np.zeros_like(junction.pdd_output['potential'])
        for i in range(len(junction.pdd_output['potential'])):
            E_field[i] = -np.gradient(junction.pdd_output['potential'][i], junction.mesh)

        self.Ec = junction.pdd_output['Ec']
        self.Ev = junction.pdd_output['Ev']
        self.Efn = junction.pdd_output['Efe']
        self.Efp = junction.pdd_output['Efh']
        self.V = junction.voltage
        self.N = junction.pdd_output['n']/1e6 # Convert from m^-3 to cm^-3
        self.P = junction.pdd_output['p']/1e6 # Convert from m^-3 to cm^-3
        self.Efield = E_field / 1e5 # Convert from V/m to kV/cm
        self.mun = un * 1e4 # Convert from m^2/(V·s) to cm^2/(V·s)
        self.mup = up * 1e4 # Convert from m^2/(V·s) to cm^2/(V·s)

    def transfer_results_to_device(self, dx = 0.05, xmin = -2, xmax = 2):
        
        """
        Interpolate 1D simulation data onto a new 2D mesh.

        .. note::
            This method for the moment it only interpolates the 1D data onto the horizontal dimension.

        Args:
            dx (float, optional): Step size for new mesh in microns. Defaults to 0.05.
            xmin (float): Minimum x value for mesh (required).
            xmax (float): Maximum x value for mesh (required).

        Returns:
            None. Stores interpolators in self.photonicdevice.charge.
        """

        reg = self.photonicdevice.reg
        # First part is to make data into 2d and fit the wg
        x = np.arange(xmin, xmax, dx)
        y = np.array(self.mesh)  # Convert list to numpy array first

        xx, yy = np.meshgrid(x, y)

        # Initialize 2D arrays for each variable
        Ec_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        Ev_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        Efn_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        Efp_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        N_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        P_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        Efield_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        mun_2d = np.zeros(shape=(len(self.V), len(y), len(x)))
        mup_2d = np.zeros(shape=(len(self.V), len(y), len(x)))

        # Store coordinate grids
        self.x_2d = x
        self.y_2d = y
        self.xx_2d = xx
        self.yy_2d = yy

        # For each voltage, replicate 1D data across x-axis
        for i, v in enumerate(self.V):
            # Take 1D data (shape: n_y_points) and replicate across x-axis
            # Using broadcasting: 1D array becomes column, then broadcast to all x positions
            Ec_2d[i] = np.broadcast_to(self.Ec[i][:, np.newaxis], (len(y), len(x)))
            Ev_2d[i] = np.broadcast_to(self.Ev[i][:, np.newaxis], (len(y), len(x)))
            Efn_2d[i] = np.broadcast_to(self.Efn[i][:, np.newaxis], (len(y), len(x)))
            Efp_2d[i] = np.broadcast_to(self.Efp[i][:, np.newaxis], (len(y), len(x)))
            N_2d[i] = np.broadcast_to(self.N[i][:, np.newaxis], (len(y), len(x)))
            P_2d[i] = np.broadcast_to(self.P[i][:, np.newaxis], (len(y), len(x)))
            Efield_2d[i] = np.broadcast_to(self.Efield[i][:, np.newaxis], (len(y), len(x)))
            mun_2d[i] = np.broadcast_to(self.mun[i][:, np.newaxis], (len(y), len(x)))
            mup_2d[i] = np.broadcast_to(self.mup[i][:, np.newaxis], (len(y), len(x)))

        # Now we need to mask the values to include only points inside the polygons involved in the charge transport simulations
        points = shapely.points(xx, yy)
        total_mask = np.zeros(points.shape, dtype=bool)

        for poly_name in self.line_segments.keys():
            photo_poly = next((x for x in self.photonicdevice.photo_polygons if x.name == poly_name), None)
            poly = photo_poly.polygon

            # Compute mask once
            mask_inside = shapely.covers(poly, points)

            total_mask += mask_inside

        for data in [
            Ec_2d,
            Ev_2d,
            Efn_2d,
            Efp_2d,
            N_2d,
            P_2d,
            Efield_2d,
            mun_2d,
            mup_2d,
        ]:
            for i in range(len(self.V)):
                data[i] *= total_mask  # Mask out points outside polygons

        #Transform the Efield into a 3d vector field of shape (Ny, Nx, 3)
        Efield_2d = Efield_2d[..., np.newaxis]*self.sim_vector_norm
        
        #this part needs to poop out the interpolators 
        #if the interpolator is called the out of bound points should return the boundary values
            # Initialize interpolator dictionaries
        Ec_int = []
        Ev_int = []
        Efn_int = []
        Efp_int = []
        N_int = []
        P_int = []
        Efield_int = []
        mun_int = []
        mup_int = []
        
        for i ,v in enumerate(self.V):
            Ec_int.append(
                lambda x,y, arr=Ec_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.eV
            )

            Ev_int.append(
                lambda x,y, arr=Ev_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.eV
            )

            Efn_int.append(
                lambda x,y, arr=Efn_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.eV
            )

            Efp_int.append(
                lambda x,y, arr=Efp_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.eV
            )

            N_int.append(
                lambda x,y, arr=N_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.cm**-3
            )

            P_int.append(
                lambda x,y, arr=P_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.cm**-3
            )

            Efield_int.append(
                lambda x,y, arr=Efield_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.kV / reg.cm
            )

            mun_int.append(
                lambda x,y, arr=mun_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.cm**2 / reg.V / reg.s
            )

            mup_int.append(
                lambda x,y, arr=mup_2d[i]: RegularGridInterpolator(
                    (self.y_2d, self.x_2d),
                    arr,
                    method='linear',
                    bounds_error=False,
                    fill_value=None,
                )((y, x)) * reg.cm**2 / reg.V / reg.s
            )
            
        self.photonicdevice.charge = {
            "Ec": Ec_int,
            "Ev": Ev_int,
            "Efn": Efn_int,
            "Efp": Efp_int,
            "N": N_int,
            "P": P_int,
            "Efield": Efield_int,
            "mun": mun_int,
            "mup": mup_int,
            "V": self.V,
        }


    def plot_with_simulation_line(
        self,
        color_polygon="black",
        color_line="green", 
        color_junctions="blue",
        color_simulation_line="red",
        fill_polygons=False,
        linewidth_simulation=2,
        fig=None,
        ax=None,
    ):
        """
        Plot the device polygons with the simulation line overlay.
        
        Args:
            color_polygon: Color for device polygons
            color_line: Color for current calculation lines  
            color_junctions: Color for junction regions
            color_simulation_line: Color for the simulation line
            fill_polygons: Whether to fill the polygons
            linewidth_simulation: Line width for simulation line
            fig: Existing figure object (optional)
            ax: Existing axis object (optional)
            
        Returns:
            fig, ax: matplotlib figure and axis objects
        """
        if fig is None and ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        for name, poly in self.polygon_entities.items():
            if isinstance(poly, Polygon):
                ax.plot(
                    *poly.exterior.xy,
                    color=color_polygon if "junction" not in name else color_junctions,
                )
                if fill_polygons:
                    ax.fill(
                        *poly.exterior.xy,
                        color=np.random.rand(3,),
                        alpha=0.5,
                    )
                    
            elif isinstance(poly, Line):
                ax.plot(*poly.xy, color=color_line)

        
        # Plot the simulation line if it exists
        if hasattr(self, 'simulation_line') and self.simulation_line is not None:
            ax.plot(
                *self.simulation_line.xy, 
                color=color_simulation_line, 
                linewidth=linewidth_simulation,
                label="Simulation Line"
            )
            ax.legend()
            
            if len(self.line_segments) > 0:
                ax.legend()
        
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_title("PhotonicDevice with Simulation Line")
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        return fig, ax
    
    def plot_results(
            self,
            V_idx=None, 
            cmap = 'tab10',
            log_scale_carriers = True,
            plot_limits = True):
        """
        Plot simulation results in a 2x1 subplot layout.

        Args:
            V_idx (list, optional): Indices of voltages to plot. Defaults to first and last.
            cmap (str, optional): Colormap for voltage lines. Defaults to 'tab10'.
            log_scale_carriers (bool, optional): Whether to use log scale for carrier concentrations. Defaults to True.
            plot_limits (bool, optional): Whether to plot vertical lines at segment boundaries. Defaults to True.


        Returns:
            tuple: Figure and axes objects
        """
        if V_idx is None:
            V_idx = [0, len(self.V)-1]

        cmap = plt.get_cmap(cmap)
        norm = plt.Normalize(vmin=min(self.V), vmax=max(self.V))
        colors = [cmap(norm(self.V[v])) for v in V_idx]

        if V_idx == None:
            V_idx = [0,len(self.V)-1]
            
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6, 8),sharex=True)

        ax11 = ax1.twinx() #Ax to plot the vertical lines of the boundaries
        ax21 = ax2.twinx() #Ax to plot the vertical lines of the boundaries
        ax31 = ax3.twinx() #Ax to plot the vertical lines of the boundaries

        ax11.set_axis_off()
        ax21.set_axis_off()
        ax31.set_axis_off()

        for ax in [ax11, ax21, ax31]:
            ax.set_ylim(0,1)
            ax.set_yticks([])

        # ax2r = ax2.twinx()
        for i, v in enumerate(V_idx):
            ax1.plot(self.mesh, self.Ec[v], "-", color=colors[i])
            ax1.plot(self.mesh, self.Ev[v], "-", color=colors[i])
            # Plot quasi-Fermi levels
            ax1.plot(self.mesh, self.Efn[v], "-.", color=colors[i], linewidth=0.5)
            ax1.plot(self.mesh, self.Efp[v], "-.", color=colors[i], linewidth=0.5)
        # Configure first subplot
                        
            # ax2 = ax1.twinx()
            ax2.plot(self.mesh,  self.N[v],"-", color=colors[i])
            ax2.plot(self.mesh, self.P[v],"-.", color=colors[i])
            
            ax3.plot(self.mesh,self.Efield[v], color = colors[i])
            # ax3.set_ylim(-300,100)

        if plot_limits:
            for ax in [ax11, ax21, ax31]:
                for name, segment in self.line_segments.items():
                    x0, x1 = segment.xy[0]
                    y0, y1 = segment.xy[1]

                    ax.plot([y0, y0], [0,1], 'r-', label=f'Segment {name}', alpha = 0.2)
                    ax.plot([y1, y1], [0,1], 'r-', label=f'Segment {name}', alpha = 0.2)

            
        ax1.set_ylabel('Energy (eV)')
        ax1.grid(True, alpha=0.3)
        # ax1.legend(loc='best')

        ax2.set_ylabel(r"Carrier conc. ($cm^{-3}$)")
        ax2.grid(True, alpha=0.3)
        # ax2.legend(loc='best')
        if log_scale_carriers:
            ax2.set_yscale('log')
        ax2.set_ylim(1e15,1e19)

        #Add dummy lines so that we can put a legend for electron and hole
        l1, = ax2.plot([],[], 'k-', label='Electrons')
        l2, = ax2.plot([],[], 'k--', label='Holes')
        ax2.legend(handles=[l1,l2], loc='best')

        ax3.set_ylabel(r"Electric field (kV/cm)")
        ax3.grid(True, alpha=0.3)
        # ax3.legend(loc='best')

        from mpl_toolkits.axes_grid1 import make_axes_locatable

        norm_cbar = plt.Normalize(vmin=self.V[min(V_idx)], vmax=self.V[max(V_idx)])

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes('top', size='8%', pad=0.1)

        ax11.set_ylim(0,1.17) #A small compensation so that the lines are nice

        sm = plt.cm.ScalarMappable(norm=norm_cbar, cmap=cmap)
        sm.set_array([])  # required in older mpl versions

        cb = fig.colorbar(sm, cax=cax, orientation='horizontal')
        cb.ax.xaxis.set_ticks_position('top')
        cb.ax.xaxis.set_label_position('top')

        cb.set_label('Bias Voltage (V)')

        for v in [self.V[vi] for vi in V_idx]:
            cb.ax.axvline(v, color='k', lw=1.5)

    def plot_mesh(self):

        """

        Plot the mesh points along the simulation line.

        Returns:
            fig, ax: matplotlib figure and axis objects
        """
        
        fig, ax = plt.subplots(figsize=(10, 5))
        ax1 = ax.twinx()

        y0_min = self.simulation_line.xy[1][0]
        y0_max = self.simulation_line.xy[1][1]

        for name, segment in self.line_segments.items():
            x0, x1 = segment.xy[0]
            y0, y1 = segment.xy[1]

            ax.plot([y0, y0], [0,1], 'r-', label=f'Segment {name}')
            ax.plot([y1, y1], [0,1], 'r-', label=f'Segment {name}')

            ax.text(y0, 0.5, f'{name} start', rotation=90, verticalalignment='center')

        ax1.plot(np.linspace(y0_min, y0_max, len(self.mesh)), self.junction.mesh, marker= 'o', label='Mesh points')

        ax.set_ylabel('Mesh point')
        ax.set_xlabel('Real space (um)')

        return fig, ax