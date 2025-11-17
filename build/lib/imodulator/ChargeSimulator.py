from __future__ import annotations

from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from shapely.geometry import (
    Point,
    Polygon,
    MultiPolygon,
    LineString,
    MultiLineString,
    LinearRing,
)
from scipy.interpolate import RegularGridInterpolator

from shapely.ops import clip_by_rect, linemerge, unary_union
import os

import numpy as np
import pandas as pd
import copy


from collections import OrderedDict

from imodulator import PhotonicDevice
from imodulator.PhotonicPolygon import (
    SemiconductorPolygon,
    MetalPolygon,
    InsulatorPolygon,
)

PhotonicPolygon = SemiconductorPolygon | MetalPolygon | InsulatorPolygon
Line = LineString | MultiLineString | LinearRing

# from imodulator.ElectroOpticalModel import InGaAsPElectroOpticalModel
##Configured imports
from imodulator.Config import config_instance
# Get access to imported modules
nn = config_instance.get_nextnanopy()
# lumapi = config_instance.get_lumapi()
InGaAsP_models = config_instance.get_ingaasp_models()

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
                    alloy_y = self.photonicdevice.entities[segment_name]
                    alloy_x = self.photonicdevice.entities[segment_name]
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
                    name = "Ga(x)In(1-x)As(y)P(1-y)"    
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
            ax1.plot(self.grid, self.Ec[v], "-", color=colors[i], label=r"$E_c(eV) @ V={{{:.1f}}} V)$".format(v))
            ax1.plot(self.grid, self.Ev[v], "-", color=colors[i])
            # Plot quasi-Fermi levels
            ax1.plot(self.grid, self.Efn[v], "-.", color=colors[i], linewidth=0.5)
            ax1.plot(self.grid, self.Efp[v], "-.", color=colors[i], linewidth=0.5)
        # Configure first subplot
                        
            # ax2 = ax1.twinx()
            ax2.plot(self.grid,  self.N[v]*1e18,"-", color=colors[i], label=r"e conc. @ V={{{:.1f}}} V)$".format(v))
            ax2.plot(self.grid, -self.P[v]*1e18,"-.", color=colors[i], label=r"-h conc. @ V={{{:.1f}}} V)$".format(v))
            
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