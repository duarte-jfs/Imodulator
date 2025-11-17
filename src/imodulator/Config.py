import sys
from pathlib import Path
import yaml
import importlib.util
import os

class Config:
    def __init__(self):
        """
        Configuration manager for the imodulator package.

        This class handles loading configuration settings from a YAML file and provides
        methods to import and configure external simulation tools like Lumerical, 
        nextnano, and InGaAsP models.

        The configuration file should be located in the same directory as this module
        and named 'config.yaml'. It manages system paths, software paths, and licensing
        information for various photonic simulation tools.
        
        Example: 

        from imodulator.Config import config_instance as config
        nn = config.get_nextnanopy()
        lumapi = config.get_lumapi()
        InGaAsP_models = config.get_ingaasp_models()

        Attributes:
            config_dir (Path): Directory containing the configuration file
            config_file (Path): Path to the config.yaml file
            config (dict): Loaded configuration data from YAML file
            lumapi: Lumerical API module (if successfully imported)
            nn: nextnanopy module (if successfully imported)
        """
        # Get the directory where this config.py file is located
        self.config_dir = Path(__file__).resolve().parent

        # Load the configuration from YAML file
        self.config_file = self.config_dir / 'config.yaml'

        if not self.config_file.exists():
            print(f"WARNING: Configuration file not found: {self.config_file}. Using template file instead.")

            self.config_file = self.config_dir / 'config_template.yaml'

        with open(self.config_file, 'r') as file:
            self.config = yaml.safe_load(file)

        self.lumapi = None
        self.nn = None
        self._lumapi_imported = False
        self._nn_imported = False

    def get_lumapi(self):
        if self._lumapi_imported:
            return self.lumapi

        self.lumerical_api_path = self.config['lumerical_api']['path']
        try:
            spec = importlib.util.spec_from_file_location('lumapi', self.lumerical_api_path)
            lumapi = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(lumapi)
            self.lumapi = lumapi
            print("Successfully imported lumapi")
        except (ImportError, FileNotFoundError) as e:
            print(f"Failed to import lumapi: {e}")
            self.lumapi = None

        self._lumapi_imported = True
        return self.lumapi
    
    def get_nextnanopy(self):
        if self._nn_imported:
            return self.nn

        try:
            import nextnanopy as nn
            self.nn = nn  # Store the module in an instance attribute
            print("Successfully imported nextnanopy")
            
            if 'nextnano' in self.config and 'nextnano++' in self.config['nextnano']:
                nnp_config = self.config['nextnano']['nextnano++']
                OUTPUT_FOLDER = os.path.join(nnp_config['output'],'nn_output')
                nn.config.set('nextnano++', 'exe', nnp_config['exe'])
                nn.config.set('nextnano++', 'license', nnp_config['license'])
                nn.config.set('nextnano++', 'database', nnp_config['database'])
                nn.config.set('nextnano++','outputdirectory',OUTPUT_FOLDER)
                print("Successfully configured nextnano++ settings")
            
            if 'nextnano' in self.config and 'nextnano3' in self.config['nextnano']:
                nnp_config = self.config['nextnano']['nextnano3']
                nn.config.set('nextnano3', 'exe', nnp_config['exe'])
                nn.config.set('nextnano3', 'license', nnp_config['license'])
                nn.config.set('nextnano3', 'database', nnp_config['database'])
                print("Successfully configured nextnano3 settings")  # Fixed the message

        except ImportError as e:
            print(f"Failed to import nextnanopy: {e}")
            self.nn = None  

        self._nn_imported = True
        return self.nn
    
config_instance = Config()
# Create an instance of the Config class
