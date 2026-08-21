import sys
from pathlib import Path
import yaml
import importlib.util
import os


class Config:
    def __init__(self, config_dir=None):
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

        Args:
            config_dir (str | Path | None): Directory containing 'config.yaml'.
                If omitted, the directory of this module is used.

        Attributes:
            config_dir (Path): Directory containing the configuration file
            config_file (Path): Path to the config.yaml file
            config (dict): Loaded configuration data from YAML file
            lumapi: Lumerical API module (if successfully imported)
            nn: nextnanopy module (if successfully imported)
        """
        # If path is not defined get the directory where this config.py file is located
        if config_dir is None:
            self.config_dir = Path(__file__).resolve().parent
        else:
            self.config_dir = Path(config_dir).resolve()

        # Load the configuration from YAML file
        self.config_file = self.config_dir / "config.yaml"

        if not self.config_file.exists():
            print(
                f"WARNING: Configuration file not found: {self.config_file}. Using template file instead."
            )

            self.config_file = self.config_dir / "config_template.yaml"

        with open(self.config_file, "r") as file:
            self.config = yaml.safe_load(file)

        self.lumapi = None
        self.nn = None
        self._lumapi_imported = False
        self._nn_imported = False

    def get_lumapi(self):
        if self._lumapi_imported:
            return self.lumapi

        self.lumerical_api_path = self.config["lumerical_api"]["path"]
        try:
            spec = importlib.util.spec_from_file_location("lumapi", self.lumerical_api_path)
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

            if "nextnano" in self.config and "nextnano++" in self.config["nextnano"]:
                nnp_config = self.config["nextnano"]["nextnano++"]
                OUTPUT_FOLDER = os.path.join(nnp_config["output"], "nn_output")
                nn.config.set("nextnano++", "exe", nnp_config["exe"])
                nn.config.set("nextnano++", "license", nnp_config["license"])
                nn.config.set("nextnano++", "database", nnp_config["database"])
                nn.config.set("nextnano++", "outputdirectory", OUTPUT_FOLDER)
                print("Successfully configured nextnano++ settings")

            if "nextnano" in self.config and "nextnano3" in self.config["nextnano"]:
                nnp_config = self.config["nextnano"]["nextnano3"]
                nn.config.set("nextnano3", "exe", nnp_config["exe"])
                nn.config.set("nextnano3", "license", nnp_config["license"])
                nn.config.set("nextnano3", "database", nnp_config["database"])
                print("Successfully configured nextnano3 settings")  # Fixed the message

            # Confirm the configured paths actually exist before they are used
            # by a simulation, where a wrong path would otherwise fail cryptically.
            self.validate_nextnano_paths()

        except ImportError as e:
            print(f"Failed to import nextnanopy: {e}")
            self.nn = None

        self._nn_imported = True
        return self.nn

    def validate_nextnano_paths(self, verbose=True):
        """
        Verify that the nextnano paths tagged on the nextnanopy API point to
        files/directories that actually exist.

        This is a "test step" for the configuration: after ``get_nextnanopy`` has
        imported nextnanopy and tagged the API via ``nn.config.set(...)``, this
        method reads those values back from ``nn.config`` and checks each one on
        disk. It catches placeholder template paths (e.g. ``path/to/nn++.exe``)
        and typos before they cause an obscure failure inside a simulation.

        Args:
            verbose (bool): If True (default) print the outcome for each path.

        Returns:
            dict: Mapping of ``"<variant>.<key>"`` to a tuple
            ``(path, is_valid)`` for every path that was checked. The output
            directory is created if it does not exist, so it is always valid.
        """
        results = {}

        if self.nn is None:
            if verbose:
                print("Cannot validate nextnano paths: nextnanopy is not imported.")
            return results

        # Only validate the variants imodulator actually configured in
        # config.yaml (rather than every variant present in nextnanopy's global
        # config), so unused variants don't produce spurious warnings.
        configured_variants = list(self.config.get("nextnano", {}).keys())

        # exe/license/database must exist on disk; the output directory is
        # created on demand rather than required to pre-exist.
        keys = ("exe", "license", "database")

        all_valid = True
        for variant in configured_variants:
            if variant not in self.nn.config.config:
                continue

            for key in keys:
                path = self.nn.config.config[variant].get(key)
                is_valid = bool(path) and os.path.exists(path)
                results[f"{variant}.{key}"] = (path, is_valid)
                if not is_valid:
                    all_valid = False
                if verbose:
                    status = "OK" if is_valid else "MISSING"
                    print(f"  [{status}] {variant} {key}: {path}")

            output_dir = self.nn.config.config[variant].get("outputdirectory")
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                results[f"{variant}.outputdirectory"] = (output_dir, True)
                if verbose:
                    print(f"  [OK] {variant} outputdirectory: {output_dir}")

        if verbose:
            if all_valid:
                print("nextnano configuration test passed: all paths are valid.")
            else:
                print(
                    "WARNING: some nextnano paths are invalid. "
                    "Update config.yaml with the correct paths."
                )

        return results


config_instance = Config()
# Create an instance of the Config class
