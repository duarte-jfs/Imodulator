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
        methods to import and configure external simulation tools: Lumerical
        (via lumapi) and nextnano (via nextnanopy). It manages system paths,
        software paths, and licensing information for those tools.

        Which file gets loaded:

        Step 1 -- pick the directory. If ``config_dir`` is not given, it
        defaults to the folder this Config.py lives in, i.e. the imodulator
        package folder itself (``src/imodulator/`` in a source checkout or
        editable install, or the installed package directory otherwise).
        That folder is where the package ships its own config.yaml and
        config_template.yaml.

        Step 2 -- pick the file inside that directory. The filename is always
        appended by this class, which is why ``config_dir`` is a *directory*
        and not a file: pass ``".../my_project"`` and NOT
        ``".../my_project/config.yaml"`` -- the latter would be looked up as
        ``".../my_project/config.yaml/config.yaml"``. Within the chosen
        directory:

        1. ``config.yaml`` is used if it exists;
        2. otherwise ``config_template.yaml`` is used, after printing a
           warning. This fallback keeps a fresh checkout importable, but it
           means a typo in ``config_dir`` does not raise immediately -- you
           silently get the template's placeholder paths instead. Check
           ``self.config_file`` if the loaded settings look wrong.

        So with no arguments you get ``<imodulator package>/config.yaml``,
        falling back to ``<imodulator package>/config_template.yaml``.

        A ``config_dir`` that does not exist, or that holds neither file,
        raises ``FileNotFoundError`` when the YAML is opened.

        Using the default instance:

        from imodulator.Config import config_instance as config
        nn = config.get_nextnanopy()
        lumapi = config.get_lumapi()

        Pointing at your own config directory:

        import imodulator.Config as config_module
        from imodulator.Config import Config

        config = Config("D:/Repos/my_project")
        config_module.config_instance = config

        Rebinding ``config_instance`` only affects code that reads it *after*
        this point. Modules that already ran
        ``from imodulator.Config import config_instance`` hold a reference to
        the old object, so in a notebook restart the kernel and do this before
        importing the rest of imodulator. See the module-level
        ``config_instance`` note below.

        Args:
            config_dir (str | Path | None): Directory containing 'config.yaml'.
                If omitted, the directory of this module is used, i.e. the
                config.yaml shipped inside the installed package.

        Attributes:
            config_dir (Path): Resolved directory containing the configuration file
            config_file (Path): Path to the file actually loaded -- either
                config.yaml or the config_template.yaml fallback
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
            # An all-comments or zero-byte file parses to None, not {}.
            self.config = yaml.safe_load(file) or {}

        self.lumapi = None
        self.nn = None
        self._lumapi_imported = False
        self._nn_imported = False

    def get_lumapi(self):
        if self._lumapi_imported:
            return self.lumapi

        self._lumapi_imported = True
        self.lumerical_api_path = self.config.get("lumerical_api", {}).get("path")

        if not self.lumerical_api_path:
            print(
                f"WARNING: no lumerical_api.path in {self.config_file}. "
                "Lumerical-backed simulators are unavailable on this machine."
            )
            return None

        try:
            spec = importlib.util.spec_from_file_location("lumapi", self.lumerical_api_path)
            lumapi = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(lumapi)
            self.lumapi = lumapi
            print("Successfully imported lumapi")
        except (ImportError, FileNotFoundError) as e:
            print(f"Failed to import lumapi: {e}")
            self.lumapi = None

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


#: Default, package-local Config, built at import time from the config.yaml
#: that sits next to this module. Importing imodulator therefore always reads
#: (or warns about) that file, even if you later build your own Config.
#:
#: To use a different config directory, rebind this attribute *before* the
#: rest of imodulator is imported -- see Config.__init__ for the caveats:
#:
#:     import imodulator.Config as config_module
#:     config_module.config_instance = config_module.Config("D:/Repos/my_project")
config_instance = Config()
