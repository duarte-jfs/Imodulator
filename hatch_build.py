from hatchling.builders.hooks.plugin.interface import BuildHookInterface
import subprocess
import sys

class CustomBuildHook(BuildHookInterface):
    def initialize(self, version, build_data):
        print("Installing proj1 from Git (no dependencies)...")
        print("Python executable:", sys.executable)
        subprocess.check_call([sys.executable, '-m',
            "pip", "install", "--no-deps",
            "git+https://github.com/HelgeGehring/femwell.git@36e2ff1d8507e3839f29b5f14c298a091b463c49"
        ])