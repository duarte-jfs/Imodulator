Install
=======

.. _installation:


To use Imodulator, simply install via pip:

.. code-block:: console

   (.venv) $ pip install imodulator
   (.venv) $ pip install --no-deps git+https://github.com/HelgeGehring/femwell.git@36e2ff1d8507e3839f29b5f14c298a091b463c49

The second line has been added as a temporary measure. New changes have been made in FEMWELL on the mesh refinement which have not been made into a release yet. Once a new FEMWELL release is made, it shall be added to the dependencies of Imodulator and the second line can be omitted.
Alternatively you can close the repository and install it locally:

.. code-block:: console

   (.venv) $ git clone https://github.com/yourusername/imodulator.git
   (.venv) $ cd imodulator
   (.venv) $ pip install .

For development purposes, you can install the package in "editable" mode:

.. code-block:: console

   (.venv) $ git clone https://github.com/yourusername/imodulator.git
   (.venv) $ cd imodulator
   (.venv) $ pip install -e .

If you plan to use the `ChargeSimmulatorNN` or the `OpticalSimulatorMODE`, you must provide some more information relative to the licences and the paths to the api's. To do so, you must go to `imodulator/src/config_template.yaml` and create a copy named `config.yaml` where you fill in the required information. Namely:

.. code-block:: yaml

   lumerical_api: 
      path: "path to lumapi python package lumapi.py"

   nextnano:
      nextnano++:
         exe: "path to nextnano++ executable"
         license: "path to nextnano++ license file"
         database: "path to nextnano++ database file"
         output: "path to nextnano++ output folder"

