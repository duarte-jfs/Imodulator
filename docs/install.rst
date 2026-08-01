Install
=======

.. _installation:


To use Imodulator, simply install via pip:

.. code-block:: console

   (.venv) $ pip install imodulator

The simulation backends are optional dependencies, installed as *extras* so you
only pull in what you need:

.. list-table::
   :header-rows: 1

   * - Extra
     - Enables
   * - ``imodulator[femwell]``
     - Optical and RF mode solving (``OpticalSimulatorFEMWELL``, ``RFSimulatorFEMWELL``, ``ElectroOpticalSimulator``)
   * - ``imodulator[solcore]``
     - Solcore-based charge transport (``ChargeSimulatorSolcore``)
   * - ``imodulator[nextnanopy]``
     - nextnano-based charge transport (``ChargeSimulatorNN``)
   * - ``imodulator[all]``
     - All of the above

Extras combine, so you can request any subset in one command:

.. code-block:: console

   (.venv) $ pip install imodulator[femwell]
   (.venv) $ pip install imodulator[nextnanopy,femwell]
   (.venv) $ pip install imodulator[all]

Instantiating a simulator whose extra is not installed raises a clear error
telling you which ``pip install`` command to run.

.. note::

   ``solcore`` depends on ``solsesame==2.1a1``, which is a **pre-release**. Plain
   ``pip`` (including inside a conda environment) installs it automatically
   because the version is pinned exactly, but ``uv`` rejects pre-releases by
   default. If you install the ``solcore`` extra with ``uv`` and hit a resolver
   error, allow pre-releases explicitly:

   .. code-block:: console

      (.venv) $ uv pip install -e '.[solcore]' --prerelease=allow

Alternatively you can clone the repository and install it locally:

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

