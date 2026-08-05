Install
=======

.. _installation:


To use Imodulator, simply install via pip:

.. code-block:: console

   (.venv) $ pip install imodulator

Using uv (recommended)
----------------------

`uv <https://docs.astral.sh/uv/>`_ creates the virtual environment and resolves
dependencies in one step, and is considerably faster than ``pip``. Install it
once (see the `uv installation guide <https://docs.astral.sh/uv/getting-started/installation/>`_),
then:

.. code-block:: console

   $ uv venv
   $ uv pip install imodulator

``uv venv`` creates ``.venv`` in the current directory. Activate it with
``source .venv/bin/activate`` (``.venv\Scripts\activate`` on Windows), or skip
activation entirely by prefixing commands with ``uv run``.

Add the extras you need in the same command:

.. code-block:: console

   $ uv pip install imodulator[femwell]
   $ uv pip install imodulator[femwell,solcore]
   $ uv pip install imodulator[all]

Reproducible installs with uv.lock
----------------------------------

The repository ships a ``uv.lock`` file. It records the exact version of every
package in the dependency tree — direct and transitive, across all extras — so
that everyone working from a clone builds the *same* environment instead of
whatever happens to be newest on PyPI that day.

Use ``uv sync``, which installs strictly from the lockfile:

.. code-block:: console

   $ git clone https://github.com/duarte-jfs/Imodulator.git
   $ cd Imodulator
   $ uv sync --all-extras

``uv sync`` creates ``.venv`` if needed, installs the project in editable mode,
and removes anything not in the lockfile, so the environment matches exactly.
Pick a subset with ``--extra`` instead of ``--all-extras``:

.. code-block:: console

   $ uv sync --extra femwell
   $ uv sync --extra femwell --extra solcore
   $ uv sync --extra dev          # formatting toolchain (ruff)

Run commands inside that environment without activating it:

.. code-block:: console

   $ uv run python my_simulation.py

.. important::

   ``uv pip install`` does **not** read ``uv.lock`` — it resolves fresh every
   time. Only ``uv sync`` and ``uv run`` use the lockfile. If you want the
   pinned, reproducible environment, use ``uv sync``.

To verify your environment matches the lockfile without changing anything —
useful in CI, where it should fail rather than silently re-resolve:

.. code-block:: console

   $ uv sync --all-extras --locked

After editing dependencies in ``pyproject.toml``, regenerate the lockfile and
commit it along with your change:

.. code-block:: console

   $ uv lock

Note that ``uv.lock`` only governs installs from a clone. Users installing the
published package with ``pip install imodulator`` resolve dependencies normally
against the constraints in ``pyproject.toml``.

Optional extras
---------------

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

   ``solcore`` depends on ``solsesame==2.1a1``, which is a **pre-release**.
   ``pip`` installs it without complaint because the version is pinned exactly,
   and current ``uv`` resolves it too. To keep that from depending on the
   resolver's default behaviour, this project sets ``prerelease = "allow"``
   under ``[tool.uv]`` in ``pyproject.toml``, so no extra flags are needed when
   working from a clone.

   If you do hit a pre-release resolver error — an older ``uv``, or installing
   the published package outside a clone — pass the flag explicitly:

   .. code-block:: console

      (.venv) $ uv pip install 'imodulator[solcore]' --prerelease=allow

Installing from source
----------------------

Alternatively you can clone the repository and install it locally:

.. code-block:: console

   (.venv) $ git clone https://github.com/duarte-jfs/Imodulator.git
   (.venv) $ cd Imodulator
   (.venv) $ pip install .

For development purposes, you can install the package in "editable" mode:

.. code-block:: console

   (.venv) $ git clone https://github.com/duarte-jfs/Imodulator.git
   (.venv) $ cd Imodulator
   (.venv) $ pip install -e .

Configuring nextnano and Lumerical paths
----------------------------------------

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

