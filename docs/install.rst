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

Reproducible environment with Nix (flake)
------------------------------------------

The repository also ships a ``flake.nix`` that builds the same pinned
environment (femwell + solcore extras) as a self-contained Nix store path, with
the system libraries the wheels need baked in — nothing is installed
system-wide. The same commands work on **Linux** (``x86_64-linux``) and
**macOS / Apple Silicon** (``aarch64-darwin``):

.. code-block:: console

   $ cd Imodulator
   $ nix develop .          # the environment on your PATH, plus uv

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
In that case we suggest you checkout to development branch for the latest 
changes and fixes. 

.. code-block:: console

   (.venv) $ git clone https://github.com/duarte-jfs/Imodulator.git
   (.venv) $ cd Imodulator
   (.venv) $ pip install -e .

Configuring nextnano and Lumerical paths
----------------------------------------

``ChargeSimulatorNN`` (nextnano) and ``OpticalSimulatorMODE`` (Lumerical MODE)
drive third-party software that Imodulator cannot install for you. You have to
tell it where those tools live by writing a ``config.yaml``. Every other
simulator works without it.

Step 1 — find where the file goes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``config.yaml`` is read from the *installed package directory*, i.e. next to
``Config.py`` — not from your working directory or your home folder. Print that
directory with:

.. code-block:: console

   $ python -c "import imodulator, pathlib; print(pathlib.Path(imodulator.__file__).parent)"

If you installed from a clone (``uv sync`` or ``pip install -e .``) this is just
``src/imodulator/`` in the repository.

Step 2 — copy the template
^^^^^^^^^^^^^^^^^^^^^^^^^^

The package ships a ``config_template.yaml`` in that same directory. Copy it,
dropping the ``_template`` suffix:

.. code-block:: console

   $ cd "$(python -c 'import imodulator, pathlib; print(pathlib.Path(imodulator.__file__).parent)')"
   $ cp config_template.yaml config.yaml

On Windows (PowerShell):

.. code-block:: powershell

   > cd (python -c "import imodulator, pathlib; print(pathlib.Path(imodulator.__file__).parent)")
   > Copy-Item config_template.yaml config.yaml

Step 3 — fill in the paths
^^^^^^^^^^^^^^^^^^^^^^^^^^

Open ``config.yaml`` and replace the placeholders. A filled-in file looks like
this — copy it and edit the values to match your machine:

.. code-block:: yaml

   lumerical_api:
     path: "C:/Program Files/Lumerical/v242/api/python/lumapi.py"

   nextnano:
     nextnano++:
       exe: "C:/Program Files/nextnano/2024_12_16/nextnano++/bin/nextnano++_Intel_64bit.exe"
       license: "C:/Users/<you>/Documents/nextnano/License/License_nnp.lic"
       database: "C:/Program Files/nextnano/2024_12_16/nextnano++/Syntax/database_nnp.in"
       output: "C:/Users/<you>/Documents/nextnano/Output"

What each key means:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Value
   * - ``lumerical_api.path``
     - Full path **to the** ``lumapi.py`` **file itself**, not to the folder
       containing it. It is loaded directly as a module file, so a directory
       path fails.
   * - ``nextnano.nextnano++.exe``
     - The nextnano++ executable (pick the ``Intel_64bit`` build unless you
       have reason not to).
   * - ``nextnano.nextnano++.license``
     - Your ``.lic`` license file, usually under
       ``Documents/nextnano/License``.
   * - ``nextnano.nextnano++.database``
     - The ``database_nnp.in`` material database shipped with nextnano.
   * - ``nextnano.nextnano++.output``
     - A folder where simulation results are written. Imodulator creates an
       ``nn_output`` subfolder inside it, so point this at a parent directory
       you own — it does not need to exist beforehand.

Use forward slashes ``/`` even on Windows, and quote every path — backslashes
inside a quoted YAML string are escape characters and will corrupt the path.

You only need to fill in the section for the tool you actually use; the
Lumerical and nextnano blocks are read independently, and only at the moment
that tool is first requested. Keep both top-level keys present, though — remove
the ``lumerical_api`` key entirely and importing ``lumapi`` raises a
``KeyError`` instead of a helpful message.

Step 4 — check it worked
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from imodulator.Config import config_instance as config

   config.get_lumapi()       # -> "Successfully imported lumapi"
   config.get_nextnanopy()   # -> "Successfully imported nextnanopy"
                             #    "Successfully configured nextnano++ settings"

If instead you see::

   WARNING: Configuration file not found: .../imodulator/config.yaml. Using template file instead.

the file is missing or misnamed — it fell back to the template, whose
placeholder paths point nowhere, and the import will fail afterwards. Re-check
Step 1 and Step 2, and make sure the name is exactly ``config.yaml``.

.. note::

   ``config.yaml`` is git-ignored, so your local paths and license location are
   never committed. The flip side: because it lives inside the package
   directory, reinstalling or upgrading Imodulator into a fresh environment
   does not carry it over. Keep a copy somewhere outside the package if you
   rebuild environments often.

