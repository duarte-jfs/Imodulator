# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

**Imodulator** is a Python package for simulating electro-optic phase modulators. It integrates multiple FEM-based solvers (optical modes, RF modes, charge transport, electro-optic interaction) for arbitrary cross-sections, currently with full support for InGaAsP alloys lattice-matched to InP.

## Build & Install

```bash
# Install in editable mode (development)
pip install -e .

# Build distribution archives
python -m build
```

There are no automated tests in this repository currently.

## Documentation

```bash
# Build HTML docs (run from the docs/ directory)
python -m sphinx-build -M html . ./_build -a

# Doc dependencies
pip install sphinxcontrib-bibtex sphinx_autodoc_typehints sphinx-rtd-theme sphinx nbsphinx
```

## Architecture Overview

The simulation pipeline flows from geometry definition → device assembly → individual simulators → electro-optic integration:

### 1. Geometry: `PhotonicPolygon.py`
Three polygon types hold per-region material and mesh properties:
- `SemiconductorPolygon` — semiconductor regions; carries `rf_eps`, `electro_optic_module`, doping/alloy kwargs for charge transport, and separate mesh settings per solver (`charge_mesh_settings`, `eo_mesh_settings`, `rf_mesh_settings`, `optical_mesh_settings`)
- `MetalPolygon` — conductor regions
- `InsulatorPolygon` — dielectric regions

### 2. Device: `PhotonicDevice.py`
`PhotonicDevice(photo_polygons: list)` assembles all polygons into a single device. **Polygon order matters** — it determines hierarchy in the meshing algorithm (first = highest priority). It auto-generates a background bounding box if no `"background"` polygon is provided. Results from each simulator are stored here (`device.mode`, `device.charge`).

### 3. Simulators (each takes a `PhotonicDevice`)
| Class | File | Backend | Purpose |
|---|---|---|---|
| `OpticalSimulatorFEMWELL` | `OpticalSimulator.py` | femwell + scikit-fem | FEM optical mode solving |
| `OpticalSimulatorMODE` | `OpticalSimulator.py` | Lumerical MODE (via lumapi) | Commercial optical solver |
| `RFSimulatorFEMWELL` | `RFSimulator.py` | femwell + scikit-fem | RF mode + small-signal analysis |
| `ChargeSimulatorSolcore` | `ChargeSimulator.py` | solcore + openbandparams | 1D Poisson-drift-diffusion |
| `ChargeSimulatorNN` | `ChargeSimulator.py` | NextNano++ | 2D/3D charge transport |
| `ElectroOpticalSimulator` | `ElectroOpticalSimulator.py` | femwell + scikit-fem | EO overlap integral |

All FEM simulators use **femwell** for mesh generation (`mesh_from_OrderedDict`) and **scikit-fem** for solving. Meshes are built from `OrderedDict` of Shapely geometries.

### 4. Electro-Optic Models: `ElectroOpticalModel.py`
Contains material models (e.g., `InGaAsPElectroOpticalModel`) that compute `Δε̄(V, Ec, Ev, Efp, Efv, μn, μp, E, N, P)`. These are assigned to `SemiconductorPolygon.electro_optic_module` and consumed by `ElectroOpticalSimulator`. Currently only InGaAsP models are implemented, but the interface is extensible.

### 5. Configuration: `Config.py`
Loads `config.yaml` (or falls back to `config_template.yaml`) to provide paths to external tools. A module-level singleton `config_instance` is used by simulators at import time:
```python
from imodulator.Config import config_instance as config
lumapi = config.get_lumapi()        # Lumerical API
nn = config.get_nextnanopy()        # nextnanopy
```
To configure external tools, copy `src/imodulator/config_template.yaml` to `src/imodulator/config.yaml` and fill in paths. **`config.yaml` is not committed.**

## Key Dependencies

- **femwell** (0.1.11) — mesh generation and waveguide mode computation
- **scikit-fem** — FEM basis and assembly
- **shapely** — all geometry operations (polygons, intersections, clipping)
- **solcore** — Poisson-drift-diffusion solver and mobility models (limited to InGaAs, InGaP, AlGaAs, InAlAs, InGaAsP)
- **openbandparams** — III-V alloy band parameters for arbitrary compositions
- **nextnanopy** — Python interface to NextNano++ for charge simulations
- **pint** — physical units throughout the codebase

## Branching Model

Uses [Gitflow](http://nvie.com/posts/a-successful-git-branching-model/): feature branches merge into `develop`, releases cut from `develop` to `main`. When releasing: bump version in both `pyproject.toml` and `docs/conf.py`.
