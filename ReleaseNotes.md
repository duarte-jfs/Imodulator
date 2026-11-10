# 1.1.0

## Breaking changes

- `ElectroOpticalSimulator` now takes `wavelength` in µm (default `1.55`) instead of nm, consistent with the rest of the package (#36).
- `transfer_results_to_device` of `ChargeSimulatorSolcore` and `ChargeSimulatorNN` no longer takes `xmin`/`xmax`; the mapping region is derived from the polygons. The new `polygons_to_map` argument lets you explicitly map results to more polygons than those used for the simulation, e.g. several mesas (#5).
- `ChargeSimulatorNN`: the third entry of `bias_start_stop_step` is now the number of bias points (nextnano's `steps` is points − 1). The default changed from `[0, 1, 1]` to `[0, 1, 2]`.
- `OpticalSimulatorMODE` now saves the simulation file by default (`save_sim=True`).
- Python 3.11 is required (`>=3.11,<3.12`).
- `femwell`, `solcore` and `nextnanopy` are now optional extras. Install what you need with `pip install imodulator[femwell]`, `[solcore]`, `[nextnanopy]` or `[all]`. Using a simulator without its backend raises an error naming the install command.

## Model changes (results differ from 1.0.x)

- Added intraband absorption for electron charges on the `InGaAsPElectroOpticalModel`. It now includes jumps from the Gamma to the X and L valleys.
- Added wavelength dependence to `InGaAsPElectroOpticalModel` and on the `ElectroOpticalSimulator`. The `InGaAsPElectroOpticalModel` can now be used for the simulation across the telecommunications bands.
- Re-fitted the Franz-Keldysh electro-absorption and Kerr coefficients of `InGaAsPElectroOpticalModel` so the model reproduces the Maat (2001) benchmark.
- Fixed the optical mode power normalization in `ElectroOpticalSimulator.calculate_EO_response` (spurious ½ factor removed, following Jia-Ming Liu, *Photonic Devices*). This changes the magnitude of the returned κ_ab(V).
- Fixed the crystal rotation matrix in `ElectroOpticalSimulator`: the x and z rotation angles were swapped.
- `RFSimulatorFEMWELL` and `OpticalSimulatorFEMWELL` now use a piecewise-constant (P0) basis for the permittivity, so material interfaces stay sharp instead of being smeared over the neighbouring elements.
- `ChargeSimulatorNN`: carrier densities were read without nextnano's 1e18 cm⁻³ unit factor.

## New features

- `ChargeSimulatorNN` supports graded (non-constant) quaternary alloys through the `alloy_type`, `alloy_x` and `alloy_y` keys of `SemiconductorPolygon` (#15). Requires nextnano 2026-02-05 or newer.
- Improved data transfer from `ChargeSimulatorSolcore` to `PhotonicDevice`.
- Improved data transfer of optical modes from `OpticalSimulatorFEMWELL` to `PhotonicDevice`.
- Added material dispersion compatibility in `PhotonicPolygon.optical_material` and subsequent compatibility with `OpticalSimulatorFEMWELL` and `OpticalSimulatorMODE` (new `wavelength` argument, in µm).
- `PhotonicDevice` now allows the user to specify `PhotonicPolygon`s to combine into a `MultiPolygon` for current line integral calculations, through the argument `line_integral_multi_polygons`. This is useful if a user has a complex metallic contact which might allow for different meshing resolutions for more efficient computation, while the current calculation remains over the whole metallic contact. This also fixes the current being integrated twice over the boundary between two touching metals (#31).
- Added confinement factor calculation to the `RFSimulator` (`RFSimulatorFEMWELL.confinement_factor(mode, polygons)`).
- `RFSimulatorFEMWELL.plot_mode` accepts `vector_density`, `axes` and `fig`; `plot_mesh` accepts `ax`.
- `plot_polygons` of `PhotonicDevice` and all simulators accepts `equal_aspect` for a 1:1 aspect ratio.
- `plot_results` of `ChargeSimulatorNN` and `ChargeSimulatorSolcore` accept `log_scale_carriers` and `plot_limits`; `ChargeSimulatorSolcore` adds `plot_design_doping` to overlay the design doping profile.
- `OpticalSimulatorMODE`:
	- the FDE region uses the `simulation_window` bounds when given, and `plot_mode` plots over the same bounds;
	- TE/TM propagation loss (dB/cm) is stored after `select_modes`;
	- `create_fde` can be called again on an existing simulation (#7);
	- empty polygons are skipped with a warning.
- `Config` accepts an optional `config_dir`, checks the nextnano paths on load and warns if they are wrong (`Config.validate_nextnano_paths()`, #10). A missing Lumerical API path now gives a warning instead of a crash.
- Utility functions have been added in `imodulator.utils` which include some functions and models that can generate broadband refractive index values for arbitrary doping, As concentration and wavelength for InGaAsP alloys lattice matched to InP.
- `SemiconductorPolygons`, `MetalPolygons` and `InsulatorPolygons` have a precision of 1nm set on `shapely`.

## Bug fixes

- Fixed the `epsilon_rf` generation. The real part was not respecting the polygon hierarchy.
- `OpticalSimulatorFEMWELL` skips polygons that are not in the mesh subdomains (#13).
- `OpticalSimulatorFEMWELL` and `RFSimulatorFEMWELL` no longer try to remove the same polygon several times when the hierarchy is passed more than once (#24).
- `OpticalSimulatorMODE`: the complex refractive index is now the square root of the complex permittivity (the real and imaginary parts were previously square-rooted separately). Automatic TE/TM selection no longer crashes when no indices are given.
- `ChargeSimulatorNN.load_output_data` no longer crashes on the `Structure` folder or on folders without `bandedges.dat` (#19, #26).

## Installation and documentation

- `femwell` is pinned to a git commit because the PyPI releases lack fixes Imodulator depends on.
- uv installs allow pre-releases (`solcore` depends on a `solsesame` pre-release); lockfile-based installs and the config file are documented in the installation guide.
- Added a Nix flake with a development shell.
- Benchmarks have been included in the documentation (Maat 2001 electro-optic benchmark).
- All tutorial notebooks were re-run against this release.
