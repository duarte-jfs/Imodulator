# 1.1.0

- Added intraband absorption for electron charges on the `InGaAsPElectroOpticalModel`. It now includes jumps from the Gamma to the X and L valleys.
- `SemiconductorPolygons`, `MetalPolygons` and `InsulatorPolygons` have a precision of 1nm set on `shapely`.
- `ChargeSimulatorSolcore` now allows you to explicitly map results to more polygons than those it used for simulation. 
- Improved data transfer from `ChargeSimulatorSolcore` to `PhotonicDevice`. 
- Improved data transfer of optical modes from `OpticalSimulatorFEMWELL` to `PhotonicDevice`.
- Added wavelength dependence to `InGaAsPElectroOpticalModel` and on the `ElectroOpticalSimulator`. The `InGaAsPElectroOpticalModel` can now be used for the simulation accross the telecomunications bands. 
- Fixed the `epsilon_rf` generation. The real part was not respecting the polygon hierarchy.+
- Added material dispersion compatibility in `PhotonicPolygon.optical_material` and subsequent compatibility with `OpticalSimulatorFEMWELL`.
- `PhotonicDevice` now allows the user to specify `PhotonicPolygon`s to combine into a `MultiPolygon` for current line integral calculations, through the argument `line_integral_multi_polygons`. This is useful if a user has a complex metalic contact which might allow for different meshing resolutions for more efficient computation, while the current calculation remains over the whole metalic contact.