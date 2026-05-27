# 1.1.0

- Added intraband absorption for electron charges on the `InGaAsPElectroOpticalModel`. It now includes jumps from the Gamma to the X and L valleys.
- `SemiconductorPolygons`, `MetalPolygons` and `InsulatorPolygons` have a precision of 1nm set on `shapely`.
- `ChargeSimulatorSolcore` now allows you to explicitly map results to more polygons than those it used for simulation. 
- Improved data transfer from `ChargeSimulatorSolcore` to `PhotonicDevice`. 
- Improved data transfer of optical modes from `OpticalSimulatorFEMWELL` to `PhotonicDevice`.
- Added wavelength dependence to `InGaAsPElectroOpticalModel` and on the `ElectroOpticalSimulator`. The `InGaAsPElectroOpticalModel` can now be used for the simulation accross the telecomunications bands. 
- Fixed the `epsilon_rf` generation. The real part was not respecting the polygon hierarchy.