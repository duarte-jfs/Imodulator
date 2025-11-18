# Welcome to the Imodulator's repository

**Imodulator** is an all-in-one tool for the simulation of electro-optic phase modulators.  
Simply define your geometry and materials, and then send the information to the various solvers available — including optical mode solver, RF mode solver with small-signal analysis, charge transport simulations, and electro-optic interaction simulations.

<p align="center">
  <img src="docs\architecture.png" width="600">
</p>

Check out the [docs](https://imodulator.readthedocs.io/en/latest/) to see how to install and use this package.

---

## Current Limitations

For the moment, the full functionality of the **Imodulator** package is limited to InGaAsP alloys lattice-matched to InP. However, when we consider the different parts of the simulator, we have various limitations that may or may not be relevant.

- **`OpticalSimulatorFEMWELL`**:  
  There is virtually no limitation here. As long as you provide a refractive index for each polygon, you're good to go.

- **`OpticalSimulatorMODE`**:  
  There is virtually no limitation here. As long as you provide a refractive index for each polygon, you're good to go.

- **`RFSimulatorFEMWELL`**:  
  There is also no limitation in this solver. You need only to input the material properties and it will work.

- **`ChargeSimulatorSolcore`**:  
  [`Solcore`](https://github.com/qpv-research-group/solcore5) has been developed with solar cells in mind, and we have found that the internal library of material parameters was limiting for the purpose of this package.  
  Therefore, we have made a connection between `Solcore` and `openbandparams` so that we can use arbitrary III–V alloys (excluding strain effects) in solving the Poisson–drift–diffusion equations.  
  The limitation here is that we must work with III–V alloys only.  
  Furthermore, the mobility values are calculated through `Solcore` via [mobility_solcore](https://github.com/qpv-research-group/solcore5/blob/develop/solcore/material_data/mobility.py), and we are therefore limited to:

  - InGaAs  
  - InGaP  
  - AlGaAs  
  - InAlAs  
  - InGaAsP  

- **`ChargeSimulatorNN`**:  
  There isn’t really a limitation here. We only need to provide materials supported by [NextNano](https://www.nextnano.com/error_pages/404.php).

- **`ElectroOpticalModels`**:  
  We have only included a model compatible with electro-optical effects that take place in InGaAsP alloys lattice-matched to InP. However, the software has been written to allow for any model, as long as we provide a  
  
  ![formula](https://latex.codecogs.com/svg.image?\Delta\bar{\epsilon}(V,E_c,E_v,E_{fp},E_{fv},\mu_n,\mu_p,\vec{E},N,P))  
  
  function.

---

## Where You Can Contribute

- Generalization of [openbandparams](https://github.com/duarte-jfs/openbandparams) to include other semiconductor compounds such as Si and SiGe.  
- Generalization of [mobility_solcore](https://github.com/qpv-research-group/solcore5/blob/develop/solcore/material_data/mobility.py) to include other mobility models explored in [Sotoodeh *et al.*, 2000](https://doi.org/10.1109/16.826799).  
  Alternatively, one could explore the inclusion of those models directly in `openbandparams`.
- Include more electro-optic models.  
- Include surface impedance boundary conditions in the RF mode solver. 
- Include a 2D PDD solver based on [sesame](https://sesame.readthedocs.io/en/latest/) integrated with `openbandparams`
- Improve the documentation (**help wanted!**).

If you have any questions please reach out to the **Discussions** tab and we can brainstorm some ideas.

## Acknowledgements

This work was funded by the European Union through the QuGANTIC project and the Dutch National Growth Fund and PhotonDelta. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Innovation Council. Neither the European Union nor the granting authority can be held responsible for them.

## Contributors

- Duarte Silva (Eindhoven University of Technology)
- Ali Kaan Sünnetçioğlu (Eindhoven University of Technology)

## Was this package useful for your work?

The DOI tracking for this project is handled by Zenodo, and a single DOI shall be generated per official release. 

|  Version | DOI  |
|----------|------|
|  all     | https://doi.org/10.5281/zenodo.17633733 |
|  v1.0.2  | https://doi.org/10.5281/zenodo.17639012 |



#### Which DOI should I use in citations?

You should normally always use the DOI for the specific version of your record in citations. This is to ensure that other researchers can access the exact research artefact you used for reproducibility. By default, we use the specific version to generate citations.

You can use the Concept DOI representing all versions in citations when it is desirable to cite an evolving research artifact, without being specific about the version.

When translating to a bibtex entry, please consider using the following formatting to display the author names properly:

```
@MISC{,
  title     = "Imodulator",
  author    = "{Duarte J.F. da Silva, Ali Kaan Sünnetçioğlu}",
  publisher = "Zenodo",
  year      =  2025,
  url       = {https://github.com/duarte-jfs/Imodulator},
  doi       = {https://doi.org/10.5281/zenodo.17633733}
}
```