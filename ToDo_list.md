## Tasks:

- Talk with Xaveer for help in losses per polygon.

- +A config file for file paths and api's
- +Look at \_\_getter and setter fuctions for sweep capabilities

### Duarte

- Currently there are some issues with femwell. First there is no mesh refinement that preserves named boundaries. This needs to be developed. The trick is to find all the facets that are split into 2 and then insert the new facet in the correct place so as to preserve the orientation of the line. Another issue is that the orientation of lines provided with the shapely.Polygon.exterior is not correctly inserted https://github.com/kinnala/scikit-fem/issues/1150
- Make it an option for the RF simulator to self.get_epsilon_rf() to either use charge transport sim data or not.

#### Charge

- What is the most versatile architecture to interface nextnano?
- For large signal simulation it will be needed to simulate both arms of a CPS;
- For metals you can create a polygon (https://www.nextnano.com/documentation/tools/nextnanoplus/keywords/structure/region_shapes.html#polygon) and define it as a contact;
- 2D simulation of the entire structure would be easiest. Allow for 1D in NextNano would be more efficient. But how to do it properly. We can do arbitrary geometries in nextnano by importing material properties in a 2D grid (https://www.nextnano.com/documentation/tools/nextnanoplus/tutorials/basics/structure/import-dat_1D_2D_3D.html);

#### Kaan

-Note somewhere the background polygon should cover the total area of the simulation
+ Mesh per solid
- Interpolators output to device

#### Comments

The mobility model used by NextNano is not to be trusted blindly. It returns, for InGaAs n doped with 1e19 a mobility of ~7400 cm^2/(Vs), yet a low field empirical model given by [1] establishes that such a value for InP lattice matched InGaAs will only happen at room temperature for very low doping. It also drops quite steeply for heavy doping. Their model predicts 1815 cm^2/(Vs) for 1e19 n doped InGaAs, which is far lower.

[1] - Sotoodeh, M., A. H. Khalid, and A. A. Rezazadeh. 2000. “Empirical Low-Field Mobility Model for III–V Compounds Applicable in Device Simulation Codes.” Journal of Applied Physics 87 (6): 2890–2900. https://doi.org/10.1063/1.372274.
