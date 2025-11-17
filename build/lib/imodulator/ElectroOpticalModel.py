from __future__ import annotations
import numpy as np
from pint import UnitRegistry, Quantity
from scipy.interpolate import interp1d, make_interp_spline
from scipy.special import exp1, expi


def _get_n(E, y=0):
    """
    This file contains the code to generate the real refractive index as well as the
    absorption coefficient according to [1].

    [1] -Sten Seifert and Patrick Runge, "Revised refractive index and absorption of In1-xGaxAsyP1-y lattice-matched to InP in transparent and absorption IR-region," Opt. Mater. Express 6, 629-639 (2016)

    IT REQUIRES VALIDATION!!

    """
    Eg = 1.35 - 0.72 * y + 0.12 * y**2

    R = -0.00115 + 0.00191 * Eg
    Gamma = -0.000691 + 0.00433 * Eg
    A = -0.0453 + 2.1103 * Eg
    a = 72.32 + 12.78 * Eg
    b = 4.84 + 4.66 * Eg
    c = -0.015 + 0.02 * Eg
    d = -0.178 + 1.042 * Eg

    def cot(x):
        return 1 / np.tan(x)

    n = np.sqrt(
        1
        + a / (b - (E + 1j * Gamma) ** 2)
        + A
        * np.sqrt(R)
        / (E + 1j * Gamma) ** 2
        * (
            np.log(Eg**2 / (Eg**2 - (E + 1j * Gamma) ** 2))
            + np.pi
            * (
                2 * cot(np.pi * np.sqrt(R / Eg))
                - cot(np.pi * np.sqrt(R / (Eg - (E + 1j * Gamma))))
                - cot(np.pi * np.sqrt(R / (Eg + (E + 1j * Gamma))))
            )
        )
    )

    if type(E) == np.ndarray and type(y) != np.ndarray:
        idx = np.where(E < Eg)
        k = np.sqrt(n**2 + a / (b - E**2)).imag + c * (E - Eg) + d * (E - Eg) ** 2
        k[idx] = np.sqrt(n[idx] ** 2 + a / (b - E[idx] ** 2)).imag

    elif type(E) != np.ndarray and type(y) != np.ndarray:
        if E < Eg:
            c = d = 0

        k = np.sqrt(n**2 + a / (b - E**2)).imag + c * (E - Eg) + d * (E - Eg) ** 2

    elif type(E) != np.ndarray and type(y) == np.ndarray:
        idx = np.where(E < Eg)
        c[idx] = 0
        d[idx] = 0

        k = np.sqrt(n**2 + a / (b - E**2)).imag + c * (E - Eg) + d * (E - Eg) ** 2

    # print(n)
    return n.real + 1j * k


"""
These arguments are MANDATORY for each ElectroOpticalModule

mup: Quantity,
mun: Quantity,
Ec: Quantity,
Ev: Quantity,
Efn: Quantity,
Efp: Quantity,
N: Quantity,
P: Quantity,
Efield: Quantity,
reg: UnitRegistry,

"""
class ElectroOpticalModel:
    """
    A base class for all the electro optical models
    """

class InGaAsPElectroOpticalModel(ElectroOpticalModel):

    n_effects: int = 5  # This is a MUST
    
    def __init__(
        self,
        mup: Quantity,
        mun: Quantity,
        Ec: Quantity,
        Ev: Quantity,
        Efn: Quantity,
        Efp: Quantity,
        N: Quantity,
        P: Quantity,
        Efield: Quantity,
        reg: UnitRegistry,
        T: float = 300,
        y: float = 0,
        bandgap_model: str = "jain",
        BF_model: str = "vinchant",
        growing_direction: str = "z",
    ):
        """
        Base model for In_{1-x}Ga_{x}As_{y}P_{1-y}. It gives the changes in absorption and
        refractive index for charge carrier effects.

        [!!] It expects Quantities for a single voltage
        The kwargs are to accomodate interoperability with other ElectroOpticalModels
        NOTE: the equations used for the modification of the bandstructure are applied assuming validity for compensated material, but that is not the case.
        It  may be a source of error. Beware.

        Ec: Conduction band energy in eV. Assumed to have shape (N)
        Ev: Valence band energy in eV. Assumed to have shape (N)
        Efn: Conduction quasi fermi level energy in eV. Assumed to have shape (N)
        Ep: Valence quasi fermi level energy in eV. Assumed to have shape (N)
        Efield: induced electric field. Assumed to have shape (N, 3).
        n: Electron concentration in cm^-3. Assumed to have shape (N)
        p: Hole concentration in cm^-3. Assumed to have shape (N)
        T: Temperature of operation in kelvin
        y: Concentration. Must be between 0 and 1. Assumed to have shape (Nx).
        growing direction: the axis of growth. Must be x,y or z. Default is z.
        BF_model: if 'vinchant' it will use the data calculated by Vinchant 1992 [7]. If 'BGN' it will use the result from the model while including BGN, and 'no BGN' will not include it.

        References:
        1) http://www.ioffe.ru/SVA/NSM/Semicond/InP/bandstr.html#Masses
        2) Bennett, B.R., R.A. Soref, and J.A. Del Alamo. “Carrier-Induced Change in Refractive Index of InP, GaAs and InGaAsP.” IEEE Journal of Quantum Electronics 26, no. 1 (January 1990): 113–22. https://doi.org/10.1109/3.44924.
        3) Adachi, Sadao. Properties of Group-IV, III-V and II-VI Semiconductors. Wiley Series in Materials for Electronic and Optoelectronic Applications. Chichester, West Sussex, England: John Wiley & Sons, Ltd, 2006.
        4) Sze, S. M., and Kwok Kwok Ng. Physics of Semiconductor Devices. 3rd ed. Hoboken, N.J: Wiley-Interscience, 2007.
        5) Fiedler, F., and A. Schlachetzki. “Optical Parameters of InP-Based Waveguides.” Solid-State Electronics 30, no. 1 (January 1987): 73–83. https://doi.org/10.1016/0038-1101(87)90032-3.
        6) Moss, T. S., Geoffrey John Burrell, and Brian Ellis. Semiconductor Opto-Electronics. London: Butterworths, 1973.
        7) Vinchant, J.-F., J.A. Cavailles, M. Erman, P. Jarry, and M. Renaud. “InP/GaInAsP Guided-Wave Phase Modulators Based on Carrier-Induced Effects: Theory and Experiment.” Journal of Lightwave Technology 10, no. 1 (January 1992): 63–70. https://doi.org/10.1109/50.108738.


        """

        self.n_effects = 5  # This is a MUST

        self.reg = reg

        self.wl = 1550 * self.reg.nanometer
        self.y = y

        self.x = self.y / (2.2020 - 0.0659 * self.y)  # Taken from [5]

        self.bandgap_model = bandgap_model

        self.growing_direction = growing_direction

        self.Ec = Ec
        self.Ev = Ev
        self.Efn = Efn
        self.Efp = Efp
        self.N = N
        self.P = P
        self.Efield = Efield

        self.T = T * self.reg.kelvin

        self.e = 1.602176634e-19 * self.reg.coulomb  # Coulombs
        self.kb = (
            1.380649e-23
            * self.reg.meter**2
            * self.reg.kg
            * self.reg.second**-2
            * self.reg.kelvin**-1
        )  # m^2 kg s^-2 K^-1
        self.e0 = 8.854e-12 * self.reg.farad * self.reg.meter**-1  # F m^-1
        self.h = 6.62607015e-34 * self.reg.joule * self.reg.second  # J Hz^-1
        self.hbar = self.h / (2 * np.pi)
        self.c = 3e8 * self.reg.meter * self.reg.second**-1  # m s^-1
        self.m0 = 9.10e-31 * self.reg.kg  # kg

        self.energy = (2 * np.pi / self.wl * self.hbar * self.c).to(self.reg.eV)

        # Taken from  [5]
        self.me = (0.07 - 0.0308 * self.y) * self.m0
        self.mhh = (0.6 - 0.218 * self.y + 0.07 * self.y**2) * self.m0
        self.mhl = (0.12 - 0.078 * self.y + 0.002 * self.y**2) * self.m0

        # Formulas [2]
        self.Nc = 2 * (
            (self.me * self.kb * self.T / (2 * np.pi * self.hbar**2)) ** 1.5
        ).to(
            self.reg.centimeter**-3
        )  # cm^-3
        self.Nv = 2 * (
            (
                (self.mhh**1.5 + self.mhl**1.5) ** (2 / 3)
                * self.kb
                * self.T
                / (2 * np.pi * self.hbar**2)
            )
            ** 1.5
        ).to(
            self.reg.centimeter**-3
        )  # cm^-3

        self.so = (
            0.119 + 0.30 * self.y - 0.107 * self.y**2
        ) * self.reg.eV  # eV. Taken from [5]

        self.eps_s = self.get_eps_s()

        # Taken from [2]
        self.C = (
            4.4e12
            * self.reg.centimeter**-1
            * self.reg.second**-0.5
            * np.sqrt(self.hbar)
        )  # Taken from Bennet 1990
        # The sqrt(hbar) comes from the fact that C comes from an earlier paper that fits an absorption curve to frequency rather than energy

        # Adapted using formulas from [5] and theoretical predictions from [6]

        mr_InP_hh = (1 / 0.07 + 1 / 0.6) ** -1 * self.m0
        mr_InP_hl = (1 / 0.07 + 1 / 0.12) ** -1 * self.m0
        mr_hh = (1 / self.me + 1 / self.mhh) ** -1
        mr_hl = (1 / self.me + 1 / self.mhl) ** -1

        n0_InP = _get_n(self.energy.magnitude, y=0).real
        n0 = _get_n(self.energy.magnitude, y=y).real

        self.Chh = self.C * (mr_InP_hh**1.5 / (mr_InP_hh**1.5 + mr_InP_hl**1.5))
        self.Chl = self.C * (mr_InP_hl**1.5 / (mr_InP_hh**1.5 + mr_InP_hl**1.5))
        # print('p',self.Chh, self.Chl)

        self.Chh = (mr_hh / mr_InP_hh) ** 1.5 * n0_InP / n0 * self.Chh
        self.Chl = (mr_hl / mr_InP_hh) ** 1.5 * n0_InP / n0 * self.Chl

        self.mue, self.muh = self.get_mobility()

        # Parameters for piezo effects. Taken from [3]
        self.S11 = (
            1.639e-12 * self.reg.centimeter**2 * self.reg.dyne**-1
        )  # mechanical compliance
        self.S12 = -0.589e-12 * self.reg.centimeter**2 * self.reg.dyne**-1
        self.S44 = 2.26e-12 * self.reg.centimeter**2 * self.reg.dyne**-1
        self.e14 = (
            -0.083 * self.reg.coulomb * self.reg.meter**-2
        )  # piezoelectric stress constant.

        # Load data for bandfilling effect
        if BF_model == "BGN":
            data_BF = np.asarray(
                [
                    [0.0, -4.935845233448187e-21, -1.339378005082977e-21],
                    [0.1, -6.052742734846166e-21, -1.5824393501361544e-21],
                    [0.2, -7.588619428674746e-21, -1.8843286756997448e-21],
                    [0.3, -9.73845412951748e-21, -2.265714433315684e-21],
                    [0.4, -1.2962568457500875e-20, -2.7582981460228194e-21],
                    [0.53, -2.0482999504509693e-20, -3.656520005343032e-21],
                    [0.6, -2.826656837654254e-20, -4.340072388303156e-21],
                    [0.7, -5.843171447512438e-20, -5.78945579548802e-21],
                ]
            )
        elif BF_model == "no BGN":
            data_BF = np.asarray(
                [
                    [0.0, -3.519944131135538e-21, -1.2578195480077319e-21],
                    [0.1, -4.283763658235055e-21, -1.4795311461745242e-21],
                    [0.2, -5.276454725508618e-21, -1.7520509288891255e-21],
                    [0.3, -6.597394990306974e-21, -2.0916355997477367e-21],
                    [0.4, -8.40993170474557e-21, -2.522768453275956e-21],
                    [0.53, -1.2013320880787854e-20, -3.2895688460532474e-21],
                    [0.6, -1.4978905917602458e-20, -3.853518305029724e-21],
                    [0.7, -2.167838137390019e-20, -4.982386824381391e-21],
                    [0.8, -3.5522589772741646e-20, -6.935918911360957e-21],
                ]
            )
        elif BF_model == "vinchant":

            data_BF = np.asarray(
                [
                    [-0.003, -5.625e-21, 0.0],
                    [0.047, -6.192e-21, 0.0],
                    [0.089, -6.637e-21, 0.0],
                    [0.127, -7.122000000000001e-21, 0.0],
                    [0.166, -7.446e-21, 0.0],
                    [0.2, -7.932e-21, 0.0],
                    [0.235, -8.255e-21, 0.0],
                    [0.27, -8.984e-21, 0.0],
                    [0.307, -9.469e-21, 0.0],
                    [0.352, -1.085e-20, 0.0],
                    [0.398, -1.206e-20, 0.0],
                    [0.44, -1.3350000000000001e-20, 0.0],
                    [0.473, -1.4650000000000002e-20, 0.0],
                    [0.515, -1.643e-20, 0.0],
                    [0.557, -1.8860000000000002e-20, 0.0],
                    [0.603, -2.1930000000000002e-20, 0.0],
                    [0.636, -2.533e-20, 0.0],
                    [0.669, -2.8250000000000005e-20, 0.0],
                    [0.692, -3.148e-20, 0.0],
                    [0.72, -3.6340000000000004e-20, 0.0],
                    [0.742, -4.0870000000000006e-20, 0.0],
                    [0.762, -4.735e-20, 0.0],
                    [0.778, -5.2850000000000003e-20, 0.0],
                    [0.793, -5.835e-20, 0.0],
                ]
            )

        slope_interp_n = interp1d(data_BF[:, 0], data_BF[:, 1], kind="quadratic")
        slope_interp_p = interp1d(data_BF[:, 0], data_BF[:, 2], kind="quadratic")

        self.slope_P_BF = slope_interp_p(self.y) * self.reg.centimeter**3
        self.slope_N_BF = slope_interp_n(self.y) * self.reg.centimeter**3

    def get_eps_s(self):
        """
        Returns the relative permeability as a function of wavelength and adjusted for the different concentrations.
        It uses the modified single oscillator model and the conversion formulas from the single oscilator model formula as outlined in [1] (Appendix).

        References:
        1) Fiedler, F., and A. Schlachetzki. “Optical Parameters of InP-Based Waveguides.” Solid-State Electronics 30, no. 1 (January 1987): 73–83. https://doi.org/10.1016/0038-1101(87)90032-3.

        """

        bandgap_model = self.bandgap_model
        y = self.y
        wl = self.wl

        x = y / (2.2020 - 0.0659 * y)

        Ed = (28.91 - 9.278 * y + 5.626 * y**2) * self.reg.eV
        E0 = (3.391 - 1.652 * y + 0.863 * y**2 - 0.123 * y**3) * self.reg.eV
        Eg = self.Ec - self.Ev
        E = (2 * np.pi / wl * self.hbar * self.c).to(self.reg.eV)

        a1 = Ed / E0
        a2 = Ed * E**2 / E0**3
        a3 = (
            Ed
            / (2 * E0**3 * (E0**2 - Eg**2))
            * E**4
            * np.log((2 * E0**2 - Eg**2 - E**2) / (Eg**2 - E**2))
        )

        n_sq = 1 + a1 + a2 + a3
        n_sq = n_sq.to(self.reg.dimensionless)

        return n_sq

    def get_BGN(self):
        """
        Return the bandgap narrowing energy of InGaAsP.
        Temperature dependence was taken from [1] page 121.

        BGN is calculated for uncompensated InGaAsP.

        if model='jain' it makes use of the model from [2]. The BGN adaptation for InGaAsP is based on the Vegard's law taken from [5].

        References:
        1) Adachi, Sadao. Properties of Group-IV, III-V and II-VI Semiconductors. Wiley Series in Materials for Electronic and Optoelectronic Applications. Chichester, West Sussex, England: John Wiley & Sons, Ltd, 2006.
        2) Jain, S. C., J. M. McGregor, and D. J. Roulston. “Band‐gap Narrowing in Novel III‐V Semiconductors.” Journal of Applied Physics 68, no. 7 (October 1990): 3747–49. https://doi.org/10.1063/1.346291.
        3) Bennett, B.R., R.A. Soref, and J.A. Del Alamo. “Carrier-Induced Change in Refractive Index of InP, GaAs and InGaAsP.” IEEE Journal of Quantum Electronics 26, no. 1 (January 1990): 113–22. https://doi.org/10.1109/3.44924.
        4) http://www.ioffe.ru/SVA/NSM/Semicond/InP
        5) Fiedler, F., and A. Schlachetzki. “Optical Parameters of InP-Based Waveguides.” Solid-State Electronics 30, no. 1 (January 1987): 73–83. https://doi.org/10.1016/0038-1101(87)90032-3.

        """

        model = self.bandgap_model

        P = self.P
        N = self.N
        x = self.x
        y = self.y
        T = self.T

        if model == "jain":
            constants = {
                "GaAs": {"A": 9.83e-9, "B": 3.90e-7, "C": 3.90e-12},
                "GaP": {"A": 12.7e-9, "B": 5.85e-7, "C": 3.90e-12},
                "InP": {"A": 10.3e-9, "B": 4.43e-7, "C": 3.38e-12},
                "InAs": {"A": 8.34e-9, "B": 2.91e-7, "C": 4.53e-12},
            }

            A = (
                x * y * constants["GaAs"]["A"]
                + x * (1 - y) * constants["GaP"]["A"]
                + y * (1 - x) * constants["InAs"]["A"]
                + (1 - x) * (1 - y) * constants["InP"]["A"]
            )

            B = (
                x * y * constants["GaAs"]["B"]
                + x * (1 - y) * constants["GaP"]["B"]
                + y * (1 - x) * constants["InAs"]["B"]
                + (1 - x) * (1 - y) * constants["InP"]["B"]
            )

            C = (
                x * y * constants["GaAs"]["C"]
                + x * (1 - y) * constants["GaP"]["C"]
                + y * (1 - x) * constants["InAs"]["C"]
                + (1 - x) * (1 - y) * constants["InP"]["C"]
            )

            BGN_p = (
                A * P.magnitude ** (1 / 3)
                + B * P.magnitude**0.25
                + C * P.magnitude**0.5
            ) * self.reg.eV

        elif model == "none":
            BGN_p = 0

        # Calculation for BGN_n

        if model == "jain":
            constants = {
                "GaAs": {"A": 16.5e-9, "B": 2.39e-7, "C": 91.4e-12},
                "GaP": {"A": 10.7e-9, "B": 3.45e-7, "C": 9.97e-12},
                "InP": {"A": 17.2e-9, "B": 2.62e-7, "C": 98.4e-12},
                "InAs": {"A": 14.0e-9, "B": 1.97e-7, "C": 57.9e-12},
            }

            A = (
                x * y * constants["GaAs"]["A"]
                + x * (1 - y) * constants["GaP"]["A"]
                + y * (1 - x) * constants["InAs"]["A"]
                + (1 - x) * (1 - y) * constants["InP"]["A"]
            )

            B = (
                x * y * constants["GaAs"]["B"]
                + x * (1 - y) * constants["GaP"]["B"]
                + y * (1 - x) * constants["InAs"]["B"]
                + (1 - x) * (1 - y) * constants["InP"]["B"]
            )

            C = (
                x * y * constants["GaAs"]["C"]
                + x * (1 - y) * constants["GaP"]["C"]
                + y * (1 - x) * constants["InAs"]["C"]
                + (1 - x) * (1 - y) * constants["InP"]["C"]
            )

            BGN_n = (
                A * N.magnitude ** (1 / 3)
                + B * N.magnitude**0.25
                + C * N.magnitude**0.5
            ) * self.reg.eV

        elif model == "none":
            BGN_n = 0

        return BGN_n + BGN_p

    def get_mobility(self):
        """
        Returns the electron and hole mobility for In_{1-x}Ga_{x}As_{y}P_{1-y} based on an interpolation scheme reported by [1]

        References:
        1- Sotoodeh, M., A. H. Khalid, and A. A. Rezazadeh. “Empirical Low-Field Mobility Model for III–V Compounds Applicable in Device Simulation Codes.” Journal of Applied Physics 87, no. 6 (March 15, 2000): 2890–2900. https://doi.org/10.1063/1.372274.

        """

        x = self.x

        y = self.y

        x_values = np.asarray([0, 0.47, 1])
        values_InGaAs = {
            "mu_max": {"n": [34000, 14000, 9400], "p": [530, 320, 491.5]},
            "mu_min": {"n": [1000, 300, 500], "p": [20, 10, 20]},
            "Nref": {"n": [1.1e18, 1.3e17, 6.0e16], "p": [1.1e17, 4.9e17, 1.48e17]},
            "lambda": {"n": [0.32, 0.48, 0.394], "p": [0.46, 0.403, 0.38]},
            "theta1": {"n": [1.57, 1.59, 2.1], "p": [2.3, 1.59, 2.2]},
            "theta2": {"n": [3.0, 3.68, 3.0], "p": [3.0, 3.0, 3.0]},
        }
        for key1 in values_InGaAs.keys():
            for key2 in values_InGaAs[key1].keys():
                values_InGaAs[key1][key2] = np.asarray(values_InGaAs[key1][key2])

        for key in values_InGaAs.keys():
            if key != "Nref":
                values_InGaAs[key]["n_out"] = make_interp_spline(
                    x_values, values_InGaAs[key]["n"], k=2
                )(x)
                values_InGaAs[key]["p_out"] = make_interp_spline(
                    x_values, values_InGaAs[key]["p"], k=2
                )(x)
            else:
                values_InGaAs[key]["n_out"] = 10 ** make_interp_spline(
                    x_values, np.log10(values_InGaAs[key]["n"]), k=2
                )(x)
                values_InGaAs[key]["p_out"] = 10 ** make_interp_spline(
                    x_values, np.log10(values_InGaAs[key]["p"]), k=2
                )(x)

        x_values = np.asarray([0, 0.51, 1])
        values_InGaP = {
            "mu_max": {"n": [5200, 4300, 152], "p": [170, 150, 147]},
            "mu_min": {"n": [400, 400, 10], "p": [10, 15, 10]},
            "Nref": {"n": [3.0e17, 2.0e16, 4.4e18], "p": [4.87e17, 1.5e17, 1.0e18]},
            "lambda": {"n": [0.47, 0.70, 0.80], "p": [0.62, 0.8, 0.85]},
            "theta1": {"n": [2.0, 1.66, 1.60], "p": [2.0, 2.0, 1.98]},
            "theta2": {"n": [3.25, 0.71], "p": [3.0, 0]},
        }

        for key1 in values_InGaP.keys():
            for key2 in values_InGaP[key1].keys():
                values_InGaP[key1][key2] = np.asarray(values_InGaP[key1][key2])

        for key in values_InGaP.keys():
            if key not in ["Nref", "theta2"]:
                values_InGaP[key]["n_out"] = make_interp_spline(
                    x_values, values_InGaP[key]["n"], k=2
                )(x)
                values_InGaP[key]["p_out"] = make_interp_spline(
                    x_values, values_InGaP[key]["p"], k=2
                )(x)
            elif key == "Nref":
                values_InGaP[key]["n_out"] = 10 ** make_interp_spline(
                    x_values, np.log10(values_InGaP[key]["n"]), k=2
                )(x)
                values_InGaP[key]["p_out"] = 10 ** make_interp_spline(
                    x_values, np.log10(values_InGaP[key]["p"]), k=2
                )(x)
            else:
                values_InGaP[key]["n_out"] = make_interp_spline(
                    [0, 1], values_InGaP[key]["n"], k=1
                )(x)
                values_InGaP[key]["p_out"] = make_interp_spline(
                    [0, 1], values_InGaP[key]["p"], k=1
                )(x)

        values = {
            "mu_max": {
                "n": (
                    y * values_InGaAs["mu_max"]["n_out"]
                    + (1 - y) * values_InGaP["mu_max"]["n_out"]
                )
                / (1 + 6 * y * (1 - y)),
                "p": (
                    y * values_InGaAs["mu_max"]["p_out"]
                    + (1 - y) * values_InGaP["mu_max"]["p_out"]
                )
                / (1 + 6 * y * (1 - y)),
            },
            "mu_min": {
                "n": (
                    y * values_InGaAs["mu_min"]["n_out"]
                    + (1 - y) * values_InGaP["mu_min"]["n_out"]
                )
                / (1 + 6 * y * (1 - y)),
                "p": (
                    y * values_InGaAs["mu_min"]["p_out"]
                    + (1 - y) * values_InGaP["mu_min"]["p_out"]
                ),
            },
            "Nref": {
                "n": 10
                ** (
                    y * np.log10(values_InGaAs["Nref"]["n_out"])
                    + (1 - y) * np.log10(values_InGaP["Nref"]["n_out"])
                ),
                "p": 10
                ** (
                    y * np.log10(values_InGaAs["Nref"]["p_out"])
                    + (1 - y) * np.log10(values_InGaP["Nref"]["p_out"])
                ),
            },
            "lambda": {
                "n": (
                    y * values_InGaAs["lambda"]["n_out"]
                    + (1 - y) * values_InGaP["lambda"]["n_out"]
                ),
                "p": (
                    y * values_InGaAs["lambda"]["p_out"]
                    + (1 - y) * values_InGaP["lambda"]["p_out"]
                ),
            },
            "theta1": {
                "n": (
                    y * values_InGaAs["theta1"]["n_out"]
                    + (1 - y) * values_InGaP["theta1"]["n_out"]
                )
                / (1 + 1 * y * (1 - y)),
                "p": (
                    y * values_InGaAs["theta1"]["p_out"]
                    + (1 - y) * values_InGaP["theta1"]["p_out"]
                )
                / (1 + 1 * y * (1 - y)),
            },
            "theta2": {
                "n": (
                    y * values_InGaAs["theta2"]["n_out"]
                    + (1 - y) * values_InGaP["theta2"]["n_out"]
                ),
                "p": (
                    y * values_InGaAs["theta2"]["p_out"]
                    + (1 - y) * values_InGaP["theta2"]["p_out"]
                ),
            },
        }

        T = self.T.to(self.reg.kelvin).magnitude
        N = self.N.to(self.reg.centimeter**-3).magnitude
        P = self.P.to(self.reg.centimeter**-3).magnitude
        
        mobility_n = values["mu_min"]["n"] + (
            values["mu_max"]["n"] * (300 / T) ** values["theta1"]["n"]
            - values["mu_min"]["n"]
        ) / (
            1
            + (N / (values["Nref"]["n"] * (300 / T) ** values["theta2"]["n"]))
            ** values["lambda"]["n"]
        )

        mobility_p = values["mu_min"]["p"] + (
            values["mu_max"]["p"] * (300 / T) ** values["theta1"]["p"]
            - values["mu_min"]["p"]
        ) / (
            1
            + (P / (values["Nref"]["p"] * (300 / T) ** values["theta2"]["p"]))
            ** values["lambda"]["p"]
        )

        return (
            mobility_n * self.reg.centimeter**2 / (self.reg.volt * self.reg.second),
            mobility_p * self.reg.centimeter**2 / (self.reg.volt * self.reg.second),
        )

    def FD(self, E, Ef=1):
        return 1 / (1 + np.exp((E - Ef) / (self.kb * self.T)))

    def get_dalpha_BF(self):
        """
        Gives the change in absorption due to the bandfilling effect.
        The treatment if based on [1] but we neglect the quasi-fermi levels and use the fermi level numerically calculated instead.

        Note that if the bandgap_model is not 'none' the BGN effect is taken into consideration within the bandfilling.

        E must be a Pint quantity in energy using the parent object's register!!

        References:
        1) Bennett, B.R., R.A. Soref, and J.A. Del Alamo. “Carrier-Induced Change in Refractive Index of InP, GaAs and InGaAsP.” IEEE Journal of Quantum Electronics 26, no. 1 (January 1990): 113–22. https://doi.org/10.1109/3.44924.

        """

        E = np.asarray([self.energy.magnitude])[..., None, None] * self.energy.units

        Eg = (
            self.Ec - self.Ev
        )  - self.get_BGN()

        Eg = Eg[None, ...]

        Eah = (Eg - E) * (self.me / (self.me + self.mhh)) - Eg + self.Ec
        Eal = (Eg - E) * (self.me / (self.me + self.mhl)) - Eg + self.Ec
        Ebh = (E - Eg) * (self.mhh / (self.me + self.mhh)) + self.Ec
        Ebl = (E - Eg) * (self.mhl / (self.me + self.mhl)) + self.Ec

        # Eah=Eah[None, ...]
        # Eal=Eal[None, ...]
        # Ebh=Ebh[None, ...]
        # Ebl=Ebl[None, ...]
        # print(E.shape, Eg.shape)
        alpha0 = self.Chh / E * np.sqrt(E - Eg) + self.Chl / E * np.sqrt(E - Eg)
        alpha = self.Chh / E * np.sqrt(E - Eg) * (
            self.FD(Eah, Ef=self.Efp[None, ...]) - self.FD(Ebh, Ef=self.Efn[None, ...])
        ) + self.Chl / E * np.sqrt(E - Eg) * (
            self.FD(Eal, Ef=self.Efp[None, ...]) - self.FD(Ebl, Ef=self.Efn[None, ...])
        )
        # print(alpha.shape, alpha0.shape)
        # print(E.shape, Eg.shape)
        alpha = np.nan_to_num(alpha)
        alpha0 = np.nan_to_num(alpha0)

        return np.squeeze((alpha - alpha0)).to(self.reg.centimeter**-1)

    def get_alpha_sqrt(self, E=None, bandgap_model=None):
        E = self.energy
        bandgap_model = self.bandgap_model

        Eg = self.Ec - self.Ev - self.get_BGN()

        alpha0 = self.Chh / E * np.sqrt(E - Eg) + self.Chl / E * np.sqrt(E - Eg)

        return alpha0.to(self.reg.centimeter**-1)

    def get_dn_BF(self):
        """
        Returns the change in refractive index based only on the band filling effect based on the kramers kronig relations.


        """

        return self.slope_P_BF * self.P + self.slope_N_BF * self.N

    def get_dalpha_plasma(self):
        """
        Gives the change in absorption due to the plasma effect as reported in [1]. Simply put, the model is based on a second order perturbation theory and considers 3 scattering mechanisms: electron - optical phonon, electron - acoustical phonon and electron - ionized impurity. It assumes room temperature of 300K. It is limited to concentrations below 6e18 cm^-3 and above 1116cm^-3

        E must be a Pint quantity in energy using the parent object's register!!

        The values predicted by the below formula are consistent with the experimental results from [2].

        References:
        1) Walukiewicz, W., J. Lagowski, L. Jastrzebski, P. Rava, M. Lichtensteiger, C. H. Gatos, and H. C. Gatos. “Electron Mobility and Free-Carrier Absorption in InP; Determination of the Compensation Ratio.” Journal of Applied Physics 51, no. 5 (May 1, 1980): 2659–68. https://doi.org/10.1063/1.327925.

        2) Dumke, W. P., M. R. Lorenz, and G. D. Pettit. “Intra- and Interband Free-Carrier Absorption and the Fundamental Absorption Edge in n -Type InP.” Physical Review B 1, no. 12 (June 15, 1970): 4668–73. https://doi.org/10.1103/PhysRevB.1.4668.

        """

        E = self.energy

        N1 = np.copy(self.N.magnitude) * self.N.units

        # Account for limitation on doping
        idx1 = np.where(N1.to(self.reg.centimeter**-3).magnitude > 6e18)
        idx2 = np.where(N1.to(self.reg.centimeter**-3).magnitude < 1e16)

        N1[idx1] = 6e18 * self.reg.centimeter**-3
        N1[idx2] = 1e16 * self.reg.centimeter**-3

        wv = self.hbar * 2 * np.pi * self.c / E

        dopings = np.asarray(
            [1e16, 1.5e16, 2e16, 3e16, 4e16, 5e16, 6e16, 7e16, 8e16, 9e16]
            + [1e17, 1.5e17, 2e17, 3e17, 4e17, 5e17, 6e17, 7e17, 8e17, 9e17]
            + [1e18, 1.5e18, 2e18, 3e18, 4e18, 5e18, 6e18]
        )

        alpha_imp = np.asarray(
            [
                0.004,
                0.008,
                0.014,
                0.031,
                0.056,
                0.086,
                0.123,
                0.167,
                0.217,
                0.273,
                0.314,
                0.690,
                1.201,
                2.602,
                4.474,
                6.790,
                9.510,
                12.64,
                16.13,
                20.00,
                24.22,
                50.28,
                93.91,
                170.3,
                276.7,
                396.1,
                522.8,
            ]
        )

        alpha_ac = np.asarray(
            [
                0.034,
                0.052,
                0.069,
                0.104,
                0.139,
                0.173,
                0.208,
                0.243,
                0.278,
                0.313,
                0.325,
                0.491,
                0.660,
                1.005,
                1.360,
                1.726,
                2.100,
                2.488,
                2.879,
                3.285,
                3.699,
                5.912,
                8.354,
                13.87,
                20.12,
                26.98,
                34.34,
            ]
        )

        alpha_op = np.asarray(
            [
                0.623,
                0.932,
                1.239,
                1.850,
                2.456,
                3.051,
                3.646,
                4.240,
                4.815,
                5.397,
                5.578,
                8.227,
                10.79,
                15.75,
                20.52,
                25.16,
                29.65,
                34.11,
                38.44,
                42.75,
                47.01,
                67.80,
                88.02,
                127.0,
                164.1,
                199.6,
                233.4,
            ]
        )

        alpha_imp_interp = lambda x: interp1d(
            np.log10(dopings), alpha_imp, kind="linear"
        )(np.log10(x))
        alpha_ac_interp = lambda x: interp1d(
            np.log10(dopings), alpha_ac, kind="linear"
        )(np.log10(x))
        alpha_op_interp = lambda x: interp1d(
            np.log10(dopings), alpha_op, kind="linear"
        )(np.log10(x))

        lam0 = 10e-6 * self.reg.meter

        wv_ratio = (wv / lam0).to(self.reg.dimensionless).magnitude

        doping = N1.to(self.reg.centimeter**-3).magnitude

        alpha = (
            alpha_imp_interp(doping) * (wv_ratio) ** 3.5
            + alpha_op_interp(doping) * (wv_ratio) ** 2.5
            + alpha_ac_interp(doping) * (wv_ratio) ** 1.5
        )

        return alpha * self.reg.centimeter**-1

    def get_dn_plasma(self, E=None):
        """
        Returns the change in refractive index based on the plasma effect. It makes use of [1], page 79, eq.5.2.14

        E: float. Pint quantity in energy using the parent object's register.

        References:
        1) Hunsperger, Robert G. Integrated Optics: Theory and Technology. 6th ed. Advanced Texts in Physics. New York London: Springer, 2009.

        """
        if E is None:
            E = self.energy

        n_contribution = (
            -1
            / 2
            * (
                self.N
                * self.e**2
                / (self.me * self.e0 * E**2 / self.hbar**2 * np.sqrt(self.eps_s))
            ).to(self.reg.dimensionless)
        )
        # mass=(self.mhh**-2+self.mhl**-2)**-0.5
        # p_contribution=(-1/2 * (self.P*self.e**2/(mass*self.e0*E**2/self.hbar**2*np.sqrt(self.eps_s)))).to(self.reg.dimensionless)

        return n_contribution

    def get_dalpha_iv(self, E=None):
        """
        Returns the intervalence absorption component. This is calculated from [1], eq. 16.

        This formula may give worse results at low dopings and underpredict absorption. [2]

        E: float. Pint quantity in energy using the parent object's register.

        References:
        1) Weber, J.-P. “Optimization of the Carrier-Induced Effective Index Change in InGaAsP Waveguides-Application to Tunable Bragg Filters.” IEEE Journal of Quantum Electronics 30, no. 8 (August 1994): 1801–16. https://doi.org/10.1109/3.301645.
        2) Casey, H. C., and P. L. Carter. “Variation of Intervalence Band Absorption with Hole Concentration in p -Type InP.” Applied Physics Letters 44, no. 1 (January 1, 1984): 82–83. https://doi.org/10.1063/1.94561.


        """

        if E is None:
            E = self.energy

        return (
            (
                4.252e-20
                * np.exp(E.to(self.reg.eV).magnitude * -3.657)
                * self.P.to(self.reg.meter**-3).magnitude
            )
            * self.reg.meter**-1
        ).to(self.reg.centimeter**-1)

    def get_dn_iv(self, E=None):
        """
        Returns the intervalence absorption component. This is calculated from [1], eq. 17 and 18.

        E must be a Pint quantity in energy using the parent object's register!!

        References:
        1) Weber, J.-P. “Optimization of the Carrier-Induced Effective Index Change in InGaAsP Waveguides-Application to Tunable Bragg Filters.” IEEE Journal of Quantum Electronics 30, no. 8 (August 1994): 1801–16. https://doi.org/10.1109/3.301645.

        """

        if E is None:
            E = self.energy

        alpha0 = 4.252e-20 * self.reg.meter**2
        b = 3.657 * self.reg.eV**-1

        # note that there is no e in the denominator because that is accounted for in the Pint quantity.
        return (
            -self.hbar
            * self.c
            / np.pi
            * alpha0
            * 1
            / (2 * E)
            * (
                np.exp(-(b * E).to(self.reg.dimensionless).magnitude)
                * expi((b * E).to(self.reg.dimensionless).magnitude)
                + np.exp((b * E).to(self.reg.dimensionless).magnitude)
                * exp1((b * E).to(self.reg.dimensionless).magnitude)
            )
            * self.P
        ).to(self.reg.dimensionless)

    def get_dperm_pockels(self):
        """
        Returns the change in refractive index owing to the pockels effect.

        E: energy of the field's photons. Pint quantity
        Efield: Electric field components. It is assumed to have dimensions (Nx, 3)

        References:
        1) Adachi, Sadao, and Kunishige Oe. “Internal Strain and Photoelastic Effects in Ga  1− x  Al  x  As/GaAs and In  1− x  Ga  x  As  y  P  1− y  /InP Crystals.” Journal of Applied Physics 54, no. 11 (November 1983): 6620–27. https://doi.org/10.1063/1.331898.
        2) Adachi, Sadao, and Kunishige Oe. “Linear Electro‐optic Effects in Zincblende‐type Semiconductors: Key Properties of InGaAsP Relevant to Device Design.” Journal of Applied Physics 56, no. 1 (July 1984): 74–80. https://doi.org/10.1063/1.333731.

        """
        E = self.energy
        Efield = (
            self.Efield.T
        )  # The transpose stems from the fact that I'm reusing old code that requires an Efield of shape (3, N)
        bandgap_model = self.bandgap_model

        def g(chi):
            return 1 / chi**2 * (2 - (1 + chi) ** -0.5 - (1 - chi) ** -0.5)

        def f(chi):
            return 1 / chi**2 * (2 - (1 + chi) ** 0.5 - (1 - chi) ** 0.5)

        def symmetric(x):
            """
            Returns an array in the standard symmetric form instead of the voigt notation. Assumes shape of [6, Nx]

            """
            y = np.zeros((3, 3) + x.shape[-1:]) * x.units

            y[0, 0] = x[0]
            y[0, 1] = y[1, 0] = x[5]
            y[0, 2] = y[2, 0] = x[4]
            y[1, 1] = x[1]
            y[1, 2] = y[2, 1] = x[3]
            y[2, 2] = x[2]

            return y

        constants = {
            "InP": {
                "E0": -42.06e-12 * self.reg.meter * self.reg.volt**-1,
                "F0": 91.32e-12 * self.reg.meter * self.reg.volt**-1,
                "C": -0.36e-10 * self.reg.meter**2 * self.reg.newton**-1,
                "D": 2.60e-10 * self.reg.meter**2 * self.reg.newton**-1,
            },
            "GaP": {
                "E0": -83.31e-12 * self.reg.meter * self.reg.volt**-1,
                "F0": 16.60e-12 * self.reg.meter * self.reg.volt**-1,
                "C": -0.06e-10 * self.reg.meter**2 * self.reg.newton**-1,
                "D": 1.92e-10 * self.reg.meter**2 * self.reg.newton**-1,
            },
            "GaAs": {
                "E0": -71.48e-12 * self.reg.meter * self.reg.volt**-1,
                "F0": 123.16e-12 * self.reg.meter * self.reg.volt**-1,
                "C": -0.21e-10 * self.reg.meter**2 * self.reg.newton**-1,
                "D": 2.12e-10 * self.reg.meter**2 * self.reg.newton**-1,
            },
            "InAs": {
                "E0": -30.23e-12 * self.reg.meter * self.reg.volt**-1,
                "F0": 197.88e-12 * self.reg.meter * self.reg.volt**-1,
                "C": -1.48e-10 * self.reg.meter**2 * self.reg.newton**-1,
                "D": 2.32e-10 * self.reg.meter**2 * self.reg.newton**-1,
            },
        }

        E0 = (
            constants["InP"]["E0"] * (1 - self.x) * (1 - self.y)
            + constants["GaP"]["E0"] * self.x * (1 - self.y)
            + constants["GaAs"]["E0"] * self.x * self.y
            + constants["InAs"]["E0"] * (1 - self.x) * self.y
        )

        F0 = (
            constants["InP"]["F0"] * (1 - self.x) * (1 - self.y)
            + constants["GaP"]["F0"] * self.x * (1 - self.y)
            + constants["GaAs"]["F0"] * self.x * self.y
            + constants["InAs"]["F0"] * (1 - self.x) * self.y
        )

        C = (
            constants["InP"]["C"] * (1 - self.x) * (1 - self.y)
            + constants["GaP"]["C"] * self.x * (1 - self.y)
            + constants["GaAs"]["C"] * self.x * self.y
            + constants["InAs"]["C"] * (1 - self.x) * self.y
        )

        D = (
            constants["InP"]["D"] * (1 - self.x) * (1 - self.y)
            + constants["GaP"]["D"] * self.x * (1 - self.y)
            + constants["GaAs"]["D"] * self.x * self.y
            + constants["InAs"]["D"] * (1 - self.x) * self.y
        )

        Eg = self.Ec - self.Ev - self.get_BGN()

        r41_free = -1 / self.eps_s**2 * (E0 * g(E / Eg) + F0)
        r41_piezo = (
            -1
            / self.eps_s**2
            * (
                C
                * (
                    -g(E / Eg)
                    + 4
                    * Eg
                    / self.so
                    * (f(E / Eg) - (Eg / (Eg + self.so)) ** 1.5 * f(E / (Eg + self.so)))
                )
                + D
            )
            * self.e14
        )

        pockels_tensor = (
            np.zeros((6, 3, *Efield.shape[1:])) * r41_free / r41_free.magnitude
        )

        pockels_tensor[3, 0] = r41_free + r41_piezo
        pockels_tensor[4, 1] = r41_free + r41_piezo
        pockels_tensor[5, 2] = r41_free + r41_piezo

        # Evaluate the change in impermeability
        # print(Efield.shape, pockels_tensor.shape)
        deta = np.einsum("ilj,lj->ij", pockels_tensor, Efield)
        deta = symmetric(deta)  # return to symmetric matrix

        # build the permitivity tensor
        perm = np.zeros((3, 3, *Efield.shape[1:])) * self.e0.units
        perm[0, 0] = self.eps_s * self.e0
        perm[1, 1] = self.eps_s * self.e0
        perm[2, 2] = self.eps_s * self.e0

        # find the change in permitivity

        dperm = np.einsum("ikt,klt,ljt->ijt", perm, deta, perm) * -1 / self.e0

        return dperm.to(self.e0.units)

    def get_dperm_kerr(self):
        """
        Returns the change in refractive index owing to the Kerr effect.

        E: energy of the field's photons. Pint quantity
        Efield: Electric field components. It is assumed to have dimensions (Nx, 3)

        References:
        [1] - Maat, Derk Hendrik Pieter. InP-Based Integrated MZI Switches for Optical Communication, 2001.


        """

        E = self.energy
        Efield = self.Efield.T

        def symmetric(x):
            """
            Returns an array in the standard symmetric form instead of the voigt notation. Assumes shape of [6, Nx]

            """
            y = np.zeros((3, 3) + x.shape[-1:]) * x.units

            y[0, 0] = x[0]
            y[0, 1] = y[1, 0] = x[5]
            y[0, 2] = y[2, 0] = x[4]
            y[1, 1] = x[1]
            y[1, 2] = y[2, 1] = x[3]
            y[2, 2] = x[2]

            return y

        # Taken from [1]
        A_TE = 0.25e3
        A_TM = 0.20e3
        B_TE = 0.71e9
        B_TM = 0.48e9

        # C_TE=-1.79e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2
        # C_TM=-1.82e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2

        C_TE = -3.10e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2
        C_TM = -5.60e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2

        A_mat = (
            np.asarray(
                [[A_TE, 0   , 0], 
                 [0   , A_TM, 0], 
                 [0   , 0   , A_TE]]
            )
            * self.reg.eV
            / (self.reg.volt * self.reg.meter)
        )

        B_mat = (
            np.asarray(
                [[B_TE, 0   , 0], 
                 [0   , B_TM, 0], 
                 [0   , 0   , B_TE]]
            )
            * self.reg.eV ** (-3 / 2)
            * self.reg.volt
            / self.reg.meter
        )

        freq = E / self.hbar / (2 * np.pi)
        wl = self.c / freq

        Eg = self.Ec - self.Ev - self.get_BGN()

        # print((2*np.pi*self.c/(Eg/self.hbar)).to(self.reg.nanometer))
        Efield_mag_sq = np.einsum("ij,ij -> j", Efield.conjugate(), Efield).real
        
        dalpha = (
            A_mat[..., None]
            * wl
            * np.sqrt(Efield_mag_sq)[None, None, ...]
            / (Eg[None, None, ...] - E)
            * 10
            ** (
                -(
                    B_mat[..., None]
                    * (Eg[None, None, ...] - E) ** (3 / 2)
                    / np.sqrt(Efield_mag_sq)[None, None, ...]
                )
            )
        )
        dalpha = dalpha.to(self.reg.meter**-1)

        dperm_imag = (
            np.sqrt(self.eps_s) * self.c / (2 * np.pi * freq) * dalpha * self.e0
        )

        S11 = C_TE * E**2 / (self.eps_s**2 * (Eg**2 - E**2) ** 2)
        S12 = C_TM * E**2 / (self.eps_s**2 * (Eg**2 - E**2) ** 2)
        #         print(S11.min(), S12.min())

        #         print(E)
        #         print(Eg.max(), Eg.min())

        S11 = S11.to(self.reg.meter**2 / self.reg.volt**2).magnitude
        S12 = S12.to(self.reg.meter**2 / self.reg.volt**2).magnitude
        S00 = np.zeros(S11.shape)

        S_mat = (
            np.asarray(
                [
                    [S11, S12, S12, S00, S00, S00],
                    [S12, S11, S12, S00, S00, S00],
                    [S12, S12, S11, S00, S00, S00],
                    [S00, S00, S00, S11 - S12, S00, S00],
                    [S00, S00, S00, S00, S11 - S12, S00],
                    [S00, S00, S00, S00, S00, S11 - S12],
                ]
            )
            * self.reg.meter**2
            / self.reg.volt**2
        )

        Efield_voigt = np.zeros((6, *Efield.shape[1:])) * Efield.units**2
        Efield_voigt[0] = Efield[0] * Efield[0]
        Efield_voigt[1] = Efield[1] * Efield[1]
        Efield_voigt[2] = Efield[2] * Efield[2]
        Efield_voigt[3] = Efield[1] * Efield[2]
        Efield_voigt[4] = Efield[0] * Efield[2]
        Efield_voigt[5] = Efield[0] * Efield[1]

        # Return to symmetric shape
        deta_real = symmetric(np.einsum("ijk,jk->ik", S_mat, Efield_voigt))

        # build the permitivity tensor
        perm = np.zeros((3, 3, *Efield.shape[1:])) * self.e0.units
        perm[0, 0] = self.eps_s * self.e0
        perm[1, 1] = self.eps_s * self.e0
        perm[2, 2] = self.eps_s * self.e0

        # find the change in permitivity
        dperm = (
            np.einsum("ikt,klt,ljt->ijt", perm, deta_real, perm) * -1 / self.e0
            + 1j * dperm_imag
        )
        
        np.nan_to_num(
            dperm, 0
        )  # This is to avoid some nan values that happen where the field is 0. This causes some numerical errors

        return dperm.to(self.e0.units)

    def get_dperm(
            self,
            fractions: bool = False,
        ) -> np.ndarray:
        """This function returns the change in permitivity tensor"""

        eps_s = self.get_eps_s()

        alpha = self.get_alpha_sqrt()
        alpha = np.nan_to_num(alpha, nan=0)

        n0 = np.sqrt(eps_s)
        k0 = (alpha / 2 / (2 * np.pi / self.wl)).to(self.reg.dimensionless)

        dalpha_BF = self.get_dalpha_BF()
        dalpha_plasma = self.get_dalpha_plasma()
        dalpha_iv = self.get_dalpha_iv()

        dk_BF = (dalpha_BF / 2 / (2 * np.pi / self.wl)).to(self.reg.dimensionless)
        dk_plasma = (dalpha_plasma / 2 / (2 * np.pi / self.wl)).to(
            self.reg.dimensionless
        )
        dk_iv = (dalpha_iv / 2 / (2 * np.pi / self.wl)).to(self.reg.dimensionless)

        dn_BF = self.get_dn_BF()
        dn_plasma = self.get_dn_plasma()
        dn_iv = self.get_dn_iv()

        perm0 = (
            np.asarray(
                [
                    [
                        eps_s.to(self.reg.dimensionless).magnitude,
                        np.zeros(eps_s.shape),
                        np.zeros(eps_s.shape),
                    ],
                    [
                        np.zeros(eps_s.shape),
                        eps_s.to(self.reg.dimensionless).magnitude,
                        np.zeros(eps_s.shape),
                    ],
                    [
                        np.zeros(eps_s.shape),
                        np.zeros(eps_s.shape),
                        eps_s.to(self.reg.dimensionless).magnitude,
                    ],
                ]
            )
            * self.reg.dimensionless
        )

        dperm_BF_1d = 2*n0*dn_BF + 1j*dk_BF

        dperm_plasma_1d = 2*n0*dn_plasma + 1j*dk_plasma

        dperm_iv_1d = 2*n0*dn_iv + 1j*dk_iv

        dperm_BF = (
            np.asarray(
                [
                        [
                            dperm_BF_1d.to(self.reg.dimensionless).magnitude,
                            np.zeros(dperm_BF_1d.shape),
                            np.zeros(dperm_BF_1d.shape),
                        ],
                        [
                            np.zeros(dperm_BF_1d.shape),
                            dperm_BF_1d.to(self.reg.dimensionless).magnitude,
                            np.zeros(dperm_BF_1d.shape),
                        ],
                        [
                            np.zeros(dperm_BF_1d.shape),
                            np.zeros(dperm_BF_1d.shape),
                            dperm_BF_1d.to(self.reg.dimensionless).magnitude,
                        ],
                ]
            )
            * self.reg.dimensionless
        )

        dperm_plasma = (
            np.asarray(
                    [
                        [
                            dperm_plasma_1d.to(self.reg.dimensionless).magnitude,
                            np.zeros(dperm_plasma_1d.shape),
                            np.zeros(dperm_plasma_1d.shape),
                        ],
                        [
                            np.zeros(dperm_plasma_1d.shape),
                            dperm_plasma_1d.to(self.reg.dimensionless).magnitude,
                            np.zeros(dperm_plasma_1d.shape),
                        ],
                        [
                            np.zeros(dperm_plasma_1d.shape),
                            np.zeros(dperm_plasma_1d.shape),
                            dperm_plasma_1d.to(self.reg.dimensionless).magnitude,
                        ],
                    ]
            )
            * self.reg.dimensionless
        )

        dperm_iv = (
            np.asarray(
                
                    [
                        [
                            dperm_iv_1d.to(self.reg.dimensionless).magnitude,
                            np.zeros(dperm_iv_1d.shape),
                            np.zeros(dperm_iv_1d.shape),
                        ],
                        [
                            np.zeros(dperm_iv_1d.shape),
                            dperm_iv_1d.to(self.reg.dimensionless).magnitude,
                            np.zeros(dperm_iv_1d.shape),
                        ],
                        [
                            np.zeros(dperm_iv_1d.shape),
                            np.zeros(dperm_iv_1d.shape),
                            dperm_iv_1d.to(self.reg.dimensionless).magnitude,
                        ],
                    ]
                
            )
            * self.reg.dimensionless
        )

        dperm_pockels = self.get_dperm_pockels() / self.e0
        dperm_kerr = self.get_dperm_kerr() / self.e0

        dperm_BF = np.nan_to_num(
            dperm_BF, nan = 0
        ) 

        dperm_plasma = np.nan_to_num(
            dperm_plasma, nan = 0
        ) 

        dperm_iv = np.nan_to_num(
            dperm_iv, nan = 0
        ) 

        dperm_pockels = np.nan_to_num(
            dperm_pockels, nan = 0
        ) 

        dperm_kerr = np.nan_to_num(
            dperm_kerr, nan = 0
        ) 


        # print(np.isnan(dperm_pockels.imag))
        if not fractions:
            return dperm_BF + dperm_plasma + dperm_iv + dperm_pockels + dperm_kerr
        else:
            return (
                ['Bandfilling', 'Plasma', 'Intervalence', 'Pockels', 'Kerr'],
                dperm_BF,
                dperm_plasma,
                dperm_iv,
                dperm_pockels,
                dperm_kerr,
            )


    #def get_dperm_fractions(self) -> np.ndarray:


        # eps_s = self.get_eps_s()

        # alpha = self.get_alpha_sqrt()
        # alpha = np.nan_to_num(alpha, nan=0)

        # n0 = np.sqrt(eps_s)
        # k0 = (alpha / 2 / (2 * np.pi / self.wl)).to(self.reg.dimensionless)

        # dalpha_BF = self.get_dalpha_BF()
        # dalpha_plasma = self.get_dalpha_plasma()
        # dalpha_iv = self.get_dalpha_iv()

        # dk_BF = (dalpha_BF / 2 / (2 * np.pi / self.wl)).to(self.reg.dimensionless)
        # dk_plasma = (dalpha_plasma / 2 / (2 * np.pi / self.wl)).to(
        #     self.reg.dimensionless
        # )
        # dk_iv = (dalpha_iv / 2 / (2 * np.pi / self.wl)).to(self.reg.dimensionless)

        # dn_BF = self.get_dn_BF()
        # dn_plasma = self.get_dn_plasma()
        # dn_iv = self.get_dn_iv()

        # dperm_BF_1d = 2*n0*dn_BF + 1j*dk_BF

        # dperm_plasma_1d = 2*n0*dn_plasma + 1j*dk_plasma

        # dperm_iv_1d = 2*n0*dn_iv + 1j*dk_iv

        # dperm_BF = (
        #     np.asarray(
        #         [
        #                 [
        #                     dperm_BF_1d.to(self.reg.dimensionless).magnitude,
        #                     np.zeros(dperm_BF_1d.shape),
        #                     np.zeros(dperm_BF_1d.shape),
        #                 ],
        #                 [
        #                     np.zeros(dperm_BF_1d.shape),
        #                     dperm_BF_1d.to(self.reg.dimensionless).magnitude,
        #                     np.zeros(dperm_BF_1d.shape),
        #                 ],
        #                 [
        #                     np.zeros(dperm_BF_1d.shape),
        #                     np.zeros(dperm_BF_1d.shape),
        #                     dperm_BF_1d.to(self.reg.dimensionless).magnitude,
        #                 ],
        #         ]
        #     )
        #     * self.reg.dimensionless
        # )

        # dperm_plasma = (
        #     np.asarray(
        #             [
        #                 [
        #                     dperm_plasma_1d.to(self.reg.dimensionless).magnitude,
        #                     np.zeros(dperm_plasma_1d.shape),
        #                     np.zeros(dperm_plasma_1d.shape),
        #                 ],
        #                 [
        #                     np.zeros(dperm_plasma_1d.shape),
        #                     dperm_plasma_1d.to(self.reg.dimensionless).magnitude,
        #                     np.zeros(dperm_plasma_1d.shape),
        #                 ],
        #                 [
        #                     np.zeros(dperm_plasma_1d.shape),
        #                     np.zeros(dperm_plasma_1d.shape),
        #                     dperm_plasma_1d.to(self.reg.dimensionless).magnitude,
        #                 ],
        #             ]
        #     )
        #     * self.reg.dimensionless
        # )

        # dperm_iv = (
        #     np.asarray(
                
        #             [
        #                 [
        #                     dperm_iv_1d.to(self.reg.dimensionless).magnitude,
        #                     np.zeros(dperm_iv_1d.shape),
        #                     np.zeros(dperm_iv_1d.shape),
        #                 ],
        #                 [
        #                     np.zeros(dperm_iv_1d.shape),
        #                     dperm_iv_1d.to(self.reg.dimensionless).magnitude,
        #                     np.zeros(dperm_iv_1d.shape),
        #                 ],
        #                 [
        #                     np.zeros(dperm_iv_1d.shape),
        #                     np.zeros(dperm_iv_1d.shape),
        #                     dperm_iv_1d.to(self.reg.dimensionless).magnitude,
        #                 ],
        #             ]
                
        #     )
        #     * self.reg.dimensionless
        # )

        # dperm_pockels = self.get_dperm_pockels() / self.e0
        # dperm_kerr = self.get_dperm_kerr() / self.e0

        # dperm_BF = np.nan_to_num(
        #     dperm_BF, nan = 0
        # ) 

        # dperm_plasma = np.nan_to_num(
        #     dperm_plasma, nan = 0
        # ) 

        # dperm_iv = np.nan_to_num(
        #     dperm_iv, nan = 0
        # ) 

        # dperm_pockels = np.nan_to_num(
        #     dperm_pockels, nan = 0
        # ) 

        # dperm_kerr = np.nan_to_num(
        #     dperm_kerr, nan = 0
        # ) 

        # # print(np.isnan(dperm_pockels.imag))
        # return (dperm_BF, dperm_plasma, dperm_iv, dperm_pockels, dperm_kerr)



