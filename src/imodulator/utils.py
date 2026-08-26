"""
This is simply a place to store some helper functions which we have found useful while working with the Imodulator.

Have fun
"""

import numpy as np
import pint
from pint import UnitRegistry
from scipy.integrate import quad, quad_vec, simpson
from scipy.interpolate import interp1d, make_interp_spline
from scipy.optimize import fsolve, root_scalar
from scipy.special import exp1, expi, roots_legendre
import os


def get_n(E, y=0):

    """
    code to generate the real refractive index as well as the 
absorption coefficient according to [1].

[1] -Sten Seifert and Patrick Runge, "Revised refractive index and absorption of In1-xGaxAsyP1-y lattice-matched to InP in transparent and absorption IR-region," Opt. Mater. Express 6, 629-639 (2016)

    """
    Eg=1.35-0.72*y+0.12*y**2

    R=-0.00115+0.00191*Eg
    Gamma=-0.000691+0.00433*Eg
    A=-0.0453+2.1103*Eg
    a=72.32+12.78*Eg
    b=4.84+4.66*Eg
    c=-0.015+0.02*Eg
    d=-0.178+1.042*Eg
    
    
    def cot(x):
        return 1/np.tan(x)

    n=np.sqrt(1+
              a/(b-(E+1j*Gamma)**2)+
              A*np.sqrt(R)/(E+1j*Gamma)**2*(
                                            np.log(Eg**2/(Eg**2-(E+1j*Gamma)**2))+np.pi*(
                                                                                        2*cot(np.pi*np.sqrt(R/Eg))
                                                                                        -cot(np.pi*np.sqrt(R/(Eg-(E+1j*Gamma))))
                                                                                        -cot(np.pi*np.sqrt(R/(Eg+(E+1j*Gamma)))))))
    
    if type(E) == np.ndarray and type(y) !=np.ndarray:
        idx=np.where(E<Eg)
        k=np.sqrt(n**2+a/(b-E**2)).imag+c*(E-Eg)+d*(E-Eg)**2
        k[idx]=np.sqrt(n[idx]**2+a/(b-E[idx]**2)).imag
        
    elif type(E) != np.ndarray and type(y) !=np.ndarray:
        if E<Eg:
            c=d=0
        
        k=np.sqrt(n**2+a/(b-E**2)).imag+c*(E-Eg)+d*(E-Eg)**2
        
    elif type(E) != np.ndarray and type(y) ==np.ndarray:
        idx=np.where(E<Eg)
        c[idx]=0
        d[idx]=0
        
        k=np.sqrt(n**2+a/(b-E**2)).imag+c*(E-Eg)+d*(E-Eg)**2
        
    # print(n)
    return n.real + 1j*k


def get_broadband_index(model, wl_min=1500, wl_max=1600, N_samples = 100):
    '''
    Returns the complex refractive index of the model from wl_min to wl_max.
    
    model: InGaAsP model
    wl_min: minimum wavelength in nanometer;
    wl_max: maximum wavelength in nanometers;
    N_samples: number of samples in the refractive index;
    
    Returns:
    n: complex refractive index. dtype of np.complex128, of shape (N_samples, 2)
    
    '''
    
    wl = np.linspace(wl_min,wl_max, N_samples) * model.reg.nanometer
    E = 2*np.pi*model.hbar*model.c/wl

    eps_s = model.get_eps_s(wl=wl)
    alpha = model.get_alpha_sqrt(E=E)
    alpha = np.nan_to_num(alpha, nan=0)

    n0 = np.sqrt(eps_s)
    k0 = (alpha/2/(2*np.pi/wl)).to(model.reg.dimensionless)
    
    if type(model) == n_InGaAsP:
   
        dalpha_BF = model.get_dalpha_BF(E=E)
        dalpha_plasma = model.get_dalpha_plasma(E=E)

        dk_BF = (dalpha_BF/2/(2*np.pi/wl)).to(model.reg.dimensionless)
        dk_plasma = (dalpha_plasma/2/(2*np.pi/wl)).to(model.reg.dimensionless)

        dn_BF = model.get_dn_BF(E=E)
        dn_plasma = model.get_dn_plasma(E=E)

        n = (n0+dn_BF+dn_plasma) - 1j*(k0+dk_BF+dk_plasma)
        
    elif type(model) == p_InGaAsP:
        
        dalpha_BF = model.get_dalpha_BF(E=E)
        dalpha_iv = model.get_dalpha_iv(E=E)

        dk_BF = (dalpha_BF/2/(2*np.pi/wl)).to(model.reg.dimensionless)
        dk_iv = (dalpha_iv/2/(2*np.pi/wl)).to(model.reg.dimensionless)

        dn_BF = model.get_dn_BF(E=E)
        dn_iv = model.get_dn_iv(E=E)

        n = (n0+dn_BF+dn_iv) - 1j*(k0+dk_BF+dk_iv)
        
    else:
        print("WHAT THE HELL DID YOU FEED INTO THIS??")
        
    out = np.zeros((N_samples, 2), dtype = np.complex128)

    out[:, 0] = (model.c/wl).to(model.reg.hertz).magnitude
    out[:, 1] = (n**2).to(model.reg.dimensionless).magnitude
    
    return out

class n_InGaAsP(object):
    def __init__(self, N, T, y, wl, bandgap_model="jain"):
        """
        Base model for an n-type In_{1-x}Ga_{x}As_{y}P_{1-y}.
        The module returns the effects on refractive index and absorption

        N: Electron concentration in cm^-3. Float.
        T: Temperature in kelvin. Float.
        y: Concentration. Must be between 0 and 1. Float.
        wl: Wavelength in nanometers. Float.
        bandgap_model: Specifies the default bandgap model to use in the calculation of the BGN. See self.get_bandgap() for more details. String.


        References:
        1) http://www.ioffe.ru/SVA/NSM/Semicond/InP/bandstr.html#Masses
        2) Bennett, B.R., R.A. Soref, and J.A. Del Alamo. “Carrier-Induced Change in Refractive Index of InP, GaAs and InGaAsP.” IEEE Journal of Quantum Electronics 26, no. 1 (January 1990): 113–22. https://doi.org/10.1109/3.44924.
        3) Adachi, Sadao. Properties of Group-IV, III-V and II-VI Semiconductors. Wiley Series in Materials for Electronic and Optoelectronic Applications. Chichester, West Sussex, England: John Wiley & Sons, Ltd, 2006.
        4) Sze, S. M., and Kwok Kwok Ng. Physics of Semiconductor Devices. 3rd ed. Hoboken, N.J: Wiley-Interscience, 2007.
        5) Fiedler, F., and A. Schlachetzki. “Optical Parameters of InP-Based Waveguides.” Solid-State Electronics 30, no. 1 (January 1987): 73–83. https://doi.org/10.1016/0038-1101(87)90032-3.
        6) Moss, T. S., Geoffrey John Burrell, and Brian Ellis. Semiconductor Opto-Electronics. London: Butterworths, 1973.

        """
        self.reg = UnitRegistry()

        self.wl = wl * self.reg.nanometer
        self.y = y

        self.x = self.y / (2.2020 - 0.0659 * self.y)  # Taken from [5]

        self.N = N * self.reg.centimeter**-3
        self.T = T * self.reg.kelvin
        self.bandgap_model = bandgap_model

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

        # Formula from [4]
        self.ni = np.sqrt(self.Nc * self.Nv) * np.exp(
            -self.get_bandgap(model="none") / (2 * self.kb * self.T)
        )

        if self.N < self.ni:
            self.N = self.ni

        self.P = self.ni**2 / self.N  # Assumes non degenerate semiconductor

        self.so = (
            0.119 + 0.30 * self.y - 0.107 * self.y**2
        ) * self.reg.eV  # eV. Taken from [5]

        self.eps_s = self.get_eps_s(self.y, self.wl, bandgap_model=self.bandgap_model)

        # Taken from [2]
        self.C = (
            4.4e12
            * self.reg.centimeter**-1
            * self.reg.second**-0.5
            * np.sqrt(self.hbar)
        )  # Taken from Bennet 1990
        # The sqrt(hbar) comes from the fact that C comes from an earlier paper that fits an absorption curve to frequency rather than energy

        # Parameters needed for absorption adjustment
        self.a = 1
        self.r = 1

        # Adapted using formulas from [5] and

        mr_InP_hh = (1 / 0.07 + 1 / 0.6) ** -1 * self.m0
        mr_InP_hl = (1 / 0.07 + 1 / 0.12) ** -1 * self.m0
        mr_hh = (1 / self.me + 1 / self.mhh) ** -1
        mr_hl = (1 / self.me + 1 / self.mhl) ** -1

        n0_InP = get_n(self.energy.magnitude, y=0).real
        n0 = get_n(self.energy.magnitude, y=y).real
        self.n0 = n0
        
        self.Chh = self.C * (mr_InP_hh**1.5 / (mr_InP_hh**1.5 + mr_InP_hl**1.5))
        self.Chl = self.C * (mr_InP_hl**1.5 / (mr_InP_hh**1.5 + mr_InP_hl**1.5))
        # print('n', self.Chh, self.Chl)

        self.Chh = (mr_hh / mr_InP_hh) ** 1.5 * n0_InP / n0 * self.Chh
        self.Chl = (mr_hl / mr_InP_hh) ** 1.5 * n0_InP / n0 * self.Chl

        self.mue, self.muh = self.get_mobility()

        self.Ef = self.get_fermi_level()

        # Parameters for piezo effects. Taken from [3]
        self.S11 = (
            1.639e-12 * self.reg.centimeter**2 * self.reg.dyne**-1
        )  # mechanical compliance
        self.S12 = -0.589e-12 * self.reg.centimeter**2 * self.reg.dyne**-1
        self.S44 = 2.26e-12 * self.reg.centimeter**2 * self.reg.dyne**-1
        self.e14 = (
            -0.083 * self.reg.coulomb * self.reg.meter**-2
        )  # piezoelectric stress constant.

        #Density
        self.rho = (4.81+0.74*y) * self.reg.g*self.reg.centimeter**-3  # g cm^-3
        self.s = (5.2-0.372*y-0.144*y**2)*1e5 * self.reg.cm/self.reg.second  # cm s^-1 

        #Energy transitions
        self.E10 = (0.61+0.182*y+0.105*y**2) * self.reg.eV
        self.E20 = (3.38+0.549*y-0.208*y**2) * self.reg.eV

        #Phonon energies
        self.Eac = (24-2.84*y+1.57*y**2) * self.reg.meV
        self.Eop = (42.6-21.1*y+2.87*y**2) * self.reg.meV

        #Deformation potential
        self.Edef = (7.95-2.04*y+0.839*y**2) * self.reg.eV

    def get_eps_s(self, y=None, wl=None, N=None, bandgap_model=None, regime="optical"):
        """
        Returns the relative permeability considering only the bandgap narrowing effect.
        It uses the modified single oscillator model and the conversion formulas from the single oscilator model formula as outlined in [1] (Appendix).

        if regime is 'optical' it calculates the optical permeability.
        if regime is 'static' it calculates the DC permeability
        if regime is 'high frequency' it calculates the high frequency permeability.

        References:
        1) Fiedler, F., and A. Schlachetzki. “Optical Parameters of InP-Based Waveguides.” Solid-State Electronics 30, no. 1 (January 1987): 73–83. https://doi.org/10.1016/0038-1101(87)90032-3.

        """

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        if N is None:
            N = self.N

        if y is None:
            y = self.y

        if wl is None:
            wl = self.wl

        #         A = (7.255 + 1.15*y + 0.489*y**2)
        #         B = (2.316 + 0.604*y - 0.493*y**2)
        #         C = (0.3922 + 0.396*y + 0.158*y**2) * self.reg.angstrom**2

        #         return (A + B * wl**2 / (wl**2 - C*1e8)).to(self.reg.dimensionless)
        x = y / (2.2020 - 0.0659 * y)

        if regime == "optical":
            Ed = (28.91 - 9.278 * y + 5.626 * y**2) * self.reg.eV
            E0 = (3.391 - 1.652 * y + 0.863 * y**2 - 0.123 * y**3) * self.reg.eV
            Eg = self.get_bandgap(x=x, y=y, model=bandgap_model, N=N)
            E = (2 * np.pi / wl * self.hbar * self.c).to(self.reg.eV)

            a1 = Ed / E0
            a2 = Ed * E**2 / E0**3
            a3 = (
                Ed
                / (2 * E0**3 * (E0**2 - Eg**2))
                * E**4
                * np.log((2 * E0**2 - Eg**2 - E**2) / (Eg**2 - E**2))
            )
            # print(a1.to(self.reg.dimensionless))
            # print(a2.to(self.reg.dimensionless))
            # print(a3.to(self.reg.dimensionless))
            n_sq = 1 + a1 + a2 + a3
            n_sq = n_sq.to(self.reg.dimensionless)

            return n_sq

        elif regime == "static":
            return 12.35 + 1.62 * y - 0.055 * y**2

        elif regime == "high - frequency":
            return 9.52 + 2.06 * y - 0.205 * y**2

    #     def get_mobility_InGaAs(self, x):
    #         """
    #         Returns the electron and hole mobility for In_{1-x}Ga_{x}As based on an interpolation scheme reported by [1]

    #         References:
    #         1- Sotoodeh, M., A. H. Khalid, and A. A. Rezazadeh. “Empirical Low-Field Mobility Model for III–V Compounds Applicable in Device Simulation Codes.” Journal of Applied Physics 87, no. 6 (March 15, 2000): 2890–2900. https://doi.org/10.1063/1.372274.

    #         """
    #         # InAs In_0.53Ga_0.47As GaAs
    #         x_values=np.asarray([0,0.47,1])
    #         values={'mu_max': {'n':[3400, 14000, 9400], 'p':[530, 320, 491.5]},
    #                 'mu_min': {'n':[1000, 300, 500], 'p':[20,10,20]},
    #                 'Nref': {'n':[1.1e18, 1.3e17, 6.0e16], 'p':[1.1e17, 4.9e17, 1.48e17]},
    #                 'lambda': {'n':[0.32, 0.48, 0.394], 'p':[0.46, 0.403, 0.38]},
    #                 'theta1': {'n':[1.57, 1.59, 2.1], 'p':[2.3, 2.2, 1.59]},
    #                 'theta2': {'n':[3.0, 3.68, 3.0], 'p':[3.0, 3.0, 3.0]}}
    #         for key1 in values.keys():
    #             for key2 in values[key1].keys():
    #                 values[key1][key2]=np.asarray(values[key1][key2])

    #         for key in values.keys():
    #             if key != 'Nref':
    #                 values[key]['n_out']=make_interp_spline(x_values, values[key]['n'], k=2)(x)
    #                 values[key]['p_out']=make_interp_spline(x_values, values[key]['p'], k=2)(x)
    #             else:
    #                 values[key]['n_out']=10**make_interp_spline(x_values, np.log10(values[key]['n']), k=2)(x)
    #                 values[key]['p_out']=10**make_interp_spline(x_values, np.log10(values[key]['p']), k=2)(x)

    #         T=self.T.to(self.reg.kelvin).magnitude
    #         N=self.N.to(self.reg.centimeter**-3).magnitude
    #         P=self.P.to(self.reg.centimeter**-3).magnitude

    #         mobility_n=values['mu_min']['n_out'] + (values['mu_max']['n_out'] * (300/T)**values['theta1']['n_out']-values['mu_min']['n_out'])/(1+(N/values['Nref']['n_out']*(300/T)**values['theta2']['n_out'])**values['lambda']['n_out'])
    #         mobility_p=values['mu_min']['p_out'] + (values['mu_max']['p_out'] * (300/T)**values['theta1']['p_out']-values['mu_min']['p_out'])/(1+(P/values['Nref']['p_out']*(300/T)**values['theta2']['p_out'])**values['lambda']['p_out'])

    #         return (mobility_n * self.reg.centimeter**2/(self.reg.volt*self.reg.second),
    #                 mobility_p * self.reg.centimeter**2/(self.reg.volt*self.reg.second))

    #     def get_mobility_InGaP(self, x):
    #         """
    #         Returns the electron and hole mobility for In_{1-x}Ga_{x}P based on an interpolation scheme reported by [1]

    #         References:
    #         1- Sotoodeh, M., A. H. Khalid, and A. A. Rezazadeh. “Empirical Low-Field Mobility Model for III–V Compounds Applicable in Device Simulation Codes.” Journal of Applied Physics 87, no. 6 (March 15, 2000): 2890–2900. https://doi.org/10.1063/1.372274.

    #         """

    #         # InP In_0.49Ga_0.51P GaP
    #         x_values=np.asarray([0,0.51,1])
    #         values={'mu_max': {'n':[5200, 4300, 152], 'p':[170, 150, 147]},
    #                 'mu_min': {'n':[400, 400, 10], 'p':[10,15,10]},
    #                 'Nref': {'n':[3.0e17, 2.0e16, 4.4e18], 'p':[4.87e17, 1.5e17, 1.0e18]},
    #                 'lambda': {'n':[0.47, 0.70, 0.85], 'p':[0.62, 0.8, 0.85]},
    #                 'theta1': {'n':[2.0, 1.66, 1.60], 'p':[2.0, 2.0, 1.98]},
    #                 'theta2': {'n':[3.25, 0.71], 'p':[3.0, 0]}}

    #         for key1 in values.keys():
    #             for key2 in values[key1].keys():
    #                 values[key1][key2]=np.asarray(values[key1][key2])

    #         for key in values.keys():
    #             if key not in ['Nref', 'theta2']:
    #                 values[key]['n_out']=make_interp_spline(x_values, values[key]['n'], k=2)(x)
    #                 values[key]['p_out']=make_interp_spline(x_values, values[key]['p'], k=2)(x)
    #             elif key=='Nref':
    #                 values[key]['n_out']=10**make_interp_spline(x_values, np.log10(values[key]['n']), k=2)(x)
    #                 values[key]['p_out']=10**make_interp_spline(x_values, np.log10(values[key]['p']), k=2)(x)
    #             else:
    #                 values[key]['n_out']=make_interp_spline([0,1], values[key]['n'], k=1)(x)
    #                 values[key]['p_out']=make_interp_spline([0,1], values[key]['p'], k=1)(x)

    #         T=self.T.to(self.reg.kelvin).magnitude
    #         N=self.N.to(self.reg.centimeter**-3).magnitude
    #         P=self.P.to(self.reg.centimeter**-3).magnitude

    #         mobility_n=values['mu_min']['n_out'] + (values['mu_max']['n_out'] * (300/T)**values['theta1']['n_out']-values['mu_min']['n_out'])/(1+(N/values['Nref']['n_out']*(300/T)**values['theta2']['n_out'])**values['lambda']['n_out'])
    #         mobility_p=values['mu_min']['p_out'] + (values['mu_max']['p_out'] * (300/T)**values['theta1']['p_out']-values['mu_min']['p_out'])/(1+(P/values['Nref']['p_out']*(300/T)**values['theta2']['p_out'])**values['lambda']['p_out'])

    #         return (mobility_n * self.reg.centimeter**2/(self.reg.volt*self.reg.second),
    #                 mobility_p * self.reg.centimeter**2/(self.reg.volt*self.reg.second))

    def get_mobility(self, x=None, y=None):
        """
        Returns the electron and hole mobility for In_{1-x}Ga_{x}As_{y}P_{1-y} based on an interpolation scheme reported by [1]

        References:
        1- Sotoodeh, M., A. H. Khalid, and A. A. Rezazadeh. “Empirical Low-Field Mobility Model for III–V Compounds Applicable in Device Simulation Codes.” Journal of Applied Physics 87, no. 6 (March 15, 2000): 2890–2900. https://doi.org/10.1063/1.372274.

        """

        if x is None:
            x = self.x
        if y is None:
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

    def get_bandgap(self, x=None, y=None, N=None, T=None, model=None, get_BGN=False):
        """
        Return the bandgap energy of InP.
        Temperature dependence was taken from [1] page 121.

        BGN is calculated for uncompensated InGaAsP.



        if model='jain' it makes use of the model from [2]. The BGN adaptation for InGaAsP is based on the Vegard's law taken from [5].
        if model='bennet' it makes use of the model from [3]. The BGN adaptation for InGaAsP is done via the fundamental constants alone.

        References:
        1) Adachi, Sadao. Properties of Group-IV, III-V and II-VI Semiconductors. Wiley Series in Materials for Electronic and Optoelectronic Applications. Chichester, West Sussex, England: John Wiley & Sons, Ltd, 2006.
        2) Jain, S. C., J. M. McGregor, and D. J. Roulston. “Band‐gap Narrowing in Novel III‐V Semiconductors.” Journal of Applied Physics 68, no. 7 (October 1990): 3747–49. https://doi.org/10.1063/1.346291.
        3) http://www.ioffe.ru/SVA/NSM/Semicond/InP
        4) Fiedler, F., and A. Schlachetzki. “Optical Parameters of InP-Based Waveguides.” Solid-State Electronics 30, no. 1 (January 1987): 73–83. https://doi.org/10.1016/0038-1101(87)90032-3.

        """

        if model is None:
            model = self.bandgap_model

        if N is None:
            N = self.N.to(self.reg.centimeter**-3)

        if T is None:
            T = self.T

        if y is None:
            y = self.y

        if x is None:
            x = self.x

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

            BGN = (
                A * N.magnitude ** (1 / 3)
                + B * N.magnitude**0.25
                + C * N.magnitude**0.5
            ) * self.reg.eV

        elif model == "none":
            BGN = 0

        # Eg formula taken from [5].
        if not get_BGN:
            Eg0 = 1.35 - 0.72 * y + 0.12 * y**2

            return Eg0 * self.reg.eV - BGN
        else:
            return BGN

    def get_fermi_level(self, N=None, T=None, bandgap_model=None):
        """
        Returns the fermi level calculated assuming degenerate semiconductors. It uses equation 26a from [1] which is the known Joyce-Dixon approximation. Assume Ec=0.
        It assumes that only uncompensated materials are given.
        It assumes full ionization and the carrier concentration.


        References:
        1) Sze, S. M., and Kwok Kwok Ng. Physics of Semiconductor Devices. 3rd ed. Hoboken, N.J: Wiley-Interscience, 2007.

        """

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        if N is None:
            N = self.N
        if T is None:
            T = self.T

        Eg = self.get_bandgap(N=N, T=T, model=bandgap_model)

        Ef = (np.log(N / self.Nc) + 1 / np.sqrt(8) * N / self.Nc) * self.kb * self.T

        return Ef.to(self.reg.eV)

    def FD(self, E, Ef=1):
        return 1 / (1 + np.exp((E - Ef) / (self.kb * self.T)))

    def get_dalpha_BF(self, E=None, bandgap_model=None):
        """
        Gives the change in absorption due to the bandfilling effect.
        The treatment if based on [1] but we neglect the quasi-fermi levels and use the fermi level numerically calculated instead.

        Note that if the bandgap_model is not 'none' the BGN effect is taken into consideration within the bandfilling.

        E must be a Pint quantity in energy using the parent object's register!!

        References:
        1) Bennett, B.R., R.A. Soref, and J.A. Del Alamo. “Carrier-Induced Change in Refractive Index of InP, GaAs and InGaAsP.” IEEE Journal of Quantum Electronics 26, no. 1 (January 1990): 113–22. https://doi.org/10.1109/3.44924.

        """

        if E is None:
            E = self.energy

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        Ef = self.Ef
        Eg = self.get_bandgap(model=bandgap_model)

        Eah = (Eg - E) * (self.me / (self.me + self.mhh)) - Eg
        Eal = (Eg - E) * (self.me / (self.me + self.mhl)) - Eg
        Ebh = (E - Eg) * (self.mhh / (self.me + self.mhh))
        Ebl = (E - Eg) * (self.mhl / (self.me + self.mhl))

        alpha0 = self.Chh / E * np.sqrt(E - Eg) + self.Chl / E * np.sqrt(E - Eg)
        alpha = self.Chh / E * np.sqrt(E - Eg) * (
            self.FD(Eah, Ef=Ef) - self.FD(Ebh, Ef=Ef)
        ) + self.Chl / E * np.sqrt(E - Eg) * (self.FD(Eal, Ef=Ef) - self.FD(Ebl, Ef=Ef))
        # print(self.Chh/E*self.FD(Eah, Ef=Ef1).magnitude+self.Chl/E*self.FD(Eal, Ef=Ef1).magnitude)
        alpha = np.nan_to_num(alpha)
        alpha0 = np.nan_to_num(alpha0)

        return (alpha - alpha0).to(self.reg.centimeter**-1)

    def get_alpha_sqrt(self, E=None, bandgap_model=None):
        if E is None:
            E = self.energy

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        Eg = self.get_bandgap(model=bandgap_model)

        alpha0 = self.Chh / E * np.sqrt(E - Eg) + self.Chl / E * np.sqrt(E - Eg)

        return alpha0.to(self.reg.centimeter**-1)

    def get_dalpha_plasma(self, E=None, N=None):
        """
        Gives the change in absorption due to the plasma effect as reported in [1]. Simply put, the model is based on a second order perturbation theory and considers 3 scattering mechanisms: electron - optical phonon, electron - acoustical phonon and electron - ionized impurity. It assumes room temperature of 300K. It is limited to concentrations above 6e18 cm^-3.

        It also accounts for the interband transitions happening to other valleys in the conduction band as reported in [3] through eq. 4.4.

        E must be a Pint quantity in energy using the parent object's register!!

        The values predicted by the below formula are consistent with the experimental results from [2].

        References:
        1) Walukiewicz, W., J. Lagowski, L. Jastrzebski, P. Rava, M. Lichtensteiger, C. H. Gatos, and H. C. Gatos. “Electron Mobility and Free-Carrier Absorption in InP; Determination of the Compensation Ratio.” Journal of Applied Physics 51, no. 5 (May 1, 1980): 2659–68. https://doi.org/10.1063/1.327925.

        2) Dumke, W. P., M. R. Lorenz, and G. D. Pettit. “Intra- and Interband Free-Carrier Absorption and the Fundamental Absorption Edge in n -Type InP.” Physical Review B 1, no. 12 (June 15, 1970): 4668–73. https://doi.org/10.1103/PhysRevB.1.4668.

        3) Fiedler, F. and Schlachetzki, A. (1987) “Optical parameters of InP-based waveguides,” Solid-State Electronics, 30(1), pp. 73–83. Available at: https://doi.org/10.1016/0038-1101(87)90032-3.


        """

        if N is None:
            N = self.N

        if E is None:
            E = self.energy

        # Account for limitation on doping
        if type(N) == np.ndarray:
            idx = np.where(N.to(self.reg.centimeter**-3).magnitude > 6e18)
            N[idx] = 6e18 * self.reg.centimeter**-3
        else:
            if N.to(self.reg.centimeter**-3).magnitude > 6e18:
                N = 6e18 * self.reg.centimeter**-3

            if N.to(self.reg.centimeter**-3).magnitude < 1e16:
                N = 1e16 * self.reg.centimeter**-3

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

        doping = N.to(self.reg.centimeter**-3).magnitude

        alpha = (
            alpha_imp_interp(doping) * (wv_ratio) ** 3.5
            + alpha_op_interp(doping) * (wv_ratio) ** 2.5
            + alpha_ac_interp(doping) * (wv_ratio) ** 1.5
        )

        ## Now we must also account for the interband transitions 

        # A = (1.4+1.85*self.y)*1e-5*1.3e23 * self.reg.eV**-1 * self.reg.g * self.reg.s**-2*self.reg.cm**-2

        A = (1.4+1.85*self.y)*1e-5 * (
            self.me.units**1.5 *
            self.Eac.units * 
            self.Edef.units *
            1/self.rho.units *
            1/self.s.units**2 *
            1/self.reg.eV**3 *
            self.reg.eV**2
        )**-1 * self.reg.cm**-1 * 1.3e65 #The factor of is necessary to make the calculations match the ones from fiedler paper


        if type(E.magnitude) == np.ndarray:
            tmp1 = (E + self.Eac - self.E10).to(self.reg.eV).magnitude
            idx_pos = np.where(tmp1 > 0)[0]
            idx_neg = np.where(tmp1 <= 0)[0]

            u_plus = np.zeros(E.shape) * self.reg.eV
            u_plus[idx_neg] = self.E10 - E[idx_neg] - self.Eac

            tmp2 = (E - self.Eac - self.E10).to(self.reg.eV).magnitude
            idx_pos = np.where(tmp2 > 0)[0]
            idx_neg = np.where(tmp2 <= 0)[0]
            u_minus = np.zeros(E.shape) * self.reg.eV
            u_minus[idx_neg] = self.E10 - E[idx_neg] + self.Eac

        else:
            if (E + self.Eac - self.E10).to(self.reg.eV).magnitude > 0:
                u_plus = 0 * self.reg.eV
            else:
                u_plus = self.E10 - E - self.Eac

            if (E - self.Eac - self.E10).to(self.reg.eV).magnitude > 0:
                u_minus = 0 * self.reg.eV
            else:
                u_minus = self.E10 - E + self.Eac


        a1 = A*(self.me)**1.5*self.Eac*self.Edef/self.n0/self.rho/self.s**2/(np.exp(self.Eac/(self.kb*self.T))-1)
        a2 = 1/((self.E20 - E)**2 * E)
        a3 = np.exp(self.Eac/self.kb/self.T)

        energy_integrand_a4 = np.linspace(u_plus.to(self.reg.eV).magnitude, 100, 1000) * self.reg.eV
        a4_integrand = energy_integrand_a4**0.5 * (energy_integrand_a4-self.E10+E+self.Eac)**0.5/(np.exp((energy_integrand_a4-self.Ef)/self.kb/self.T)+1)
        a4_integrand = np.nan_to_num(a4_integrand) #This removes nan values that stem from very small negative numbers inside the root like -1e-17
        a4 = simpson(a4_integrand.to(self.reg.eV).magnitude, x=energy_integrand_a4.to(self.reg.eV).magnitude, axis=0) * self.reg.eV**2


        energy_integrand_a5 = np.linspace(u_minus.to(self.reg.eV).magnitude, 100, 10000) * self.reg.eV
        a5_integrand = energy_integrand_a5**0.5 * (energy_integrand_a5-self.E10+E+self.Eac)**0.5/(np.exp((energy_integrand_a5-self.Ef)/self.kb/self.T)+1)
        a5_integrand = np.nan_to_num(a5_integrand) #This removes nan values that stem from very small negative numbers inside the root like -1e-17
        a5 = simpson(a5_integrand.to(self.reg.eV).magnitude, x=energy_integrand_a5.to(self.reg.eV).magnitude, axis=0) * self.reg.eV**2

        alpha_IB = a1*a2*(a3*a4+a5)

        Eg = self.get_bandgap(model='jain')
        alpha_VC = (3e3 * np.exp(-100*(Eg-E).to(self.reg.eV).magnitude)) * self.reg.centimeter**-1

        return (alpha) * self.reg.centimeter**-1 + alpha_VC + alpha_IB 

    def get_dn_BF(self, E=None, bandgap_model=None, h=1e-3):
        """
        Returns the change in refractive index based only on the band filling effect based on the kramers kronig relations.

        E must be a Pint quantity in energy using the parent object's register!!

        h: float. energy step in eV.

        The algorithm is an adaptation of the Mclaurin algorithm described in [1].

        References:
        [1] Ohta, Koji, and Hatsuo Ishida. “Comparison among Several Numerical Integration Methods for Kramers-Kronig Transformation.” Applied Spectroscopy 42, no. 6 (August 1988): 952–57. https://doi.org/10.1366/0003702884430380.

        """

        if E is None:
            E = self.energy

        if type(E.magnitude) is not np.ndarray:
            E = np.asarray([E.magnitude]) * E.units

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        E = E.to(self.reg.eV)
        h = h * self.reg.eV

        Eg = self.get_bandgap(
            model="none"
        )  # this is left as 'none' because it is only used for the limits of integration. That way they are always the same.

        n_points_right = ((2 * Eg - E) / h).to(self.reg.dimensionless).astype(int)
        n_points_left = (
            ((E - 0.001 * self.reg.eV) / h).to(self.reg.dimensionless).astype(int)
        )
        Nright = np.max(n_points_right)
        Nleft = np.max(n_points_left)
        # print(Nright, Nleft)

        energies_right = E[..., None] + h * np.arange(Nright)[None, ...]
        energies_left = E[..., None] - h * (1 + np.arange(Nleft)[None, ...])

        energies = np.concatenate((energies_left, energies_right), axis=1)
        # print(energies.shape, Nleft, Nright, Nleft+Nright)
        if Nleft % 2 == 0:
            energies_integrand = energies[:, 1::2]
        else:
            energies_integrand = energies[:, 0::2]

        alpha = self.get_dalpha_BF(energies_integrand, bandgap_model=bandgap_model)

        integral = (
            np.sum(
                alpha
                * (
                    1 / (energies_integrand - E[..., None])
                    + 1 / (energies_integrand + E[..., None])
                )
                * 1
                / (2 * energies_integrand),
                axis=1,
            )
            * 2
            * h
        )

        return np.squeeze((integral * self.c * self.hbar / np.pi)).to(
            self.reg.dimensionless
        )  # the reason why you dont divide by e is because of the eV dimension on the integral result!!

    def get_dn_plasma(self, E=None, N=None):
        """
        Returns the change in refractive index based on the plasma effect. It makes use of [1], page 79, eq.5.2.14

        E must be a Pint quantity in energy using the parent object's register!!

        References:
        1) Hunsperger, Robert G. Integrated Optics: Theory and Technology. 6th ed. Advanced Texts in Physics. New York London: Springer, 2009.

        """

        if N is None:
            N = self.N

        if E is None:
            E = self.energy

        return (
            -1
            / 2
            * (
                N
                * self.e**2
                / (self.me * self.e0 * E**2 / self.hbar**2 * np.sqrt(self.eps_s))
            ).to(self.reg.dimensionless)
        )

    def get_dperm_pockels(self, E, Efield, bandgap_model=None):
        """
        Returns the change in refractive index owing to the pockels effect.

        E: energy of the field's photons. Pint quantity
        Efield: Electric field components. It is assumed to have dimensions (3,Ny, Nx)

        References:
        1) Adachi, Sadao, and Kunishige Oe. “Internal Strain and Photoelastic Effects in Ga  1− x  Al  x  As/GaAs and In  1− x  Ga  x  As  y  P  1− y  /InP Crystals.” Journal of Applied Physics 54, no. 11 (November 1983): 6620–27. https://doi.org/10.1063/1.331898.
        2) Adachi, Sadao, and Kunishige Oe. “Linear Electro‐optic Effects in Zincblende‐type Semiconductors: Key Properties of InGaAsP Relevant to Device Design.” Journal of Applied Physics 56, no. 1 (July 1984): 74–80. https://doi.org/10.1063/1.333731.

        """

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        def g(chi):
            return 1 / chi**2 * (2 - (1 + chi) ** -0.5 - (1 - chi) ** -0.5)

        def f(chi):
            return 1 / chi**2 * (2 - (1 + chi) ** 0.5 - (1 - chi) ** 0.5)

        def symmetric(x):
            """
            Returns an array in the standard symmetric form instead of the voigt notation. Assumes shape of [6, Ny, Nx]

            """
            y = np.zeros((3, 3) + x.shape[-2:]) * x.units

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

        Eg = self.get_bandgap()

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

        pockels_tensor = np.zeros((6, 3)) * r41_free / r41_free.magnitude

        pockels_tensor[3, 0] = r41_free + r41_piezo
        pockels_tensor[4, 1] = r41_free + r41_piezo
        pockels_tensor[5, 2] = r41_free + r41_piezo

        # Evaluate the change in impermeability
        deta = np.einsum("il,ljk->ijk", pockels_tensor, Efield)
        deta = symmetric(deta)  # return to symmetric matrix

        # build the permitivity tensor
        perm = np.zeros((3, 3)) * self.e0.units
        perm[0, 0] = self.eps_s * self.e0
        perm[1, 1] = self.eps_s * self.e0
        perm[2, 2] = self.eps_s * self.e0

        # find the change in permitivity
        dperm = np.einsum("ik,kltp,lj->ijtp", perm, deta, perm) * -1 / self.e0

        return dperm

    def get_pockels_coeffs(self, E, bandgap_model=None):
        """
        Returns the pockels coefficients

        E: energy of the field's photons. Pint quantity
        Efield: Electric field components. It is assumed to have dimensions (3,Ny, Nx)

        References:
        1) Adachi, Sadao, and Kunishige Oe. “Internal Strain and Photoelastic Effects in Ga  1− x  Al  x  As/GaAs and In  1− x  Ga  x  As  y  P  1− y  /InP Crystals.” Journal of Applied Physics 54, no. 11 (November 1983): 6620–27. https://doi.org/10.1063/1.331898.
        2) Adachi, Sadao, and Kunishige Oe. “Linear Electro‐optic Effects in Zincblende‐type Semiconductors: Key Properties of InGaAsP Relevant to Device Design.” Journal of Applied Physics 56, no. 1 (July 1984): 74–80. https://doi.org/10.1063/1.333731.

        """

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        def g(chi):
            return 1 / chi**2 * (2 - (1 + chi) ** -0.5 - (1 - chi) ** -0.5)

        def f(chi):
            return 1 / chi**2 * (2 - (1 + chi) ** 0.5 - (1 - chi) ** 0.5)

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
        # print(E0, F0, C, D)
        Eg = self.get_bandgap(model=bandgap_model)

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
            ).to(self.reg.centimeter**2 / self.reg.dyne)
            * self.e14
        )

        return r41_free, r41_piezo

    def get_dperm_kerr(self, E, Efield, bandgap_model=None):
        """
        Returns the change in refractive index owing to the Kerr effect.

        E: energy of the field's photons. Pint quantity
        Efield: Electric field components. It is assumed to have dimensions (3,Ny, Nx)

        References:
        [1] - Maat, Derk Hendrik Pieter. InP-Based Integrated MZI Switches for Optical Communication, 2001.


        """

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        def symmetric(x):
            """
            Returns an array in the standard symmetric form instead of the voigt notation. Assumes shape of [6, Ny, Nx]

            """
            y = np.zeros((3, 3) + x.shape[-2:]) * x.units

            y[0, 0] = x[0]
            y[0, 1] = y[1, 0] = x[5]
            y[0, 2] = y[2, 0] = x[4]
            y[1, 1] = x[1]
            y[1, 2] = y[2, 1] = x[3]
            y[2, 2] = x[2]

            return y

        # Taken from [1]
        # A_TE = 0.25e3
        # A_TM = 0.20e3
        # B_TE = 0.71e9
        # B_TM = 0.48e9

        #These allow the replication of the results of [1] while using the Imodulator.
        A_TE = 0.9e3
        A_TM = 1.7e3
        B_TE = 0.42e9
        B_TM = 0.42e9

        # C_TE=-1.79e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2
        # C_TM=-1.82e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2

        C_TE = -3.10e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2
        C_TM = -5.60e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2

        A_mat = (
            np.asarray([[A_TE, 0, 0], [0, A_TM, 0], [0, 0, A_TE]])
            * self.reg.eV
            / (self.reg.volt * self.reg.meter)
        )

        B_mat = (
            np.asarray([[B_TE, 0, 0], [0, B_TM, 0], [0, 0, B_TE]])
            * self.reg.eV ** (-3 / 2)
            * self.reg.volt
            / self.reg.meter
        )

        freq = E / self.hbar / (2 * np.pi)
        wl = self.c / freq

        Eg = self.get_bandgap(model=bandgap_model)

        # print((2*np.pi*self.c/(Eg/self.hbar)).to(self.reg.nanometer))
        Efield_mag_sq = np.einsum("ijk,ijk -> jk", Efield.conjugate(), Efield).real

        dalpha = (
            A_mat[..., None, None]
            * wl
            * np.sqrt(Efield_mag_sq)[None, None, ...]
            / (Eg - E)
            * 10
            ** (
                -(
                    B_mat[..., None, None]
                    * (Eg - E) ** (3 / 2)
                    / np.sqrt(Efield_mag_sq)[None, None, ...]
                )
            )
        )
        dalpha = dalpha.to(self.reg.meter**-1)

        dperm_imag = (
            np.sqrt(self.eps_s) * self.c / (2 * np.pi * freq) * dalpha * self.e0
        )
        # print(type(E.magnitude))
        if type(E.magnitude) == np.ndarray:
            S11 = C_TE * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
            S12 = C_TM * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
            idx = np.where(E > Eg)
            S11[idx] = np.nan * self.reg.meter**2 * self.reg.volt**-2
            S12[idx] = np.nan * self.reg.meter**2 * self.reg.volt**-2

        else:
            if E > Eg:
                S11 = np.nan * self.reg.meter**2 * self.reg.volt**-2
                S12 = np.nan * self.reg.meter**2 * self.reg.volt**-2
            else:
                S11 = C_TE * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
                S12 = C_TM * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)

        S11 = S11.to(self.reg.meter**2 / self.reg.volt**2).magnitude
        S12 = S12.to(self.reg.meter**2 / self.reg.volt**2).magnitude
        S_mat = (
            np.asarray(
                [
                    [S11, S12, S12, 0, 0, 0],
                    [S12, S11, S12, 0, 0, 0],
                    [S12, S12, S11, 0, 0, 0],
                    [0, 0, 0, S11 - S12, 0, 0],
                    [0, 0, 0, 0, S11 - S12, 0],
                    [0, 0, 0, 0, 0, S11 - S12],
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
        deta_real = symmetric(np.einsum("ij,jkl->ikl", S_mat, Efield_voigt))

        # build the permitivity tensor
        perm = np.zeros((3, 3)) * self.e0.units
        perm[0, 0] = self.eps_s * self.e0
        perm[1, 1] = self.eps_s * self.e0
        perm[2, 2] = self.eps_s * self.e0

        # find the change in permitivity
        dperm = (
            np.einsum("ik,kltp,lj->ijtp", perm, deta_real, perm) * -1 / self.e0
            + 1j * dperm_imag
        )

        return dperm

    def get_kerr_coeffs(self, E, bandgap_model=None):
        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        # C_TE=-1.79e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2
        # C_TM=-1.82e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2

        C_TE = -3.10e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2
        C_TM = -5.60e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2

        Eg = self.get_bandgap(model=bandgap_model)

        if type(E.magnitude) == np.ndarray:
            S11 = C_TE * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
            S12 = C_TM * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
            idx = np.where(E > Eg)
            S11[idx] = np.nan * self.reg.meter**2 * self.reg.volt**-2
            S12[idx] = np.nan * self.reg.meter**2 * self.reg.volt**-2
        else:
            if E > Eg:
                S11 = np.nan * self.reg.meter**2 * self.reg.volt**-2
                S12 = np.nan * self.reg.meter**2 * self.reg.volt**-2
            else:
                S11 = C_TE * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
                S12 = C_TM * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)

        S11 = S11.to(self.reg.meter**2 / self.reg.volt**2).magnitude
        S12 = S12.to(self.reg.meter**2 / self.reg.volt**2).magnitude

        return S11, S12


class p_InGaAsP(object):
    def __init__(self, P, T, y, wl, bandgap_model="jain"):
        """
        Base model for an p-type In_{1-x}Ga_{x}As_{y}P_{1-y}.
        The module returns the effects on refractive index and absorption

        P: Electron concentration in cm^-3. Float.
        T: Temperature in kelvin. Float.
        y: Concentration. Must be between 0 and 1. Float.
        wl: Wavelength in nanometers. Float.
        bandgap_model: Specifies the default bandgap model to use in the calculation of the BGN. See self.get_bandgap() for more details. String.


        References:
        1) http://www.ioffe.ru/SVA/NSM/Semicond/InP/bandstr.html#Masses
        2) Bennett, B.R., R.A. Soref, and J.A. Del Alamo. “Carrier-Induced Change in Refractive Index of InP, GaAs and InGaAsP.” IEEE Journal of Quantum Electronics 26, no. 1 (January 1990): 113–22. https://doi.org/10.1109/3.44924.
        3) Adachi, Sadao. Properties of Group-IV, III-V and II-VI Semiconductors. Wiley Series in Materials for Electronic and Optoelectronic Applications. Chichester, West Sussex, England: John Wiley & Sons, Ltd, 2006.
        4) Sze, S. M., and Kwok Kwok Ng. Physics of Semiconductor Devices. 3rd ed. Hoboken, N.J: Wiley-Interscience, 2007.
        5) Fiedler, F., and A. Schlachetzki. “Optical Parameters of InP-Based Waveguides.” Solid-State Electronics 30, no. 1 (January 1987): 73–83. https://doi.org/10.1016/0038-1101(87)90032-3.
        6) Moss, T. S., Geoffrey John Burrell, and Brian Ellis. Semiconductor Opto-Electronics. London: Butterworths, 1973.

        """
        self.reg = UnitRegistry()

        self.wl = wl * self.reg.nanometer
        self.y = y

        self.x = self.y / (2.2020 - 0.0659 * self.y)  # Taken from [5]

        self.P = P * self.reg.centimeter**-3
        self.T = T * self.reg.kelvin
        self.bandgap_model = bandgap_model

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

        # Formula from [4]
        self.ni = np.sqrt(self.Nc * self.Nv) * np.exp(
            -self.get_bandgap(model="none") / (2 * self.kb * self.T)
        )

        if self.P < self.ni:
            self.P = self.ni

        self.N = self.ni**2 / self.P  # Assumes non degenerate semiconductor

        self.so = (
            0.119 + 0.30 * self.y - 0.107 * self.y**2
        ) * self.reg.eV  # eV. Taken from [5]

        self.eps_s = self.get_eps_s(self.y, self.wl, bandgap_model=self.bandgap_model)

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

        n0_InP = get_n(self.energy.magnitude, y=0).real
        n0 = get_n(self.energy.magnitude, y=y).real

        self.Chh = self.C * (mr_InP_hh**1.5 / (mr_InP_hh**1.5 + mr_InP_hl**1.5))
        self.Chl = self.C * (mr_InP_hl**1.5 / (mr_InP_hh**1.5 + mr_InP_hl**1.5))
        # print('p',self.Chh, self.Chl)

        self.Chh = (mr_hh / mr_InP_hh) ** 1.5 * n0_InP / n0 * self.Chh
        self.Chl = (mr_hl / mr_InP_hh) ** 1.5 * n0_InP / n0 * self.Chl

        # print('p',self.get_bandgap())

        self.mue, self.muh = self.get_mobility(self.x, self.y)

        self.Ef = self.get_fermi_level()

        # Parameters for piezo effects. Taken from [3]
        self.S11 = (
            1.639e-12 * self.reg.centimeter**2 * self.reg.dyne**-1
        )  # mechanical compliance
        self.S12 = -0.589e-12 * self.reg.centimeter**2 * self.reg.dyne**-1
        self.S44 = 2.26e-12 * self.reg.centimeter**2 * self.reg.dyne**-1
        self.e14 = (
            -0.083 * self.reg.coulomb * self.reg.meter**-2
        )  # piezoelectric stress constant.

        #Density
        self.rho = (4.81+0.74*y) * self.reg.g*self.reg.centimeter**-3  # g cm^-3
        self.s = (5.2-0.372*y-0.144*y**2)*1e5 * self.reg.cm/self.reg.second  # cm s^-1 

        #Energy transitions
        self.E10 = (0.61+0.182*y+0.105*y**2) * self.reg.eV
        self.E20 = (3.38+0.549*y-0.208*y**2) * self.reg.eV

        #Phonon energies
        self.Eac = (24-2.84*y+1.57*y**2) * self.reg.meV
        self.Eop = (42.6-21.1*y+2.87*y**2) * self.reg.meV

        #Deformation potential
        self.Edef = (7.95-2.04*y+0.839*y**2) * self.reg.eV


    def get_eps_s(self, y=None, wl=None, P=None, bandgap_model=None, regime="optical"):
        """
        Returns the relative permeability considering only the bandgap narrowing effect.
        It uses the modified single oscillator model and the conversion formulas from the single oscilator model formula as outlined in [1] (Appendix).

        if regime is 'optical' it calculates the optical permeability.
        if regime is 'static' it calculates the DC permeability
        if regime is 'high frequency' it calculates the high frequency permeability.

        References:
        1) Fiedler, F., and A. Schlachetzki. “Optical Parameters of InP-Based Waveguides.” Solid-State Electronics 30, no. 1 (January 1987): 73–83. https://doi.org/10.1016/0038-1101(87)90032-3.

        """

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        if P is None:
            P = self.P

        if y is None:
            y = self.y

        if wl is None:
            wl = self.wl

        #         A = (7.255 + 1.15*y + 0.489*y**2)
        #         B = (2.316 + 0.604*y - 0.493*y**2)
        #         C = (0.3922 + 0.396*y + 0.158*y**2) * self.reg.angstrom**2

        #         return (A + B * wl**2 / (wl**2 - C*1e8)).to(self.reg.dimensionless)
        x = y / (2.2020 - 0.0659 * y)

        if regime == "optical":
            Ed = (28.91 - 9.278 * y + 5.626 * y**2) * self.reg.eV
            E0 = (3.391 - 1.652 * y + 0.863 * y**2 - 0.123 * y**3) * self.reg.eV
            Eg = self.get_bandgap(x=x, y=y, model=bandgap_model, P=P)
            E = (2 * np.pi / wl * self.hbar * self.c).to(self.reg.eV)

            a1 = Ed / E0
            a2 = Ed * E**2 / E0**3
            a3 = (
                Ed
                / (2 * E0**3 * (E0**2 - Eg**2))
                * E**4
                * np.log((2 * E0**2 - Eg**2 - E**2) / (Eg**2 - E**2))
            )
            # print(a1.to(self.reg.dimensionless))
            # print(a2.to(self.reg.dimensionless))
            # print(a3.to(self.reg.dimensionless))
            n_sq = 1 + a1 + a2 + a3
            n_sq = n_sq.to(self.reg.dimensionless)

            return n_sq

        elif regime == "static":
            return 12.35 + 1.62 * y - 0.055 * y**2

        elif regime == "high - frequency":
            return 9.52 + 2.06 * y - 0.205 * y**2

    #     def get_mobility_InGaAs(self, x=None):
    #         """
    #         Returns the electron and hole mobility for In_{1-x}Ga_{x}As based on an interpolation scheme reported by [1]

    #         References:
    #         1- Sotoodeh, M., A. H. Khalid, and A. A. Rezazadeh. “Empirical Low-Field Mobility Model for III–V Compounds Applicable in Device Simulation Codes.” Journal of Applied Physics 87, no. 6 (March 15, 2000): 2890–2900. https://doi.org/10.1063/1.372274.

    #         """

    #         if x is None:
    #             x = self.x

    #         # InAs In_0.53Ga_0.47As GaAs
    #         x_values=np.asarray([0,0.47,1])
    #         values={'mu_max': {'n':[3400, 14000, 9400], 'p':[530, 320, 491.5]},
    #                 'mu_min': {'n':[1000, 300, 500], 'p':[20,10,20]},
    #                 'Nref': {'n':[1.1e18, 1.3e17, 6.0e16], 'p':[1.1e17, 4.9e17, 1.48e17]},
    #                 'lambda': {'n':[0.32, 0.48, 0.394], 'p':[0.46, 0.403, 0.38]},
    #                 'theta1': {'n':[1.57, 1.59, 2.1], 'p':[2.3, 2.2, 1.59]},
    #                 'theta2': {'n':[3.0, 3.68, 3.0], 'p':[3.0, 3.0, 3.0]}}
    #         for key1 in values.keys():
    #             for key2 in values[key1].keys():
    #                 values[key1][key2]=np.asarray(values[key1][key2])

    #         for key in values.keys():
    #             if key != 'Nref':
    #                 values[key]['n_out']=make_interp_spline(x_values, values[key]['n'], k=2)(x)
    #                 values[key]['p_out']=make_interp_spline(x_values, values[key]['p'], k=2)(x)
    #             else:
    #                 values[key]['n_out']=10**make_interp_spline(x_values, np.log10(values[key]['n']), k=2)(x)
    #                 values[key]['p_out']=10**make_interp_spline(x_values, np.log10(values[key]['p']), k=2)(x)

    #         T=self.T.to(self.reg.kelvin).magnitude
    #         N=self.N.to(self.reg.centimeter**-3).magnitude
    #         P=self.P.to(self.reg.centimeter**-3).magnitude

    #         mobility_n=values['mu_min']['n_out'] + (values['mu_max']['n_out'] * (300/T)**values['theta1']['n_out']-values['mu_min']['n_out'])/(1+(N/values['Nref']['n_out']*(300/T)**values['theta2']['n_out'])**values['lambda']['n_out'])
    #         mobility_p=values['mu_min']['p_out'] + (values['mu_max']['p_out'] * (300/T)**values['theta1']['p_out']-values['mu_min']['p_out'])/(1+(P/values['Nref']['p_out']*(300/T)**values['theta2']['p_out'])**values['lambda']['p_out'])

    #         return (mobility_n * self.reg.centimeter**2/(self.reg.volt*self.reg.second),
    #                 mobility_p * self.reg.centimeter**2/(self.reg.volt*self.reg.second))

    #     def get_mobility_InGaP(self, x):
    #         """
    #         Returns the electron and hole mobility for In_{1-x}Ga_{x}P based on an interpolation scheme reported by [1]

    #         References:
    #         1- Sotoodeh, M., A. H. Khalid, and A. A. Rezazadeh. “Empirical Low-Field Mobility Model for III–V Compounds Applicable in Device Simulation Codes.” Journal of Applied Physics 87, no. 6 (March 15, 2000): 2890–2900. https://doi.org/10.1063/1.372274.

    #         """

    #         # InP In_0.49Ga_0.51P GaP
    #         x_values=np.asarray([0,0.51,1])
    #         values={'mu_max': {'n':[5200, 4300, 152], 'p':[170, 150, 147]},
    #                 'mu_min': {'n':[400, 400, 10], 'p':[10,15,10]},
    #                 'Nref': {'n':[3.0e17, 2.0e16, 4.4e18], 'p':[4.87e17, 1.5e17, 1.0e18]},
    #                 'lambda': {'n':[0.47, 0.70, 0.85], 'p':[0.62, 0.8, 0.85]},
    #                 'theta1': {'n':[2.0, 1.66, 1.60], 'p':[2.0, 2.0, 1.98]},
    #                 'theta2': {'n':[3.25, 0.71], 'p':[3.0, 0]}}

    #         for key1 in values.keys():
    #             for key2 in values[key1].keys():
    #                 values[key1][key2]=np.asarray(values[key1][key2])

    #         for key in values.keys():
    #             if key not in ['Nref', 'theta2']:
    #                 values[key]['n_out']=make_interp_spline(x_values, values[key]['n'], k=2)(x)
    #                 values[key]['p_out']=make_interp_spline(x_values, values[key]['p'], k=2)(x)
    #             elif key=='Nref':
    #                 values[key]['n_out']=10**make_interp_spline(x_values, np.log10(values[key]['n']), k=2)(x)
    #                 values[key]['p_out']=10**make_interp_spline(x_values, np.log10(values[key]['p']), k=2)(x)
    #             else:
    #                 values[key]['n_out']=make_interp_spline([0,1], values[key]['n'], k=1)(x)
    #                 values[key]['p_out']=make_interp_spline([0,1], values[key]['p'], k=1)(x)

    #         T=self.T.to(self.reg.kelvin).magnitude
    #         N=self.N.to(self.reg.centimeter**-3).magnitude
    #         P=self.P.to(self.reg.centimeter**-3).magnitude

    #         mobility_n=values['mu_min']['n_out'] + (values['mu_max']['n_out'] * (300/T)**values['theta1']['n_out']-values['mu_min']['n_out'])/(1+(N/values['Nref']['n_out']*(300/T)**values['theta2']['n_out'])**values['lambda']['n_out'])
    #         mobility_p=values['mu_min']['p_out'] + (values['mu_max']['p_out'] * (300/T)**values['theta1']['p_out']-values['mu_min']['p_out'])/(1+(P/values['Nref']['p_out']*(300/T)**values['theta2']['p_out'])**values['lambda']['p_out'])

    #         return (mobility_n * self.reg.centimeter**2/(self.reg.volt*self.reg.second),
    #                 mobility_p * self.reg.centimeter**2/(self.reg.volt*self.reg.second))

    def get_mobility(self, x=None, y=None):
        """
        Returns the electron and hole mobility for In_{1-x}Ga_{x}As_{y}P_{1-y} based on an interpolation scheme reported by [1]

        References:
        1- Sotoodeh, M., A. H. Khalid, and A. A. Rezazadeh. “Empirical Low-Field Mobility Model for III–V Compounds Applicable in Device Simulation Codes.” Journal of Applied Physics 87, no. 6 (March 15, 2000): 2890–2900. https://doi.org/10.1063/1.372274.

        """

        if x is None:
            x = self.x
        if y is None:
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

    def get_bandgap(self, x=None, y=None, P=None, T=None, model=None, get_BGN=False):
        """
        Return the bandgap energy of InGaAsP.
        Temperature dependence was taken from [1] page 121.

        BGN is calculated for uncompensated InGaAsP.



        if model='jain' it makes use of the model from [2]. The BGN adaptation for InGaAsP is based on the Vegard's law taken from [5].

        References:
        1) Adachi, Sadao. Properties of Group-IV, III-V and II-VI Semiconductors. Wiley Series in Materials for Electronic and Optoelectronic Applications. Chichester, West Sussex, England: John Wiley & Sons, Ltd, 2006.
        2) Jain, S. C., J. M. McGregor, and D. J. Roulston. “Band‐gap Narrowing in Novel III‐V Semiconductors.” Journal of Applied Physics 68, no. 7 (October 1990): 3747–49. https://doi.org/10.1063/1.346291.
        3) http://www.ioffe.ru/SVA/NSM/Semicond/InP
        4) Fiedler, F., and A. Schlachetzki. “Optical Parameters of InP-Based Waveguides.” Solid-State Electronics 30, no. 1 (January 1987): 73–83. https://doi.org/10.1016/0038-1101(87)90032-3.

        """

        if model is None:
            model = self.bandgap_model

        if P is None:
            P = self.P.to(self.reg.centimeter**-3)

        if T is None:
            T = self.T

        if y is None:
            y = self.y

        if x is None:
            x = self.x

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

            BGN = (
                A * P.magnitude ** (1 / 3)
                + B * P.magnitude**0.25
                + C * P.magnitude**0.5
            ) * self.reg.eV

        elif model == "none":
            BGN = 0

        # Eg formula taken from [5].
        if not get_BGN:
            Eg0 = 1.35 - 0.72 * y + 0.12 * y**2

            return Eg0 * self.reg.eV - BGN
        else:
            return BGN

    def get_fermi_level(self, P=None, T=None, bandgap_model=None):
        """
        Returns the fermi level calculated assuming Boltzman statistics from equation 21 from [1]. Assume Ec=0.
        It assumes that only uncompensated materials are given.
        It assumes full ionization.


        References:
        1) Sze, S. M., and Kwok Kwok Ng. Physics of Semiconductor Devices. 3rd ed. Hoboken, N.J: Wiley-Interscience, 2007.

        """

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        if P is None:
            P = self.P
        if T is None:
            T = self.T

        Eg = self.get_bandgap(P=P, T=T, model=bandgap_model)

        Ef = (
            -(np.log(P / self.Nv) + 1 / np.sqrt(8) * P / self.Nv) * self.kb * self.T
            - Eg
        )

        return Ef.to(self.reg.eV)

    def FD(self, E, Ef=1):
        return 1 / (1 + np.exp((E - Ef) / (self.kb * self.T)))

    def get_dalpha_BF(self, E=None, bandgap_model=None):
        """
        Gives the change in absorption due to the bandfilling effect.
        The treatment if based on [1] but we neglect the quasi-fermi levels and use the fermi level numerically calculated instead.

        Note that if the bandgap_model is not 'none' the BGN effect is taken into consideration within the bandfilling.

        E must be a Pint quantity in energy using the parent object's register!!

        References:
        1) Bennett, B.R., R.A. Soref, and J.A. Del Alamo. “Carrier-Induced Change in Refractive Index of InP, GaAs and InGaAsP.” IEEE Journal of Quantum Electronics 26, no. 1 (January 1990): 113–22. https://doi.org/10.1109/3.44924.

        """

        if E is None:
            E = self.energy

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        Ef = self.Ef
        Eg = self.get_bandgap(model=bandgap_model)

        Eah = (Eg - E) * (self.me / (self.me + self.mhh)) - Eg
        Eal = (Eg - E) * (self.me / (self.me + self.mhl)) - Eg
        Ebh = (E - Eg) * (self.mhh / (self.me + self.mhh))
        Ebl = (E - Eg) * (self.mhl / (self.me + self.mhl))

        alpha0 = self.Chh / E * np.sqrt(E - Eg) + self.Chl / E * np.sqrt(E - Eg)
        alpha = self.Chh / E * np.sqrt(E - Eg) * (
            self.FD(Eah, Ef=Ef) - self.FD(Ebh, Ef=Ef)
        ) + self.Chl / E * np.sqrt(E - Eg) * (self.FD(Eal, Ef=Ef) - self.FD(Ebl, Ef=Ef))

        alpha = np.nan_to_num(alpha)
        alpha0 = np.nan_to_num(alpha0)

        return (alpha - alpha0).to(self.reg.centimeter**-1)

    def get_alpha_sqrt(self, E=None, bandgap_model=None):
        if E is None:
            E = self.energy

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        Eg = self.get_bandgap(model=bandgap_model)

        alpha0 = self.Chh / E * np.sqrt(E - Eg) + self.Chl / E * np.sqrt(E - Eg)

        return alpha0.to(self.reg.centimeter**-1)

    def get_dalpha_iv(self, E=None, P=None):
        """
        Returns the intervalence absorption component. This is calculated from [1], eq. 16.

        This formula may give worse results at low dopings and underpredict absorption. [2]

        E: float. Pint quantity in energy using the parent object's register.

        References:
        1) Weber, J.-P. “Optimization of the Carrier-Induced Effective Index Change in InGaAsP Waveguides-Application to Tunable Bragg Filters.” IEEE Journal of Quantum Electronics 30, no. 8 (August 1994): 1801–16. https://doi.org/10.1109/3.301645.
        2) Casey, H. C., and P. L. Carter. “Variation of Intervalence Band Absorption with Hole Concentration in p -Type InP.” Applied Physics Letters 44, no. 1 (January 1, 1984): 82–83. https://doi.org/10.1063/1.94561.

        """

        if P is None:
            P = self.P

        if E is None:
            E = self.energy

        return (
            (
                4.252e-20
                * np.exp(E.to(self.reg.eV).magnitude * -3.657)
                * P.to(self.reg.meter**-3).magnitude
            )
            * self.reg.meter**-1
        ).to(self.reg.centimeter**-1)

    def get_dn_BF(self, E=None, bandgap_model=None, h=1e-3):
        """
        Returns the change in refractive index based only on the band filling effect based on the kramers kronig relations.

        E must be a Pint quantity in energy using the parent object's register!!

        h: float. energy step in eV.

        The algorithm is an adaptation of the Mclaurin algorithm described in [1].

        References:
        [1] Ohta, Koji, and Hatsuo Ishida. “Comparison among Several Numerical Integration Methods for Kramers-Kronig Transformation.” Applied Spectroscopy 42, no. 6 (August 1988): 952–57. https://doi.org/10.1366/0003702884430380.

        """

        if E is None:
            E = self.energy

        if type(E.magnitude) is not np.ndarray:
            E = np.asarray([E.magnitude]) * E.units

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        E = E.to(self.reg.eV)
        h = h * self.reg.eV

        Eg = self.get_bandgap(
            model="none"
        )  # this is left as 'none' because it is only used for the limits of integration. That way they are always the same.

        n_points_right = ((2 * Eg - E) / h).to(self.reg.dimensionless).astype(int)
        n_points_left = (
            ((E - 0.001 * self.reg.eV) / h).to(self.reg.dimensionless).astype(int)
        )
        Nright = np.max(n_points_right)
        Nleft = np.max(n_points_left)
        # print(Nright, Nleft)

        energies_right = E[..., None] + h * np.arange(Nright)[None, ...]
        energies_left = E[..., None] - h * (1 + np.arange(Nleft)[None, ...])

        energies = np.concatenate((energies_left, energies_right), axis=1)
        # print(energies.shape, Nleft, Nright, Nleft+Nright)
        if Nleft % 2 == 0:
            energies_integrand = energies[:, 1::2]
        else:
            energies_integrand = energies[:, 0::2]

        alpha = self.get_dalpha_BF(energies_integrand, bandgap_model=bandgap_model)
        # print(np.isclose((energies_integrand**2-E[..., None]**2),0).any())
        # print(np.isclose((energies_integrand**2-E[..., None]**2),0))
        # print((energies_integrand**2-E[..., None]**2))
        # integral=np.sum(self.get_dalpha_BF(energies_integrand)/(energies_integrand**2-E[..., None]**2), axis=1)*2*h
        # print(energies_integrand)
        integral = (
            np.sum(
                alpha
                * (
                    1 / (energies_integrand - E[..., None])
                    + 1 / (energies_integrand + E[..., None])
                )
                * 1
                / (2 * energies_integrand),
                axis=1,
            )
            * 2
            * h
        )
        # print(integral)
        return np.squeeze((integral * self.c * self.hbar / np.pi)).to(
            self.reg.dimensionless
        )  # the reason why you dont divide by e is because of the eV dimension on the integral result!!

    def get_dn_iv(self, E=None, P=None, bandgap_model="jain"):
        """
        Returns the intervalence absorption component. This is calculated from [1], eq. 17 and 18.

        E must be a Pint quantity in energy using the parent object's register!!

        References:
        1) Weber, J.-P. “Optimization of the Carrier-Induced Effective Index Change in InGaAsP Waveguides-Application to Tunable Bragg Filters.” IEEE Journal of Quantum Electronics 30, no. 8 (August 1994): 1801–16. https://doi.org/10.1109/3.301645.

        """

        if P is None:
            P = self.P

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
            * P
        ).to(self.reg.dimensionless)

    def get_dperm_pockels(self, E, Efield, bandgap_model=None):
        """
        Returns the change in refractive index owing to the pockels effect.

        E: energy of the field's photons. Pint quantity
        Efield: Electric field components. It is assumed to have dimensions (3,Ny, Nx)

        References:
        1) Adachi, Sadao, and Kunishige Oe. “Internal Strain and Photoelastic Effects in Ga  1− x  Al  x  As/GaAs and In  1− x  Ga  x  As  y  P  1− y  /InP Crystals.” Journal of Applied Physics 54, no. 11 (November 1983): 6620–27. https://doi.org/10.1063/1.331898.
        2) Adachi, Sadao, and Kunishige Oe. “Linear Electro‐optic Effects in Zincblende‐type Semiconductors: Key Properties of InGaAsP Relevant to Device Design.” Journal of Applied Physics 56, no. 1 (July 1984): 74–80. https://doi.org/10.1063/1.333731.

        """

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        def g(chi):
            return 1 / chi**2 * (2 - (1 + chi) ** -0.5 - (1 - chi) ** -0.5)

        def f(chi):
            return 1 / chi**2 * (2 - (1 + chi) ** 0.5 - (1 - chi) ** 0.5)

        def symmetric(x):
            """
            Returns an array in the standard symmetric form instead of the voigt notation. Assumes shape of [6, Ny, Nx]

            """
            y = np.zeros((3, 3) + x.shape[-2:]) * x.units

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

        Eg = self.get_bandgap()

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

        pockels_tensor = np.zeros((6, 3)) * r41_free / r41_free.magnitude

        pockels_tensor[3, 0] = r41_free + r41_piezo
        pockels_tensor[4, 1] = r41_free + r41_piezo
        pockels_tensor[5, 2] = r41_free + r41_piezo

        # Evaluate the change in impermeability
        deta = np.einsum("il,ljk->ijk", pockels_tensor, Efield)
        deta = symmetric(deta)  # return to symmetric matrix

        # build the permitivity tensor
        perm = np.zeros((3, 3)) * self.e0.units
        perm[0, 0] = self.eps_s * self.e0  # type: ignore
        perm[1, 1] = self.eps_s * self.e0
        perm[2, 2] = self.eps_s * self.e0

        # find the change in permitivity
        dperm = np.einsum("ik,kltp,lj->ijtp", perm, deta, perm) * -1 / self.e0

        return dperm

    def get_pockels_coeffs(self, E, bandgap_model=None):
        """
        Returns the pockels coefficients

        E: energy of the field's photons. Pint quantity
        Efield: Electric field components. It is assumed to have dimensions (3,Ny, Nx)

        References:
        1) Adachi, Sadao, and Kunishige Oe. “Internal Strain and Photoelastic Effects in Ga  1− x  Al  x  As/GaAs and In  1− x  Ga  x  As  y  P  1− y  /InP Crystals.” Journal of Applied Physics 54, no. 11 (November 1983): 6620–27. https://doi.org/10.1063/1.331898.
        2) Adachi, Sadao, and Kunishige Oe. “Linear Electro‐optic Effects in Zincblende‐type Semiconductors: Key Properties of InGaAsP Relevant to Device Design.” Journal of Applied Physics 56, no. 1 (July 1984): 74–80. https://doi.org/10.1063/1.333731.

        """

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        def g(chi):
            return 1 / chi**2 * (2 - (1 + chi) ** -0.5 - (1 - chi) ** -0.5)

        def f(chi):
            return 1 / chi**2 * (2 - (1 + chi) ** 0.5 - (1 - chi) ** 0.5)

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
        # print(E0, F0, C, D)
        Eg = self.get_bandgap(model=bandgap_model)

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
            ).to(self.reg.centimeter**2 / self.reg.dyne)
            * self.e14
        )

        return r41_free, r41_piezo

    def get_dperm_kerr(self, E, Efield, bandgap_model=None):
        """
        Returns the change in refractive index owing to the Kerr effect.

        E: energy of the field's photons. Pint quantity
        Efield: Electric field components. It is assumed to have dimensions (3,Ny, Nx)

        References:
        [1] - Maat, Derk Hendrik Pieter. InP-Based Integrated MZI Switches for Optical Communication, 2001.


        """

        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        def symmetric(x):
            """
            Returns an array in the standard symmetric form instead of the voigt notation. Assumes shape of [6, Ny, Nx]

            """
            y = np.zeros((3, 3) + x.shape[-2:]) * x.units

            y[0, 0] = x[0]
            y[0, 1] = y[1, 0] = x[5]
            y[0, 2] = y[2, 0] = x[4]
            y[1, 1] = x[1]
            y[1, 2] = y[2, 1] = x[3]
            y[2, 2] = x[2]

            return y

        # # Taken from [1]
        # A_TE = 0.25e3
        # A_TM = 0.20e3
        # B_TE = 0.71e9
        # B_TM = 0.48e9

        #These allow the replication of the results of [1] while using the Imodulator.
        A_TE = 0.9e3
        A_TM = 1.7e3
        B_TE = 0.42e9
        B_TM = 0.42e9

        # C_TE=-1.79e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2
        # C_TM=-1.82e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2

        C_TE = -3.10e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2
        C_TM = -5.60e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2

        A_mat = (
            np.asarray([[A_TE, 0, 0], [0, A_TM, 0], [0, 0, A_TE]])
            * self.reg.eV
            / (self.reg.volt * self.reg.meter)
        )

        B_mat = (
            np.asarray([[B_TE, 0, 0], [0, B_TM, 0], [0, 0, B_TE]])
            * self.reg.eV ** (-3 / 2)
            * self.reg.volt
            / self.reg.meter
        )

        freq = E / self.hbar / (2 * np.pi)
        wl = self.c / freq

        Eg = self.get_bandgap(model=bandgap_model)

        # print((2*np.pi*self.c/(Eg/self.hbar)).to(self.reg.nanometer))
        Efield_mag_sq = np.einsum("ijk,ijk -> jk", Efield.conjugate(), Efield).real

        dalpha = (
            A_mat[..., None, None]
            * wl
            * np.sqrt(Efield_mag_sq)[None, None, ...]
            / (Eg - E)
            * 10
            ** (
                -(
                    B_mat[..., None, None]
                    * (Eg - E) ** (3 / 2)
                    / np.sqrt(Efield_mag_sq)[None, None, ...]
                )
            )
        )
        dalpha = dalpha.to(self.reg.meter**-1)

        dperm_imag = (
            np.sqrt(self.eps_s) * self.c / (2 * np.pi * freq) * dalpha * self.e0
        )
        # print(type(E.magnitude))
        if type(E.magnitude) == np.ndarray:
            S11 = C_TE * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
            S12 = C_TM * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
            idx = np.where(E > Eg)
            S11[idx] = np.nan * self.reg.meter**2 * self.reg.volt**-2
            S12[idx] = np.nan * self.reg.meter**2 * self.reg.volt**-2

        else:
            if E > Eg:
                S11 = np.nan * self.reg.meter**2 * self.reg.volt**-2
                S12 = np.nan * self.reg.meter**2 * self.reg.volt**-2
            else:
                S11 = C_TE * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
                S12 = C_TM * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)

        S11 = S11.to(self.reg.meter**2 / self.reg.volt**2).magnitude
        S12 = S12.to(self.reg.meter**2 / self.reg.volt**2).magnitude
        S_mat = (
            np.asarray(
                [
                    [S11, S12, S12, 0, 0, 0],
                    [S12, S11, S12, 0, 0, 0],
                    [S12, S12, S11, 0, 0, 0],
                    [0, 0, 0, S11 - S12, 0, 0],
                    [0, 0, 0, 0, S11 - S12, 0],
                    [0, 0, 0, 0, 0, S11 - S12],
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
        deta_real = symmetric(np.einsum("ij,jkl->ikl", S_mat, Efield_voigt))

        # build the permitivity tensor
        perm = np.zeros((3, 3)) * self.e0.units
        perm[0, 0] = self.eps_s * self.e0
        perm[1, 1] = self.eps_s * self.e0
        perm[2, 2] = self.eps_s * self.e0

        # find the change in permitivity
        dperm = (
            np.einsum("ik,kltp,lj->ijtp", perm, deta_real, perm) * -1 / self.e0
            + 1j * dperm_imag
        )

        return dperm

    def get_kerr_coeffs(self, E, bandgap_model=None):
        if bandgap_model is None:
            bandgap_model = self.bandgap_model

        # C_TE=-1.79e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2
        # C_TM=-1.82e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2

        C_TE = -3.10e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2
        C_TM = -5.60e-18 * self.reg.eV**2 * self.reg.meter**2 / self.reg.volt**2

        Eg = self.get_bandgap(model=bandgap_model)

        if type(E.magnitude) == np.ndarray:
            S11 = C_TE * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
            S12 = C_TM * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
            idx = np.where(E > Eg)
            S11[idx] = np.nan * self.reg.meter**2 * self.reg.volt**-2
            S12[idx] = np.nan * self.reg.meter**2 * self.reg.volt**-2
        else:
            if E > Eg:
                S11 = np.nan * self.reg.meter**2 * self.reg.volt**-2
                S12 = np.nan * self.reg.meter**2 * self.reg.volt**-2
            else:
                S11 = C_TE * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)
                S12 = C_TM * E**2 / (np.sqrt(self.eps_s) ** 4 * (Eg**2 - E**2) ** 2)

        S11 = S11.to(self.reg.meter**2 / self.reg.volt**2).magnitude
        S12 = S12.to(self.reg.meter**2 / self.reg.volt**2).magnitude

        return S11, S12
