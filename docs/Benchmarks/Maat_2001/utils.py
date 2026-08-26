import numpy as np
import imodulator
import shapely
import openbandparams as obp
from imodulator.ElectroOpticalModel import InGaAsPElectroOpticalModel

import sys

from imodulator.utils import n_InGaAsP, p_InGaAsP, get_broadband_index

class Modulator_CPW:

    def __init__(
            self,
            w_window = 400, #um
            h_bottom = 500, #um
            h_top = 500, #um
            wavelength = 1550, #nm
            **kwargs
    ):
        self.wl = 1550
        self.w_wg = 4.5 #um
        self.w_metal_sig = 10 #um
        self.metal_sep = 10 #um
        self.w_metal_gnd = 10   #um
        self.h_metal = 0.4
        self.w_window = w_window
        self.h_bottom = h_bottom
        self.h_top = h_top

        self.overetch = 0.1
        self.base_width = 15

        self.h_L1 = 1
        self.h_L2 = 1.0
        self.h_L3 = 0.550
        self.h_etchstop = 0.02
        self.h_L4 = 0.19
        self.h_L5 = 1
        self.h_L6 = 0.05

        self.dop_L1 = 5e18
        self.dop_L2 = 7.2e16
        self.dop_L3 = 7.7e16
        self.dop_etchstop = 5e15
        self.dop_L4 = 5e15
        self.dop_L5 = -4.8e17
        self.dop_L6 = -1e19

        #As concentration
        self.y_L1 = 0
        self.y_L2 = 0
        self.y_L3 = 0.53
        self.y_L4 = 0
        self.y_etchstop = 0.53
        self.y_L5 = 0
        self.y_L6 = 1


        for kwarg in kwargs:
            setattr(self, kwarg, kwargs[kwarg])


        self.e = 1.60e-19 # electron charge in C
        self.e0 = 8.85e-12 # vacuum permittivity in F/m

        self.optical_mesh_settings = {
            'substrate': {'resolution': 0.5, 'SizeMax': 5, 'distance': 0.1},
            'background': {'resolution': 1, 'SizeMax': 5, 'distance': 0.1},
            'sig_metal': {'resolution': 0.8, 'SizeMax': 0.2, 'distance': 0.1},
            'n_metal_left': {'resolution': 0.1, 'SizeMax': 0.2, 'distance': 0.1},
            'n_metal_right': {'resolution': 0.1, 'SizeMax': 0.2, 'distance': 0.1},
            'bcb': {'resolution': 0.15, 'SizeMax': 0.2, 'distance': 0.1},
            'L6': {'resolution': 0.05, 'SizeMax': 5, 'distance': 0.1},
            'L5': {'resolution': 0.1, 'SizeMax': 5, 'distance': 0.1},
            'L4': {'resolution': 0.08, 'SizeMax': 5, 'distance': 0.1},
            'L3': {'resolution': 0.07, 'SizeMax': 5, 'distance': 0.1},
            'L2': {'resolution': 0.15, 'SizeMax': 5, 'distance': 0.1},
            'L1': {'resolution': 0.2, 'SizeMax': 5, 'distance': 0.1},
        }


        self.rf_mesh_settings = {
            'substrate': {'resolution': 7, 'SizeMax': 5, 'distance': 0.1},
            'background': {'resolution': 7, 'SizeMax': 5, 'distance': 0.1},
            'sig_metal': {'resolution': 1, 'SizeMax': 0.5, 'distance': 0.1},
            'n_metal_left': {'resolution': 1, 'SizeMax': 0.5, 'distance': 0.1},
            'n_metal_right': {'resolution': 1, 'SizeMax': 0.5, 'distance': 0.1},
            'bcb': {'resolution': 1, 'SizeMax': 1, 'distance': 0.1},
            'L6': {'resolution': 0.2, 'SizeMax': 5, 'distance': 0.1},
            'L5': {'resolution': 0.2, 'SizeMax': 5, 'distance': 0.1},
            'L4': {'resolution': 0.2, 'SizeMax': 5, 'distance': 0.1},
            'L3': {'resolution': 0.2, 'SizeMax': 5, 'distance': 0.1},
            'L2': {'resolution': 1, 'SizeMax': 5, 'distance': 0.1},
            'L1': {'resolution': 7, 'SizeMax': 5, 'distance': 0.1},
        }


        self.eo_mesh_settings = {
            'substrate': {'resolution': 0.5, 'SizeMax': 5, 'distance': 0.1},
            'background': {'resolution': 0.5, 'SizeMax': 5, 'distance': 0.1},
            'sig_metal': {'resolution': 10, 'SizeMax': 0.2, 'distance': 0.1},
            'n_metal_left': {'resolution': 10, 'SizeMax': 0.2, 'distance': 0.1},
            'n_metal_right': {'resolution': 10, 'SizeMax': 0.2, 'distance': 0.1},
            'bcb': {'resolution': 10, 'SizeMax': 0.2, 'distance': 0.1},
            'L6': {'resolution': 0.1, 'SizeMax': 5, 'distance': 0.1},
            'L5': {'resolution': 0.1, 'SizeMax': 5, 'distance': 0.1},
            'L4': {'resolution': 0.1, 'SizeMax': 5, 'distance': 0.1},
            'L3': {'resolution': 0.05, 'SizeMax': 5, 'distance': 0.1},
            'L2': {'resolution': 0.6, 'SizeMax': 5, 'distance': 0.1},
            'L1': {'resolution': 0.6, 'SizeMax': 5, 'distance': 0.1}
        }

        self.charge_mesh_settings = {
            'substrate': {'resolution': 0.05},
            'background': {'resolution': 0.05},
            'sig_metal': {'resolution': 0.01},
            'n_metal_left': {'resolution': 0.01},
            'n_metal_right': {'resolution': 0.01},
            'bcb': {'resolution': 0.01},
            'L6': {'resolution': 0.001},
            'L5': {'resolution': 0.001},
            'L4': {'resolution': 0.002},
            'L3': {'resolution': 0.002},
            'L2': {'resolution': 0.002},
            'L1': {'resolution': 0.002},
        }



        self._create_polygons()

        photo_polygons = [
                self.sig_metal,
                self.gnd_metal_left,
                self.gnd_metal_right,
                self.L6,
                self.L5,
                self.etchstop,
                self.L4,
                self.L3,
                self.L2,
                self.L1,
                self.bcbp,
                self.bcbp_left,
                self.bcbp_right,
                self.substrate,
                self.background,
             ]
        
        
        #Just in case there are empty polygons
        idxs_to_remove = []
        for i, poly in enumerate(photo_polygons):
            if np.isclose(poly.polygon.bounds[1], poly.polygon.bounds[3]):
                idxs_to_remove.append(i)
        for i in idxs_to_remove[::-1]:
            del photo_polygons[i]
        self.device = imodulator.PhotonicDevice(
            photo_polygons
        )


        
    def _create_polygons(self):
        freq = np.logspace(np.log10(0.1), np.log10(200), 150) #GHz. This will be the simulation frequency
        
        eps_rf_metal = 1 - 1j*4.5e7/(2*np.pi*freq*1e9 * self.e0)
        eps_opt_gold = 11.259-1j*115.13

        bcb_eps_real = 3.2*np.ones(len(freq))
        bcb_eps_imag = bcb_eps_real * 0.001

        bcb_eps = bcb_eps_real - 1j*bcb_eps_imag

        optical_wavelengths = np.linspace(self.wl/1e3-0.05, self.wl/1e3+0.05, 10)

        substrate_obp = obp.GaInPAs(T=300, As = 0, a = obp.InP.a())
        self.substrate = imodulator.SemiconductorPolygon(
            shapely.box(
                -self.w_window/2,
                -self.h_bottom,
                self.w_window/2,
                0
            ),
            rf_eps=np.asarray([freq,substrate_obp.dielectric(T = 300) - 1j*57.6e3/(2*np.pi*freq*1e9 * self.e0)]),
            optical_material = substrate_obp.refractive_index(T=300)**2,
            eo_mesh_settings = self.eo_mesh_settings['substrate'],
            rf_mesh_settings = self.rf_mesh_settings['substrate'],
            optical_mesh_settings = self.optical_mesh_settings['substrate'],
            name = 'substrate',
            charge_transport_simulator_kwargs={
                'sol_obp_material': substrate_obp,
                'sol_Nd': 1e15,
            }
        )

        L1_obp = obp.GaInPAs(T=300, As = self.y_L1, a = obp.InP.a())
        if self.dop_L1 >= 0:
            nsq=get_broadband_index(
                n_InGaAsP(self.dop_L1, 300, self.y_L1, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': L1_obp,
                'sol_Nd': self.dop_L1,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L1_obp.element_fraction('Ga'),
                'alloy_y': L1_obp.element_fraction('As'),
                'doping_type': 'n',
                'doping_conc': self.dop_L1,
            }
        else:
            nsq=get_broadband_index(
                p_InGaAsP(-self.dop_L1, 300, self.y_L1, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': L1_obp,
                'sol_Na': -self.dop_L1,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L1_obp.element_fraction('Ga'),
                'alloy_y': L1_obp.element_fraction('As'),
                'doping_type': 'p',
                'doping_conc': -self.dop_L1,
             }

        self.L1 = imodulator.SemiconductorPolygon(
                    shapely.box(
                        -self.w_window/2,
                        0,
                        self.w_window/2,
                        self.h_L1)
                    ,
                    rf_eps=L1_obp.dielectric(T = 300),
                    optical_material = [optical_wavelengths, nsq],
                    eo_mesh_settings = self.eo_mesh_settings['L1'],
                    rf_mesh_settings = self.rf_mesh_settings['L1'],
                    optical_mesh_settings = self.optical_mesh_settings['L1'],
                    charge_mesh_settings=self.charge_mesh_settings['L1'],
                    name = f'L1',
                    electro_optic_module=InGaAsPElectroOpticalModel,
                    electro_optic_module_kwargs = {
                            'y': self.y_L1, 
                            'T': 300, 
                            'BF_model': 'vinchant',
                        },
                    charge_transport_simulator_kwargs=charge_transport_simulator_kwargs
                )
        
        L2_obp = obp.GaInPAs(T=300, As = self.y_L2, a = obp.InP.a())
        if self.dop_L2 >= 0:
            nsq=get_broadband_index(
                n_InGaAsP(self.dop_L2, 300, self.y_L2, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': L2_obp,
                'sol_Nd': self.dop_L2,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L2_obp.element_fraction('Ga'),
                'alloy_y': L2_obp.element_fraction('As'),
                'doping_type': 'n',
                'doping_conc': self.dop_L2,
            }
        else:
            nsq=get_broadband_index(
                p_InGaAsP(-self.dop_L2, 300, self.y_L2, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': L2_obp,
                'sol_Na': -self.dop_L2,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L2_obp.element_fraction('Ga'),
                'alloy_y': L2_obp.element_fraction('As'),
                'doping_type': 'p',
                'doping_conc': -self.dop_L2,
             }

        pts = [
            [-self.w_window/2, self.h_L1],
            [-self.w_window/2, self.h_L1 + self.h_L2 - self.overetch],
            [-self.base_width/2, self.h_L1 + self.h_L2 - self.overetch],
            [-self.base_width/2, self.h_L1 + self.h_L2],
            [self.base_width/2, self.h_L1 + self.h_L2],
            [self.base_width/2, self.h_L1 + self.h_L2 - self.overetch],
            [self.w_window/2, self.h_L1 + self.h_L2 - self.overetch],
            [self.w_window/2, self.h_L1],
        ]
        self.L2 = imodulator.SemiconductorPolygon(
                    shapely.Polygon(pts),
                    rf_eps=L2_obp.dielectric(T = 300),
                    optical_material = [optical_wavelengths, nsq],
                    eo_mesh_settings = self.eo_mesh_settings['L2'],
                    rf_mesh_settings = self.rf_mesh_settings['L2'],
                    optical_mesh_settings = self.optical_mesh_settings['L2'],
                    charge_mesh_settings=self.charge_mesh_settings['L2'],
                    name = f'L2',
                    electro_optic_module=InGaAsPElectroOpticalModel,
                    electro_optic_module_kwargs = {
                            'y': self.y_L2, 
                            'T': 300, 
                            'BF_model': 'vinchant',
                        },
                    charge_transport_simulator_kwargs=charge_transport_simulator_kwargs
                )
        

        L3_obp = obp.GaInPAs(T=300, As = self.y_L3, a = obp.InP.a())
        if self.dop_L3 >= 0:
            nsq=get_broadband_index(
                n_InGaAsP(self.dop_L3, 300, self.y_L3, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': L3_obp,
                'sol_Nd': self.dop_L3,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L3_obp.element_fraction('Ga'),
                'alloy_y': L3_obp.element_fraction('As'),
                'doping_type': 'n',
                'doping_conc': self.dop_L3,
                
            }
        else:
            nsq=get_broadband_index(
                p_InGaAsP(-self.dop_L3, 300, self.y_L3, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': L3_obp,
                'sol_Na': -self.dop_L3,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L3_obp.element_fraction('Ga'),
                'alloy_y': L3_obp.element_fraction('As'),
                'doping_type': 'p',
                'doping_conc': -self.dop_L3,
             }

        pts = [
                [-self.base_width/2, self.h_L1 + self.h_L2],
                [-self.base_width/2, self.h_L1 + self.h_L2 + self.h_L3 - self.overetch],
                [-self.w_wg/2, self.h_L1 + self.h_L2 + self.h_L3 - self.overetch],
                [-self.w_wg/2, self.h_L1 + self.h_L2 + self.h_L3],
                [self.w_wg/2, self.h_L1 + self.h_L2 + self.h_L3],
                [self.w_wg/2, self.h_L1 + self.h_L2 + self.h_L3 - self.overetch],
                [self.base_width/2, self.h_L1 + self.h_L2 + self.h_L3 - self.overetch],
                [self.base_width/2, self.h_L1 + self.h_L2],
            ]
        self.L3 = imodulator.SemiconductorPolygon(
                    shapely.Polygon(pts),
                    rf_eps=L3_obp.dielectric(T = 300),
                    optical_material = [optical_wavelengths, nsq],
                    eo_mesh_settings = self.eo_mesh_settings['L3'],
                    rf_mesh_settings = self.rf_mesh_settings['L3'],
                    optical_mesh_settings = self.optical_mesh_settings['L3'],
                    charge_mesh_settings=self.charge_mesh_settings['L3'],
                    name = f'L3',
                    electro_optic_module=InGaAsPElectroOpticalModel,
                    electro_optic_module_kwargs = {
                            'y': self.y_L3, 
                            'T': 300, 
                            'BF_model': 'vinchant',
                        },
                    charge_transport_simulator_kwargs=charge_transport_simulator_kwargs
                )
        
        L4_obp = obp.GaInPAs(T=300, As = self.y_L4, a = obp.InP.a())
        if self.dop_L4 >= 0:
            nsq=get_broadband_index(
                n_InGaAsP(self.dop_L4, 300, self.y_L4, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': L4_obp,
                'sol_Nd': self.dop_L4,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L4_obp.element_fraction('Ga'),
                'alloy_y': L4_obp.element_fraction('As'),
                'doping_type': 'n',
                'doping_conc': self.dop_L4,
            }
        else:
            nsq=get_broadband_index(
                p_InGaAsP(-self.dop_L4, 300, self.y_L4, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': L4_obp,
                'sol_Na': -self.dop_L4,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L4_obp.element_fraction('Ga'),
                'alloy_y': L4_obp.element_fraction('As'),
                'doping_type': 'p',
                'doping_conc': -self.dop_L4,
             }

        self.L4 = imodulator.SemiconductorPolygon(
                    shapely.box(
                        -self.w_wg/2,
                        self.h_L1 + self.h_L2 + self.h_L3,
                        self.w_wg/2,
                        self.h_L1 + self.h_L2 + self.h_L3 + self.h_L4)
                    ,
                    rf_eps=L4_obp.dielectric(T = 300),
                    optical_material = [optical_wavelengths, nsq],
                    eo_mesh_settings = self.eo_mesh_settings['L4'],
                    rf_mesh_settings = self.rf_mesh_settings['L4'],
                    optical_mesh_settings = self.optical_mesh_settings['L4'],
                    charge_mesh_settings=self.charge_mesh_settings['L4'],
                    name = f'L4',
                    electro_optic_module=InGaAsPElectroOpticalModel,
                    electro_optic_module_kwargs = {
                            'y': self.y_L4, 
                            'T': 300, 
                            'BF_model': 'vinchant',
                        },
                    charge_transport_simulator_kwargs=charge_transport_simulator_kwargs
                )

        etchstop_obp = obp.GaInPAs(T=300, As = self.y_etchstop, a = obp.InP.a())
        if self.dop_etchstop >= 0:
            nsq=get_broadband_index(
                n_InGaAsP(self.dop_etchstop, 300, self.y_etchstop, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': etchstop_obp,
                'sol_Nd': self.dop_etchstop,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': etchstop_obp.element_fraction('Ga'),
                'alloy_y': etchstop_obp.element_fraction('As'),
                'doping_type': 'n',
                'doping_conc': self.dop_etchstop,
            }
        else:
            nsq=get_broadband_index(
                p_InGaAsP(-self.dop_etchstop, 300, self.y_etchstop, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': etchstop_obp,
                'sol_Na': -self.dop_etchstop,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': etchstop_obp.element_fraction('Ga'),
                'alloy_y': etchstop_obp.element_fraction('As'),
                'doping_type': 'p',
                'doping_conc': -self.dop_etchstop,
            }

        self.etchstop = imodulator.SemiconductorPolygon(
                    shapely.box(
                        -self.w_wg/2,
                        self.h_L1 + self.h_L2 + self.h_L3 + self.h_L4,
                        self.w_wg/2,
                        self.h_L1 + self.h_L2 + self.h_L3 + self.h_L4 + self.h_etchstop)
                    ,
                    rf_eps=etchstop_obp.dielectric(T = 300),
                    optical_material = [optical_wavelengths, nsq],
                    eo_mesh_settings = self.eo_mesh_settings['L4'],
                    rf_mesh_settings = self.rf_mesh_settings['L4'],
                    optical_mesh_settings = self.optical_mesh_settings['L4'],
                    charge_mesh_settings=self.charge_mesh_settings['L4'],
                    name = f'etchstop',
                    electro_optic_module=InGaAsPElectroOpticalModel,
                    electro_optic_module_kwargs = {
                            'y': self.y_etchstop, 
                            'T': 300, 
                            'BF_model': 'vinchant',
                        },
                    charge_transport_simulator_kwargs=charge_transport_simulator_kwargs
                    )
        
        L5_obp = obp.GaInPAs(T=300, As = self.y_L5, a = obp.InP.a())
        if self.dop_L5 >= 0:
            nsq=get_broadband_index(
                n_InGaAsP(self.dop_L5, 300, self.y_L5, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': L5_obp,
                'sol_Nd': self.dop_L5,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L5_obp.element_fraction('Ga'),
                'alloy_y': L5_obp.element_fraction('As'),
                'doping_type': 'n',
                'doping_conc': self.dop_L5,
            }
        else:
            nsq=get_broadband_index(
                p_InGaAsP(-self.dop_L5, 300, self.y_L5, self.wl),
                wl_min=optical_wavelengths[0]*1e3, wl_max=optical_wavelengths[-1]*1e3, N_samples=len(optical_wavelengths))[:,1]
            charge_transport_simulator_kwargs={
                'sol_obp_material': L5_obp,
                'sol_Na': -self.dop_L5,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L5_obp.element_fraction('Ga'),
                'alloy_y': L5_obp.element_fraction('As'),
                'doping_type': 'p',
                'doping_conc': -self.dop_L5,
            }

        self.L5 = imodulator.SemiconductorPolygon(
                    shapely.box(
                        -self.w_wg/2,
                        self.h_L1 + self.h_L2 + self.h_L3 + self.h_L4 + self.h_etchstop,
                        self.w_wg/2,
                        self.h_L1 + self.h_L2 + self.h_L3 + self.h_L4 + self.h_etchstop + self.h_L5)
                    ,
                    rf_eps=L5_obp.dielectric(T = 300),
                    optical_material = [optical_wavelengths, nsq],
                    eo_mesh_settings = self.eo_mesh_settings['L5'],
                    rf_mesh_settings = self.rf_mesh_settings['L5'],
                    optical_mesh_settings = self.optical_mesh_settings['L5'],
                    charge_mesh_settings=self.charge_mesh_settings['L5'],
                    name = f'L5',
                    electro_optic_module=InGaAsPElectroOpticalModel,
                    electro_optic_module_kwargs = {
                            'y': self.y_L5, 
                            'T': 300, 
                            'BF_model': 'vinchant',
                        },
                    charge_transport_simulator_kwargs=charge_transport_simulator_kwargs
                )
        
        L6_obp = obp.GaInPAs(T=300, As = self.y_L6, a = obp.InP.a())
        if self.dop_L6 >= 0:
            charge_transport_simulator_kwargs={
                'sol_obp_material': L6_obp,
                'sol_Nd': self.dop_L6,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L6_obp.element_fraction('Ga'),
                'alloy_y': L6_obp.element_fraction('As'),
                'doping_type': 'n',
                'doping_conc': self.dop_L6,
            }
        else:
            charge_transport_simulator_kwargs={
                'sol_obp_material': L6_obp,
                'sol_Na': -self.dop_L6,
                'material_definition':"Ga(x)In(1-x)As(y)P(1-y)" ,
                'alloy_type': 'quaternary_constant',
                'alloy_x': L6_obp.element_fraction('Ga'),
                'alloy_y': L6_obp.element_fraction('As'),
                'doping_type': 'p',
                'doping_conc': -self.dop_L6,
             }

        self.L6 = imodulator.SemiconductorPolygon(
                    shapely.box(
                        -self.w_wg/2,
                        self.h_L1 + self.h_L2 + self.h_L3 + self.h_L4 + self.h_etchstop + self.h_L5,
                        self.w_wg/2,
                        self.h_L1 + self.h_L2 + self.h_L3 + self.h_L4 + self.h_etchstop + self.h_L5 +self.h_L6)
                    ,
                    rf_eps=L6_obp.dielectric(T = 300),
                    optical_material = (3.53-1j*0.076399)**2,
                    eo_mesh_settings = self.eo_mesh_settings['L6'],
                    rf_mesh_settings = self.rf_mesh_settings['L6'],
                    optical_mesh_settings = self.optical_mesh_settings['L6'],
                    charge_mesh_settings=self.charge_mesh_settings['L6'],
                    name = f'L6',
                    electro_optic_module=InGaAsPElectroOpticalModel,
                    electro_optic_module_kwargs = {
                            'y': self.y_L6, 
                            'T': 300, 
                            'BF_model': 'vinchant',
                        },
                    charge_transport_simulator_kwargs=charge_transport_simulator_kwargs
                )
        

        h_mesa = (self.h_L1 + self.h_L2 + self.h_L3 + self.h_L4 + self.h_L5 + self.h_L6 + self.h_etchstop)
        self.sig_metal = imodulator.MetalPolygon(
            shapely.box(
                -self.w_metal_sig/2,
                h_mesa,
                self.w_metal_sig/2,
                h_mesa + self.h_metal
            ),
            rf_eps=np.asarray([freq,eps_rf_metal]),
            optical_material = eps_opt_gold,
            eo_mesh_settings = self.eo_mesh_settings['sig_metal'],
            rf_mesh_settings = self.rf_mesh_settings['sig_metal'],
            optical_mesh_settings = self.optical_mesh_settings['sig_metal'],
            name = 'sig_metal',
            calculate_current=True,
            d_buffer_current = 0.005
        )

        self.gnd_metal_left = imodulator.MetalPolygon(
                    shapely.box(
                        -self.w_metal_sig/2 - self.metal_sep - self.w_metal_gnd,
                        self.h_L1 + self.h_L2 - self.overetch,
                        -self.w_metal_sig/2 - self.metal_sep,
                        self.h_L1 + self.h_L2 + self.h_metal - self.overetch
                    ),
                    rf_eps=np.asarray([freq,eps_rf_metal]),
                    optical_material = eps_opt_gold,
                    eo_mesh_settings = self.eo_mesh_settings['sig_metal'],
                    rf_mesh_settings = self.rf_mesh_settings['sig_metal'],
                    optical_mesh_settings = self.optical_mesh_settings['sig_metal'],
                    name = 'gnd_metal_left',
                    calculate_current=True,
                    d_buffer_current = np.min([self.h_metal/30, self.w_metal_sig/30])
                )

        self.gnd_metal_right = imodulator.MetalPolygon(
                    shapely.box(
                        self.w_metal_sig/2 + self.metal_sep,
                        self.h_L1 + self.h_L2 - self.overetch,
                        self.w_metal_sig/2 + self.metal_sep + self.w_metal_gnd,
                        self.h_L1 + self.h_L2 + self.h_metal - self.overetch
                    ),
                    rf_eps=np.asarray([freq,eps_rf_metal]),
                    optical_material = eps_opt_gold,
                    eo_mesh_settings = self.eo_mesh_settings['sig_metal'],
                    rf_mesh_settings = self.rf_mesh_settings['sig_metal'],
                    optical_mesh_settings = self.optical_mesh_settings['sig_metal'],
                    name = 'gnd_metal_right',
                    calculate_current=True,
                    d_buffer_current = np.min([self.h_metal/30, self.w_metal_sig/30])
                )

        self.bcbp_left = imodulator.InsulatorPolygon(
            shapely.box(
                -self.w_window/2,
                self.h_L1,
                -self.w_metal_sig/2 - self.metal_sep-self.w_metal_gnd,
                h_mesa
            ),
            rf_eps=np.asarray([freq, bcb_eps]),
            optical_material = 1.56**2,
            eo_mesh_settings = self.eo_mesh_settings['bcb'],
            rf_mesh_settings = self.rf_mesh_settings['bcb'],
            optical_mesh_settings = self.optical_mesh_settings['bcb'],
            name = 'bcbp_left',
        )

        poly = shapely.box(
                         - self.w_metal_sig/2 - self.metal_sep,
                        self.h_L1 + self.h_L2 - self.overetch,
                         self.w_metal_sig/2 + self.metal_sep,
                        h_mesa
                    )

        for poly_to_remove in [
            self.L1.polygon,
            self.L2.polygon,
            self.L3.polygon,
            self.L4.polygon,
            self.etchstop.polygon,
            self.L5.polygon,
            self.L6.polygon,
        ]:
            poly = poly.difference(poly_to_remove)

        self.bcbp = imodulator.InsulatorPolygon(
            poly,
            rf_eps=np.asarray([freq, bcb_eps]),
            optical_material = 1.56**2,
            eo_mesh_settings = self.eo_mesh_settings['bcb'],
            rf_mesh_settings = self.rf_mesh_settings['bcb'],
            optical_mesh_settings = self.optical_mesh_settings['bcb'],
            name = 'bcbp',
        )

        self.bcbp_right = imodulator.InsulatorPolygon(
            shapely.box(
                self.w_metal_sig/2 + self.metal_sep+self.w_metal_gnd,
                self.h_L1,
                self.w_window/2,
                h_mesa
            ),
            rf_eps=np.asarray([freq, bcb_eps]),
            optical_material = 1.56**2,
            eo_mesh_settings = self.eo_mesh_settings['bcb'],
            rf_mesh_settings = self.rf_mesh_settings['bcb'],
            optical_mesh_settings = self.optical_mesh_settings['bcb'],
            name = 'bcbp_right',
        )
        

        self.background = imodulator.InsulatorPolygon(
            shapely.box(
                -self.w_window/2,
                -self.h_bottom,
                self.w_window/2,
                (h_mesa+
                 self.h_metal +
                 self.h_top)
            ),

            rf_eps = 1,
            optical_material = 1,
            eo_mesh_settings = self.eo_mesh_settings['background'],
            rf_mesh_settings = self.rf_mesh_settings['background'],
            optical_mesh_settings = self.optical_mesh_settings['background'],
            name = 'background'
        )



if __name__ == "__main__":

    from matplotlib import pyplot as plt
    
    mod = Modulator_CPW()
    mod.device.plot_polygons(fill_polygons=True)
    plt.show()