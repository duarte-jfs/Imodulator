import numpy as np
from matplotlib import pyplot as plt
import enlighten
import shapely

from scipy.interpolate import interp1d

from utils import Modulator_CPS


from imodulator.ChargeSimulator import ChargeSimulatorSolcore
from imodulator.OpticalSimulator import OpticalSimulatorFEMWELL
from imodulator.ElectroOpticalSimulator import ElectroOpticalSimulator
from imodulator.RFSimulator import RFSimulatorFEMWELL as RFSimulator

from scipy.special import erfc

from utils_data_analysis import get_EO
import skrf as rf

import os
import json

def truncate(x, decimals=2):
    factor = 10 ** decimals
    return np.trunc(x * factor) / factor


def count_zero_crossings(profile, threshold_ratio=0.05):
    """
    Count robust zero crossings in 1D profile.
    Ignores tiny values near zero.
    """
    profile = np.real(profile).copy()

    threshold = threshold_ratio * np.max(np.abs(profile))
    profile[np.abs(profile) < threshold] = 0

    # remove zeros
    profile = profile[profile != 0]

    if len(profile) < 2:
        return 0

    signs = np.sign(profile)
    return np.sum(signs[:-1] != signs[1:])

from skfem import ElementVector, ElementTriP1

def select_fundamental_TE_mode(modes, mode, modulator):
    te_frac_all = []
    pwr_wg2_all = []
    n_crossings_h_all = []
    n_crossings_v_all = []

    if isinstance(modulator, Modulator_CPS):
        x_ref = -modulator.w_metal_sig/2 - modulator.metal_sep/2
    y_ref = modulator.h_L2+modulator.h_L3/2


    for mode_j in modes:

        elements = mode.mesh.subdomains['L3']
        pwr_wg2 = np.abs(mode_j.calculate_confinement_factor(elements))**2
        te_frac = mode_j.te_fraction

        basis_j = mode_j.basis
        mesh_j = basis_j.mesh
        E_j = mode_j.E
        H_j = mode_j.H

        new_basis = basis_j.with_element(ElementVector(ElementTriP1()))

        (et, et_basis), (ez, ez_basis) = basis_j.split(E_j)
        E_j = new_basis.project(et_basis.interpolate(et), dtype=np.complex64)
        (et_x, et_x_basis), (et_y, et_y_basis) = new_basis.split(E_j)

        interp_et_x = et_x_basis.interpolator(et_x)

        x_coords_window = []
        y_coords_window = []
        for coord in mode.simulation_window.exterior.coords:
            x_coords_window.append(coord[0])
            y_coords_window.append(coord[1])

        # xmin, xmax = min(x_coords_window), max(x_coords_window)
        # ymin, ymax = min(y_coords_window), max(y_coords_window)
        
        # xmin += 0.01
        # xmax -= 0.01
        # ymin += 0.01
        # ymax -= 0.01

        xmin = -modulator.metal_sep/2 - modulator.w_metal_sig/2 - modulator.w_wg/2
        xmax = -modulator.metal_sep/2 - modulator.w_metal_sig/2 + modulator.w_wg/2
        ymin = -modulator.h_L1
        ymax = modulator.h_L2 + modulator.h_L3 + modulator.h_L4

        N=100
        line_v = np.asarray([np.ones(N)*x_ref, np.linspace(ymin, ymax, N)])
        line_h = np.asarray([np.linspace(xmin, xmax, N), np.ones(N)*y_ref])


        Efield_x_line_v = interp_et_x(line_v)
        Efield_x_line_h = interp_et_x(line_h)

        ## Now that we have the horizontal and the vertical cuts of the x component of the electric field we can calculate the amount of crossings the Re(Ex) has on the vertical and horizontal directions. 

        zero_crossings_v = count_zero_crossings(Efield_x_line_v)
        zero_crossings_h = count_zero_crossings(Efield_x_line_h)

        te_frac_all.append(te_frac)
        pwr_wg2_all.append(pwr_wg2)
        n_crossings_h_all.append(zero_crossings_h)
        n_crossings_v_all.append(zero_crossings_v)

    total_crossings_all = np.array(n_crossings_h_all) + np.array(n_crossings_v_all)

    #Now we select those that have 0 crossings all together
    fundamental_candidates = np.where(total_crossings_all == 0)[0]

    #Of these, we select the one with the highest confinement factor in the waveguide
    if len(fundamental_candidates) == 0:
        return np.argmax(np.asarray(te_frac_all)>0.5) #If there are no candidates with 0 crossings, we select the one with the highest TE fraction (which should be higher than 0.5 for a TE mode)
    
    best = fundamental_candidates[0]

    for candidate in fundamental_candidates:
        if pwr_wg2_all[candidate] > pwr_wg2_all[best]:
            best = candidate

    return best

def get_fom(
        modulator, 
        wavelength = 1.55,
        verbose = True,
        ):
    
    ## Retrieving optical modes

    data_to_return = {
        'ng': None,
        'n_opt': None,
        'loss_opt': None,
        'phi_V': None,
        'dalpha_V': None,
        'VpiL': None,
        'ER': None,
        '3dB_ee_bandwidth': None,
        '3dB_eo_bandwidth': None,
        'neff_RF': None,
        'loss_RF': None,
        'Z0_RF': None,
    }

    if verbose:
        print("Calculating TE mode...")

    sim_window = shapely.box(-10, -3, -5,3)

    mode = OpticalSimulatorFEMWELL(
        device=modulator.device,
        simulation_window=sim_window,
        include_metals=True
    )

    if verbose:
        print("Meshing...")

    # try:
    mode.make_mesh(verbose=True)
    # except Exception as e:
    #     print("ERROR: Error during meshing for optical mode:", e)
    #     return 1000, # Return a high FOM in case of meshing failure

    if verbose:
        print("Total mesh elements:", mode.mesh.nelements)

    if verbose:
        print("Computing modes...")

    neff_values = []
    wavelengths = [1.54,1.55,1.56]

    for wl in wavelengths:
        if verbose:
            print(f"Computing modes at wavelength {wl} um...")
        m = mode.compute_modes(
            wavelength=wl,
            num_modes=6,
            order=1,
            metallic_boundaries=True,
            n_guess=3.5,
            return_modes=True
        )

        fundamental_TE_mode_index = select_fundamental_TE_mode(m, mode, modulator)

        if verbose:
            print(f"Fundamental TE mode selected: {fundamental_TE_mode_index}")
        neff_values.append(m[fundamental_TE_mode_index].n_eff.real)

    # print("neff values: ", neff_values)
    # print("wavelengths: ", wavelengths)
    data_to_return['ng'] = neff_values[1] - wavelengths[1] * np.gradient(np.real(neff_values), wavelengths)[1]

    modes = mode.compute_modes(
                wavelength=1.55,
                num_modes=6,
                order=1,
                metallic_boundaries=True,
                n_guess=3.5,
                return_modes=True
            )

    if verbose:
        for i in range(len(modes)):
            print(f'Effective index - mode {i}: {modes[i].n_eff.real:.4f}')
            print(f'loss - mode {i}: {modes[i].calculate_propagation_loss(1e4):.2f} dB/cm')
            print(f"power: ", modes[i].calculate_power())
            print(f"TE fraction: ", modes[i].te_fraction)
            print(f"TM fraction: ", modes[i].tm_fraction)
            print(f"transversality: {modes[i].transversality}")
            print('#------------------------------------------#')


    fundamental_TE_mode_index = select_fundamental_TE_mode(modes, mode, modulator)

    TE_mode = modes[fundamental_TE_mode_index]

    if verbose:

        print(f"Fundamental TE mode selected: {fundamental_TE_mode_index}")


    if verbose:
        print("Transfering mode to device")

    mode.transfer_results_to_device(
        TE_TM_idx=[fundamental_TE_mode_index,0]
    )

    
    data_to_return['n_opt'] = TE_mode.n_eff.real
    data_to_return['loss_opt'] = TE_mode.calculate_propagation_loss(1e4) #dB/cm

    if verbose:
        print("Solving PDD...")
    
    a=shapely.LineString([
    [modulator.metal_sep/2+modulator.w_metal_sig/2,-modulator.h_L1+0.001],
    [modulator.metal_sep/2+modulator.w_metal_sig/2,modulator.h_L2+modulator.h_L3+modulator.h_L4+modulator.h_L5+modulator.h_L6-0.001]])

    charge=ChargeSimulatorSolcore(
        device=modulator.device,
        simulation_line=a,
        bias_start_stop_step=[0,6,21]
    )

    if verbose:
        print("Calculating charge distribution...")
    
    charge.solve_PDD(
        verbose = True*0,
        tol = 1e-6,
        max_iter = 1000,
        smooth_output=True,
        savgol_window = 10
    )
            
    if np.any(np.isnan(np.sum(charge.junction.pdd_output['Ec'], axis = 1))):
        print("ERROR: PDD solution contains NaN values. This may indicate a failure in the PDD solver.")

            
    if verbose:
        print("Transfer PDD results to device...")

    charge.transfer_results_to_device(
	dx = 0.1,
	polygons_to_map = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L2_right', 'L3_right', 'L4_right', 'L5_right', 'L6_right'],)


    if verbose:
        print("Calculating electro-optic response...")

    EO=ElectroOpticalSimulator(
            device=modulator.device,
            simulation_window=shapely.box(
            -modulator.w_metal_sig/2-modulator.metal_sep/2-2,
            -modulator.h_L1-1,
            -modulator.w_metal_sig/2-modulator.metal_sep/2+2,
            modulator.h_L2+modulator.h_L3+modulator.h_L4+modulator.h_L5+modulator.h_L6+1
        ),
        wavelength=1550
    )

    EO.make_mesh()

    if verbose:
        print("Creating delta epsilon")

    EO.get_epsilon_optical()

    all_res_TE = []

    for i in range(len(modulator.device.charge['V'])):
        kappa_aa = EO.calculate_EO_response(
            voltage_idx=i,
            rot_x=0,
            rot_y=np.pi/4,
            rot_z=0,
            base_epsilon_voltage_idx=0,
            optical_mode_a='TE',
            optical_mode_b='TE',
            wavelength=1.55
        )

        all_res_TE.append(kappa_aa['InGaAsPElectroOpticalModel']['results'])

    all_res_TE = np.asarray(all_res_TE)
    ## Interpolate the response to get better voltage resolution
    new_voltage_axis = np.linspace(modulator.device.charge['V'][0], np.sign(modulator.device.charge['V'][-1]) * 50, 100_000)
    interp_TE = interp1d(modulator.device.charge['V'], all_res_TE, axis=0, bounds_error=False, fill_value='extrapolate')
    all_res_TE_interp = interp_TE(new_voltage_axis)

    device_length = 1e-3
    MZM_output = 1/2 * (1 + np.exp(1j*np.sum(all_res_TE_interp*device_length, axis=1)))
    MZM_Iout = np.abs(MZM_output)**2

    ## Extract Vpi from EO response
    phase_change = np.sum(all_res_TE_interp*device_length, axis=1)
    idx_Vpi = np.argmin(np.abs(phase_change - np.pi))
    
    VpiL = new_voltage_axis[idx_Vpi]*device_length/1e-3 #Vmm
    ER = 10*np.log10(np.abs(1/2*(1+np.exp(1j*phase_change[idx_Vpi])))**2)

    data_to_return['phi_V'] = phase_change
    data_to_return['dalpha_V'] = np.imag(all_res_TE_interp*device_length)
    data_to_return['VpiL'] = VpiL
    data_to_return['ER'] = ER

    ## Run the RF simulation
    if verbose:
        print("Calculating RF response...")

    freq_refinement = 40 #GHz, we want to optimize for the EO response at this frequency
    N_mesh_refinement = 20 # Maximum number of mesh refinement steps
    N_elements_max = 10_000 #Maximum number of mesh elements. Refinement stops when either N_mesh refinement steps are reached or the number of mesh elements exceeds N_elements_max
    n_modes = 1 #We only want the fundamental mode for the RF simulation

    freq_values = np.linspace(1, 100, 12) #Values to compute the RF response 


    neff_values = np.zeros((n_modes, freq_values.shape[0]))
    loss_values = np.zeros(neff_values.shape)
    Z0_values = np.zeros(loss_values.shape, dtype = np.complex128)
    s11_values = np.zeros(loss_values.shape, dtype = np.complex128)
    s12_values = np.zeros(loss_values.shape, dtype = np.complex128)
    eo_values = np.zeros(loss_values.shape)

        
    rf_sim = RFSimulator(modulator.device, simulation_window=shapely.box(-75, -50, 0, 50))
    if verbose:
        print("meshing for RF simulation...")

    rf_sim.make_mesh(
        gmsh_algorithm=1,
        verbose=True
    )
       
    ## Start mesh refinement loop
    # Vpi_half = VpiL/2 #quadrature point
    # idx_quadrature = np.argmin(np.abs(modulator.device.charge['V'] - Vpi_half))

    boundary_conditions = ['left', 'bottom','top','right']
    for i in range(N_mesh_refinement):
                        
        #Compute modes
        modes2 = rf_sim.compute_modes(
            frequency = freq_refinement,
            metallic_boundaries = boundary_conditions,
            n_guess = 3,
            num_modes = 1,
            order = 1,
            return_modes = True,
            use_charge_transport_data=True,
            voltage_idx=int(len(modulator.device.charge['V'])/2)
        )

        idx_mode=0

        ######### Calculate characteristic impedance #########
        p0, i0_all, _ = rf_sim.get_currents(mode = rf_sim.modes[idx_mode])
        

        I0 = i0_all['sig_metal_line_current'] #<--- This is crucial!
        p0 = 2*p0 

        Z0 = p0/np.abs(I0)**2
        ######################################################

        rf_sim.refine_mesh(
            mode_for_refinement=rf_sim.modes[idx_mode], 
            N_nearest_neighbours=500
        )

        if rf_sim.mesh.nelements > N_elements_max:
            if verbose:
                print("Reached maximum number of elements. Stopping refinement.")
            break
            
        if verbose:
            a=rf_sim.mesh.nelements
            b=rf_sim.modes[idx_mode].n_eff.real
            c=rf_sim.modes[idx_mode].calculate_propagation_loss(1e4)
            d=np.abs(Z0)
            print(f'mode for refinement: {idx_mode}, neff: {b:.6f}, loss: {c:.2f} dB/cm, nelements: {a} | Z0: {d:.2f} Ohm')

    #Compute modes
    for idx_freq, freq in enumerate(freq_values):
        if verbose:
            print(f"Computing RF response at {freq:.2f} GHz...")
        modes2 = rf_sim.compute_modes(
            frequency = freq,
            metallic_boundaries = boundary_conditions,
            n_guess = 3,
            num_modes = n_modes,
            order = 1,
            return_modes = True,
            use_charge_transport_data=True,
            voltage_idx=int(len(modulator.device.charge['V'])/2)
        )


        for idx_mode in range(n_modes):
            neff_values[idx_mode, idx_freq] = modes2[idx_mode].n_eff.real
            loss_values[idx_mode, idx_freq] = modes2[idx_mode].calculate_propagation_loss(1e4) #dB/cm
            p0, i0_all, _ = rf_sim.get_currents(mode = modes2[idx_mode])

            I0 = i0_all['sig_metal_line_current'] #<--- This is crucial! 

            p0 = 2*p0 
            Z0_values[idx_mode, idx_freq] = p0/np.abs(I0)**2

            ########## Calculate S-parameters #########
            s11, s12 = rf_sim.get_S(
                rf_sim.modes[idx_mode].k,
                Z = Z0_values[idx_mode, idx_freq],
                ZL = 50,
                ZS = 50,
                L = 1e3 #um
            )

            s11_values[idx_mode, idx_freq] = s11
            s12_values[idx_mode, idx_freq] = s12
    
    eo = get_EO(
        freq_values,
        neff_values[0],
        loss_values[0],
        Z0_values[0],
        1e-3,
        N_new_freq=len(freq_values),
        no=data_to_return['ng']
    )

    eo_values[idx_mode] = eo

    freq_new_axis = np.linspace(freq_values[0], freq_values[-1], 10_000)
    interp_eo = interp1d(freq_values, eo_values[idx_mode], axis=0, bounds_error=False, fill_value='extrapolate')
    eo_values_interp = interp_eo(freq_new_axis)

    if any(eo_values_interp > -3):
        eo_bandwidth_idx = np.argmin(np.abs(eo_values_interp - -3))
    else:
        eo_bandwidth_idx = 0

    EO_bandwidth = freq_new_axis[eo_bandwidth_idx]

    interp_s12 = interp1d(freq_values, 20*np.log10(np.abs(s12_values[idx_mode])), axis=0, bounds_error=False, fill_value='extrapolate')
    s12_values_interp = interp_s12(freq_new_axis)

    if any(s12_values_interp > -3):
        s12_bandwidth_idx = np.argmin(np.abs(s12_values_interp - -3))
    else:
        s12_bandwidth_idx = 0

    EE_bandwidth = freq_new_axis[s12_bandwidth_idx]

    if verbose:
        print(f"EE bandwidth: {EE_bandwidth:.2f} GHz")

    data_to_return['3dB_ee_bandwidth'] = EE_bandwidth
    data_to_return['3dB_eo_bandwidth'] = EO_bandwidth
    data_to_return['neff_RF'] = neff_values[0]
    data_to_return['loss_RF'] = loss_values[0]
    data_to_return['Z0_RF'] = Z0_values[0]

    return data_to_return, rf_sim, mode, charge
    # return modulator, mode, None, None, None
