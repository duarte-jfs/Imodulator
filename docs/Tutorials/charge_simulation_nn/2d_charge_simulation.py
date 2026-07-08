
"""
2D DC charge transport with nextnano++
======================================

This tutorial solves the charge transport of the InP-based InGaAsP modulator
cross-section on a *full 2D* domain using ``ChargeSimulatorNN2D``.

Unlike the ``dc_electro_optical`` tutorial - which solves a 1D Poisson-drift-
diffusion problem along a vertical line with ``ChargeSimulatorSolcore`` - here we
hand nextnano++ a 2D ``simulation_window`` polygon and let it resolve the true
(Ex, Ey) field and carrier distribution over the whole cross-section, including
the lateral current path from the top p-contact to the side n-contact.

Requirements
------------
* nextnano++ installed and configured (see ``imodulator.Config``).
* Run with pre-releases allowed::

    uv run --prerelease=allow python docs/Tutorials/2D_dc_electro_optical/2D_dc_electro_optical.py
"""

# Workaround for a matplotlib 3.9.0 Tk bug on Windows:
# Win32_GetForegroundWindow() raises "PyCapsule_New called with null pointer"
# when no window has focus. Fixed upstream in matplotlib 3.9.1. Must run before
# matplotlib.pyplot is imported.
import sys
if sys.platform == "win32":
    try:
        from matplotlib import _c_internal_utils
        _orig_get_fg = _c_internal_utils.Win32_GetForegroundWindow
        def _safe_get_fg():
            try:
                return _orig_get_fg()
            except ValueError:
                return None
        _c_internal_utils.Win32_GetForegroundWindow = _safe_get_fg
    except Exception:
        pass

import numpy as np
import shapely
import matplotlib.pyplot as plt

import imodulator
from imodulator.ChargeSimulator import ChargeSimulatorNN2D


# ---------------------------------------------------------------------------
# Device geometry (all dimensions in um)
# ---------------------------------------------------------------------------
# Vertical p-i-n stack, from bottom to top: n / wg1 / wg2 / p1 / p2.
# The signal metal sits on top (biased p-contact); a ground metal sits on the
# n-layer to the left, so current flows down the stack and laterally through n.

w_wg = 1.0            # stack width
h_n = 0.4
h_wg1 = 0.5
h_wg2 = 0.3
h_p1 = 1.0
h_p2 = 0.2

w_sig_metal = 5.0     # signal (p) metal width
metal_sep = 10.0      # gap between signal and ground metals
w_gnd_metal = 10.0    # ground (n) metal width
h_metal = 4.0

# cumulative y of each stack interface
y_n0 = 0.0
y_n1 = y_n0 + h_n
y_wg1 = y_n1 + h_wg1
y_wg2 = y_wg1 + h_wg2
y_p1 = y_wg2 + h_p1
y_p2 = y_p1 + h_p2      # top of the semiconductor stack

# n-layer spans the full width so the side contact can reach it
x_n_edge = w_sig_metal / 2 + metal_sep + w_gnd_metal


# ---------------------------------------------------------------------------
# nextnano charge kwargs per region
# ---------------------------------------------------------------------------
# material_definition defaults to "Ga(x)In(1-x)As(y)P(1-y)"; alloy_x is the
# InP-lattice-matched Ga fraction for the given As fraction alloy_y.
#   alloy_y = 0.00 -> InP            (alloy_x = 0.000)
#   alloy_y = 0.53 -> Q1.3 InGaAsP   (alloy_x = 0.245)
#   alloy_y = 1.00 -> Ga0.47In0.53As (alloy_x = 0.468)

def charge_kwargs(alloy_x, alloy_y, doping_type, doping_conc):
    return {
        "material_definition": None,
        "alloy_x": alloy_x,
        "alloy_y": alloy_y,
        "doping_type": doping_type,
        "doping_conc": doping_conc,
    }


n = imodulator.SemiconductorPolygon(
    shapely.box(-x_n_edge, y_n0, x_n_edge, y_n1),
    name="n",
    charge_mesh_settings={"resolution": 0.003},
    charge_transport_simulator_kwargs=charge_kwargs(0.0, 0.0, "n", 1e18),
)

wg1 = imodulator.SemiconductorPolygon(
    shapely.box(-w_wg / 2, y_n1, w_wg / 2, y_wg1),
    name="wg1",
    charge_mesh_settings={"resolution": 0.002},
    charge_transport_simulator_kwargs=charge_kwargs(0.245, 0.53, "n", 1e16),
)

wg2 = imodulator.SemiconductorPolygon(
    shapely.box(-w_wg / 2, y_wg1, w_wg / 2, y_wg2),
    name="wg2",
    charge_mesh_settings={"resolution": 0.002},
    charge_transport_simulator_kwargs=charge_kwargs(0.0, 0.0, "n", 1e16),
)

p1 = imodulator.SemiconductorPolygon(
    shapely.box(-w_wg / 2, y_wg2, w_wg / 2, y_p1),
    name="p1",
    charge_mesh_settings={"resolution": 0.003},
    charge_transport_simulator_kwargs=charge_kwargs(0.0, 0.0, "p", 1e17),
)

p2 = imodulator.SemiconductorPolygon(
    shapely.box(-w_wg / 2, y_p1, w_wg / 2, y_p2),
    name="p2",
    charge_mesh_settings={"resolution": 0.003},
    charge_transport_simulator_kwargs=charge_kwargs(0.468, 1.0, "p", 1e19),
)

# First MetalPolygon in the device becomes contact1 (biased), second contact2
# (grounded). We list the signal (p) metal first, ground (n) metal second.
sig_metal = imodulator.MetalPolygon(
    shapely.box(-w_sig_metal / 2, y_p2, w_sig_metal / 2, y_p2 + h_metal),
    name="sig_metal",
)

n_metal_left = imodulator.MetalPolygon(
    shapely.box(-x_n_edge, y_n1, -(w_sig_metal / 2 + metal_sep), y_n1 + h_metal),
    name="n_metal_left",
)

device = imodulator.PhotonicDevice(
    [sig_metal, n_metal_left, p2, p1, wg2, wg1, n]
)


# ---------------------------------------------------------------------------
# 2D charge transport simulation
# ---------------------------------------------------------------------------
# The window captures the vertical junction, the signal metal, and the side
# ground contact - but stops short of the (absent) right metal.
window = shapely.box(-(x_n_edge + 0.5), -0.05, w_sig_metal / 2 + 0.5, y_p2 + h_metal + 0.1)

charge2d = ChargeSimulatorNN2D(
    device=device,
    simulation_window=window,
    bias_start_stop_step=[0, 5, 11],
)

# Inspect geometry and detected contacts before running.
charge2d.plot_with_window(fill_polygons=False)
plt.show()

# Run nextnano++ and load the 2D results.
charge2d.NNinputf.execute(show_log=False)
charge2d.load_output_data()

# Colour-map a couple of quantities across the bias sweep.
charge2d.plot_results(V_idx=[0, len(charge2d.V) - 1], quantity="N")
charge2d.plot_results(V_idx=[0, len(charge2d.V) - 1], quantity="Ey")
plt.show()

# Transfer the interpolated 2D fields onto device.charge for the EO stage.
charge2d.transfer_results_to_device()

print("device.charge keys:", list(device.charge.keys()))
print("bias points:", np.asarray(device.charge["V"]))
