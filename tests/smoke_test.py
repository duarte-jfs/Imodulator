"""End-to-end smoke tests for imodulator.

Run with a plain interpreter -- there is no test framework here on purpose:

    python tests/smoke_test.py

Part A checks that every module imports on a machine with no config.yaml and no
proprietary solvers. Part B runs the FOSS half of the dc_electro_optical
tutorial (femwell mode solve -> solcore PDD -> EO response) and asserts the
results are physically sane.

Part B is derived from docs/Tutorials/dc_electro_optical, reduced for CI: the
20-iteration mesh refinement loop is dropped and the bias sweep is cut from 21
points to 2, which is the minimum that still yields a non-zero EO response.
"""

import importlib
import sys

import matplotlib

# Must precede any imodulator import: several modules pull matplotlib at import
# time and CI runners have no display.
matplotlib.use("Agg")

import numpy as np
import openbandparams as obp
import shapely

import imodulator
from imodulator.ElectroOpticalModel import InGaAsPElectroOpticalModel

# ---------------------------------------------------------------------------
# Geometry helpers copied verbatim from the first code cell of
# docs/Tutorials/dc_electro_optical/dc_electro_optical.ipynb. Kept byte-for-byte
# so the smoke test tracks the tutorial rather than a re-derived device.
# ---------------------------------------------------------------------------


def tand_fitted_bcb(x):
    """
    Fitted to results from https://link.springer.com/article/10.1007/s10762-009-9552-0

    x must be in GHz
    """
    out = 0.0093839 - 0.01790336 * np.exp(-0.04773444 * (x - -4.64170761))

    if isinstance(x, (list, np.ndarray)):
        x = np.asarray(x)

        out[np.where(out < 0.001)] = 0.001
    else:
        if out < 0.001:
            out = 0.001
    return out


class InP_EOPM:
    def __init__(self, **kwargs):

        self.e = 1.60e-19  # electron charge in C
        self.e0 = 8.85e-12  # vacuum permittivity in F/m

        self.w_sig_metal = 5  # Width of signal metal in um
        self.metal_sep = 10  # Separation between signal and ground metals in um
        self.h_metal = 4  # Height of metals in um
        self.w_gnd_metal = 10

        self.w_wg = 1
        self.h_n = 0.4
        self.h_wg1 = 0.5
        self.h_wg2 = 0.3
        self.h_p1 = 1
        self.h_p2 = 0.2

        self.h_box = 4

        self.w_window = 100
        self.h_bottom = 30
        self.h_top = 30

        for kwarg, value in kwargs.items():
            if hasattr(self, kwarg):
                setattr(self, kwarg, value)

    def _make_meshes(self):
        # optical mesh
        self.optical_mesh_settings = {
            "substrate": {"resolution": 0.5, "SizeMax": 5, "distance": 0.1},
            "background": {"resolution": 0.5, "SizeMax": 5, "distance": 0.1},
            "box": {"resolution": 0.3, "SizeMax": 0.2, "distance": 0.1},
            "sig_metal": {"resolution": 10, "SizeMax": 0.2, "distance": 0.1},
            "n_metal_left": {"resolution": 10, "SizeMax": 0.2, "distance": 0.1},
            "n_metal_right": {"resolution": 10, "SizeMax": 0.2, "distance": 0.1},
            "bcb": {"resolution": 0.5, "SizeMax": 5, "distance": 0.1},
            "n": {"resolution": 0.5, "SizeMax": 5, "distance": 0.1},
            "wg1": {"resolution": 0.1, "SizeMax": 5, "distance": 0.1},
            "wg2": {"resolution": 0.1, "SizeMax": 5, "distance": 0.1},
            "p1": {"resolution": 0.5, "SizeMax": 5, "distance": 0.1},
            "p2": {"resolution": 0.5, "SizeMax": 5, "distance": 0.1},
        }

        # RF mesh
        self.rf_mesh_settings = {
            "substrate": {"resolution": 5, "SizeMax": 5, "distance": 0.1},
            "background": {"resolution": 5, "SizeMax": 5, "distance": 0.1},
            "box": {"resolution": 3, "SizeMax": 3, "distance": 0.1},
            "sig_metal": {"resolution": 3, "SizeMax": 3, "distance": 0.1},
            "n_metal_left": {"resolution": 3, "SizeMax": 3, "distance": 0.1},
            "n_metal_right": {"resolution": 3, "SizeMax": 3, "distance": 0.1},
            "bcb": {"resolution": 5, "SizeMax": 5, "distance": 0.1},
            "n": {"resolution": 5, "SizeMax": 5, "distance": 0.1},
            "wg1": {"resolution": 1, "SizeMax": 5, "distance": 0.1},
            "wg2": {"resolution": 1, "SizeMax": 5, "distance": 0.1},
            "p1": {"resolution": 5, "SizeMax": 5, "distance": 0.1},
            "p2": {"resolution": 5, "SizeMax": 5, "distance": 0.1},
        }

        # eo mesh
        self.eo_mesh_settings = {
            "substrate": {"resolution": 0.5, "SizeMax": 5, "distance": 0.1},
            "background": {"resolution": 0.5, "SizeMax": 5, "distance": 0.1},
            "box": {"resolution": 0.3, "SizeMax": 0.2, "distance": 0.1},
            "sig_metal": {"resolution": 10, "SizeMax": 0.2, "distance": 0.1},
            "n_metal_left": {"resolution": 10, "SizeMax": 0.2, "distance": 0.1},
            "n_metal_right": {"resolution": 10, "SizeMax": 0.2, "distance": 0.1},
            "bcb": {"resolution": 0.5, "SizeMax": 5, "distance": 0.1},
            "n": {"resolution": 0.1, "SizeMax": 5, "distance": 0.1},
            "wg1": {"resolution": 0.05, "SizeMax": 5, "distance": 0.1},
            "wg2": {"resolution": 0.05, "SizeMax": 5, "distance": 0.1},
            "p1": {"resolution": 0.05, "SizeMax": 5, "distance": 0.1},
            "p2": {"resolution": 0.5, "SizeMax": 5, "distance": 0.1},
        }

        self.charge_mesh_settings = {
            "substrate": {"resolution": 0.5},
            "background": {"resolution": 0.5},
            "sig_metal": {"resolution": 0.01},
            "n_metal_left": {"resolution": 0.01},
            "n_metal_right": {"resolution": 0.01},
            "n": {"resolution": 0.003},
            "wg1": {"resolution": 0.002},
            "wg2": {"resolution": 0.002},
            "p1": {"resolution": 0.003},
            "p2": {"resolution": 0.003},
        }

    def _create_polygons(self):
        # We will now set the RF properties of the metals and the BCB
        freq = np.linspace(0.1, 100, 100)  # GHz. This will be the simulation frequency

        eps_rf_metal = 1 - 1j * 6e7 / (2 * np.pi * freq * 1e9 * self.e0)
        eps_rf_metal = np.asarray([freq, eps_rf_metal])

        bcb_eps_real = 2.65 * np.ones(100)
        bcb_eps_imag = bcb_eps_real * tand_fitted_bcb(freq)

        bcb_eps = bcb_eps_real - 1j * bcb_eps_imag
        bcb_eps = np.asarray([freq, bcb_eps])

        # Now we create the PhotoPolygons
        self.substrate = imodulator.SemiconductorPolygon(
            shapely.box(
                -self.w_window / 2, -self.h_box - self.h_bottom, self.w_window / 2, -self.h_box
            ),
            rf_eps=11.7,
            name="substrate",
            optical_material=3**2,
            eo_mesh_settings=self.eo_mesh_settings["substrate"],
            rf_mesh_settings=self.rf_mesh_settings["substrate"],
            optical_mesh_settings=self.optical_mesh_settings["substrate"],
        )

        self.background = imodulator.InsulatorPolygon(
            shapely.box(
                -self.w_window / 2,
                -self.h_box - self.h_bottom,
                self.w_window / 2,
                self.h_n
                + self.h_wg1
                + self.h_wg2
                + self.h_p1
                + self.h_p2
                + self.h_metal
                + self.h_top,
            ),
            rf_eps=1,
            optical_material=1,
            eo_mesh_settings=self.eo_mesh_settings["background"],
            rf_mesh_settings=self.rf_mesh_settings["background"],
            optical_mesh_settings=self.optical_mesh_settings["background"],
            name="background",
        )

        self.box = imodulator.InsulatorPolygon(
            shapely.box(-self.w_window / 2, -self.h_box, self.w_window / 2, 0),
            rf_eps=3.9 - 1j * 3.9 * 0.001,
            optical_material=1.44**2,
            eo_mesh_settings=self.eo_mesh_settings["box"],
            rf_mesh_settings=self.rf_mesh_settings["box"],
            optical_mesh_settings=self.optical_mesh_settings["box"],
            name="box",
        )

        n_obp_material = obp.GaInPAs(T=300, As=0, a=obp.InP.a())
        self.n = imodulator.SemiconductorPolygon(
            shapely.box(
                -self.w_sig_metal / 2 - self.metal_sep - self.w_gnd_metal,
                0,
                self.w_sig_metal / 2 + self.metal_sep + self.w_gnd_metal,
                self.h_n,
            ),
            rf_eps=n_obp_material.dielectric(T=300),
            optical_material=n_obp_material.refractive_index(T=300) ** 2,
            eo_mesh_settings=self.eo_mesh_settings["n"],
            rf_mesh_settings=self.rf_mesh_settings["n"],
            optical_mesh_settings=self.optical_mesh_settings["n"],
            charge_mesh_settings=self.charge_mesh_settings["n"],
            name="n",
            electro_optic_module=InGaAsPElectroOpticalModel,
            electro_optic_module_kwargs={"y": 0, "T": 300, "BF_model": "vinchant"},
            charge_transport_simulator_kwargs={"sol_obp_material": n_obp_material, "sol_Nd": 1e18},
        )

        wg1_obp_material = obp.GaInPAs(T=300, As=0.53, a=obp.InP.a())
        self.wg1 = imodulator.SemiconductorPolygon(
            shapely.box(-self.w_wg / 2, self.h_n, self.w_wg / 2, self.h_n + self.h_wg1),
            rf_eps=wg1_obp_material.dielectric(T=300),
            optical_material=wg1_obp_material.refractive_index(T=300) ** 2,
            eo_mesh_settings=self.eo_mesh_settings["wg1"],
            rf_mesh_settings=self.rf_mesh_settings["wg1"],
            optical_mesh_settings=self.optical_mesh_settings["wg1"],
            charge_mesh_settings=self.charge_mesh_settings["wg1"],
            name="wg1",
            electro_optic_module=InGaAsPElectroOpticalModel,
            electro_optic_module_kwargs={"y": 0.53, "T": 300, "BF_model": "vinchant"},
            charge_transport_simulator_kwargs={
                "sol_obp_material": wg1_obp_material,
                "sol_Nd": 1e16,
            },
        )

        wg2_obp_material = obp.GaInPAs(T=300, As=0, a=obp.InP.a())
        self.wg2 = imodulator.SemiconductorPolygon(
            shapely.box(
                -self.w_wg / 2,
                self.h_n + self.h_wg1,
                self.w_wg / 2,
                self.h_n + self.h_wg1 + self.h_wg2,
            ),
            rf_eps=wg2_obp_material.dielectric(T=300),
            optical_material=wg2_obp_material.refractive_index(T=300) ** 2,
            eo_mesh_settings=self.eo_mesh_settings["wg2"],
            rf_mesh_settings=self.rf_mesh_settings["wg2"],
            optical_mesh_settings=self.optical_mesh_settings["wg2"],
            charge_mesh_settings=self.charge_mesh_settings["wg2"],
            name="wg2",
            electro_optic_module=InGaAsPElectroOpticalModel,
            electro_optic_module_kwargs={"y": 0, "T": 300, "BF_model": "vinchant"},
            charge_transport_simulator_kwargs={
                "sol_obp_material": wg2_obp_material,
                "sol_Nd": 1e16,
            },
        )

        p1_obp_material = obp.GaInPAs(T=300, As=0, a=obp.InP.a())
        self.p1 = imodulator.SemiconductorPolygon(
            shapely.box(
                -self.w_wg / 2,
                self.h_n + self.h_wg1 + self.h_wg2,
                self.w_wg / 2,
                self.h_n + self.h_wg1 + self.h_wg2 + self.h_p1,
            ),
            rf_eps=p1_obp_material.dielectric(T=300),
            optical_material=p1_obp_material.refractive_index(T=300) ** 2,
            eo_mesh_settings=self.eo_mesh_settings["p1"],
            rf_mesh_settings=self.rf_mesh_settings["p1"],
            optical_mesh_settings=self.optical_mesh_settings["p1"],
            charge_mesh_settings=self.charge_mesh_settings["p1"],
            name="p1",
            electro_optic_module=InGaAsPElectroOpticalModel,
            electro_optic_module_kwargs={"y": 0, "T": 300, "BF_model": "vinchant"},
            charge_transport_simulator_kwargs={
                "sol_obp_material": wg2_obp_material,
                "sol_Na": 1e17,
            },
        )

        p2_obp_material = obp.GaInAs(T=300, a=obp.InP.a())
        self.p2 = imodulator.SemiconductorPolygon(
            shapely.box(
                -self.w_wg / 2,
                self.h_n + self.h_wg1 + self.h_wg2 + self.h_p1,
                self.w_wg / 2,
                self.h_n + self.h_wg1 + self.h_wg2 + self.h_p1 + self.h_p2,
            ),
            rf_eps=p2_obp_material.dielectric(T=300),
            optical_material=p2_obp_material.refractive_index(T=300) ** 2,
            eo_mesh_settings=self.eo_mesh_settings["p2"],
            rf_mesh_settings=self.rf_mesh_settings["p2"],
            optical_mesh_settings=self.optical_mesh_settings["p2"],
            charge_mesh_settings=self.charge_mesh_settings["p2"],
            name="p2",
            electro_optic_module=InGaAsPElectroOpticalModel,
            electro_optic_module_kwargs={"y": 0, "T": 300, "BF_model": "vinchant"},
            charge_transport_simulator_kwargs={"sol_obp_material": p2_obp_material, "sol_Na": 1e19},
        )

        self.bcb_far_left = imodulator.InsulatorPolygon(
            shapely.box(
                -self.w_window / 2,
                0,
                -self.w_sig_metal / 2 - self.metal_sep - self.w_gnd_metal,
                self.h_n + self.h_wg1 + self.h_wg2 + self.h_p1 + self.h_p2,
            ),
            rf_eps=bcb_eps,
            optical_material=1.56**2,
            eo_mesh_settings=self.eo_mesh_settings["bcb"],
            rf_mesh_settings=self.rf_mesh_settings["bcb"],
            optical_mesh_settings=self.optical_mesh_settings["bcb"],
            name="bcb_far_left",
        )

        self.bcb_far_right = imodulator.InsulatorPolygon(
            shapely.box(
                self.metal_sep + self.w_gnd_metal + self.w_sig_metal / 2,
                0,
                self.w_window / 2,
                self.h_n + self.h_wg1 + self.h_wg2 + self.h_p1 + self.h_p2,
            ),
            rf_eps=bcb_eps,
            optical_material=1.56**2,
            eo_mesh_settings=self.eo_mesh_settings["bcb"],
            rf_mesh_settings=self.rf_mesh_settings["bcb"],
            optical_mesh_settings=self.optical_mesh_settings["bcb"],
            name="bcb_far_right",
        )

        self.bcb_left = imodulator.InsulatorPolygon(
            shapely.box(
                -self.w_sig_metal / 2 - self.metal_sep,
                self.h_n,
                -self.w_wg / 2,
                self.h_n + self.h_wg1 + self.h_wg2 + self.h_p1 + self.h_p2,
            ),
            rf_eps=bcb_eps,
            optical_material=1.56**2,
            eo_mesh_settings=self.eo_mesh_settings["bcb"],
            rf_mesh_settings=self.rf_mesh_settings["bcb"],
            optical_mesh_settings=self.optical_mesh_settings["bcb"],
            name="bcb_left",
        )

        self.bcb_right = imodulator.InsulatorPolygon(
            shapely.box(
                self.w_wg / 2,
                self.h_n,
                self.w_sig_metal / 2 + self.metal_sep,
                self.h_n + self.h_wg1 + self.h_wg2 + self.h_p1 + self.h_p2,
            ),
            rf_eps=bcb_eps,
            optical_material=1.56**2,
            eo_mesh_settings=self.eo_mesh_settings["bcb"],
            rf_mesh_settings=self.rf_mesh_settings["bcb"],
            optical_mesh_settings=self.optical_mesh_settings["bcb"],
            name="bcb_right",
        )

        self.sig_metal = imodulator.MetalPolygon(
            shapely.box(
                -self.w_sig_metal / 2,
                self.h_n + self.h_wg1 + self.h_wg2 + self.h_p1 + self.h_p2,
                self.w_sig_metal / 2,
                self.h_n + self.h_wg1 + self.h_wg2 + self.h_p1 + self.h_p2 + self.h_metal,
            ),
            rf_eps=eps_rf_metal,
            eo_mesh_settings=self.eo_mesh_settings["sig_metal"],
            rf_mesh_settings=self.rf_mesh_settings["sig_metal"],
            optical_mesh_settings=self.optical_mesh_settings["sig_metal"],
            name="sig_metal",
            calculate_current=True,
            d_buffer_current=min(self.w_sig_metal / 20, self.h_metal / 20, 0.05),
        )

        self.n_metal_left = imodulator.MetalPolygon(
            shapely.box(
                -self.w_sig_metal / 2 - self.metal_sep - self.w_gnd_metal,
                self.h_n,
                -self.w_sig_metal / 2 - self.metal_sep,
                self.h_n + self.h_metal,
            ),
            rf_eps=eps_rf_metal,
            eo_mesh_settings=self.eo_mesh_settings["n_metal_left"],
            rf_mesh_settings=self.rf_mesh_settings["n_metal_left"],
            optical_mesh_settings=self.optical_mesh_settings["n_metal_left"],
            name="n_metal_left",
            calculate_current=False,
        )

        self.n_metal_right = imodulator.MetalPolygon(
            shapely.box(
                self.w_sig_metal / 2 + self.metal_sep,
                self.h_n,
                self.w_sig_metal / 2 + self.metal_sep + self.w_gnd_metal,
                self.h_n + self.h_metal,
            ),
            rf_eps=eps_rf_metal,
            eo_mesh_settings=self.eo_mesh_settings["n_metal_right"],
            rf_mesh_settings=self.rf_mesh_settings["n_metal_right"],
            optical_mesh_settings=self.optical_mesh_settings["n_metal_right"],
            name="n_metal_right",
            calculate_current=False,
        )

    def _initialize_device(self):
        photo_polygons = [
            self.sig_metal,
            self.n_metal_left,
            self.n_metal_right,
            self.p2,
            self.p1,
            self.wg2,
            self.wg1,
            self.n,
            self.box,
            self.bcb_left,
            self.bcb_right,
            self.bcb_far_left,
            self.bcb_far_right,
            self.substrate,
            self.background,
        ]

        # Just in case there are empty polygons
        idxs_to_remove = []
        for i, poly in enumerate(photo_polygons):
            if np.isclose(poly.polygon.bounds[1], poly.polygon.bounds[3]):
                idxs_to_remove.append(i)
        for i in idxs_to_remove[::-1]:
            del photo_polygons[i]
        self.device = imodulator.PhotonicDevice(photo_polygons)


MODULES = [
    "Config",
    "_optional_deps",
    "PhotonicPolygon",
    "PhotonicDevice",
    "ElectroOpticalModel",
    "OpticalSimulator",
    "ChargeSimulator",
    "ElectroOpticalSimulator",
    "RFSimulator",
]


def test_imports():
    """Every module must import without a config.yaml present.

    This is not as trivial as it looks. Config.py instantiates config_instance at
    module scope, OpticalSimulator and ChargeSimulator call get_lumapi() /
    get_nextnanopy() at module scope, and ChargeSimulator imports gmsh at module
    scope -- so a successful import also proves the gmsh native-library fixup in
    flake.nix worked. All of these are expected to warn, none to raise.
    """
    for name in MODULES:
        importlib.import_module(f"imodulator.{name}")

    for name in imodulator.__all__:
        assert hasattr(imodulator, name), f"__all__ exports missing name {name!r}"


def test_electro_optic_pipeline():
    """femwell mode solve -> solcore PDD -> electro-optic response."""
    from imodulator.ChargeSimulator import ChargeSimulatorSolcore
    from imodulator.ElectroOpticalSimulator import ElectroOpticalSimulator
    from imodulator.OpticalSimulator import OpticalSimulatorFEMWELL

    eopm = InP_EOPM()
    eopm._make_meshes()
    eopm._create_polygons()
    eopm._initialize_device()

    mode = OpticalSimulatorFEMWELL(
        device=eopm.device,
        simulation_window=shapely.box(-4, -1, 4, 4),
        include_metals=False,
    )
    mode.make_mesh()
    modes = mode.compute_modes(
        wavelength=1.55,
        num_modes=5,
        order=1,
        metallic_boundaries=False,
        n_guess=3.2,
        return_modes=True,
    )

    n_eff = max(m.n_eff.real for m in modes)
    assert 2.8 < n_eff < 3.5, f"guided-mode n_eff out of physical range: {n_eff}"

    mode.transfer_results_to_device(TE_TM_idx=[1, 0])

    line = shapely.LineString(
        [
            [0, 0.01],
            [0, eopm.h_n + eopm.h_wg1 + eopm.h_wg2 + eopm.h_p1 + eopm.h_p2 - 0.01],
        ]
    )
    charge = ChargeSimulatorSolcore(
        device=eopm.device,
        simulation_line=line,
        bias_start_stop_step=[0, 6, 2],
    )
    charge.solve_PDD(verbose=False, tol=1e-6, max_iter=1000, smooth_output=False)
    charge.transfer_results_to_device(dx=0.01)

    voltages = np.asarray(eopm.device.charge["V"])
    assert len(voltages) == 2, f"expected a 2-point bias sweep, got {voltages}"

    EO = ElectroOpticalSimulator(device=eopm.device, simulation_window=shapely.box(-3, -1, 3, 5))
    EO.make_mesh()
    EO.get_epsilon_optical()

    dperms = EO.epsilon_optical["InGaAsPElectroOpticalModel"]["dperms"]
    labels = EO.epsilon_optical["InGaAsPElectroOpticalModel"]["labels"]
    assert dperms.shape[:2] == (3, 3), f"dperms is not a 3x3 tensor field: {dperms.shape}"
    assert dperms.shape[-1] == len(labels), (
        f"dperms has {dperms.shape[-1]} effects but {len(labels)} labels"
    )

    kappa = EO.calculate_EO_response(
        voltage_idx=1,
        rot_x=0,
        rot_y=np.pi / 4,
        rot_z=0,
        base_epsilon_voltage_idx=0,
        optical_mode_a="TE",
        optical_mode_b="TE",
    )
    results = np.asarray(kappa["InGaAsPElectroOpticalModel"]["results"])
    assert np.all(np.isfinite(results)), "EO response contains NaN or inf"
    assert np.any(results != 0), "EO response is identically zero"


def main():
    for test in (test_imports, test_electro_optic_pipeline):
        print(f"--- {test.__name__} ---", flush=True)
        test()
        print(f"--- {test.__name__} PASSED ---", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
