from __future__ import annotations

import uuid
from dataclasses import dataclass, field
from typing import (
    Callable, Dict, Optional, Union)

import numpy as np
from scipy.interpolate import interp1d
from shapely.geometry import MultiPolygon, Polygon


@dataclass
class SemiconductorPolygon:
    """
    Base class for a polygon with semiconductor properties.

    Attributes:
        polygon (Polygon | Multipolygon): The polygon geometry. It must be a `shapely`_ object, either a `Polygon` or a `MultiPolygon`. 
        
        name (str | None): The name of the polygon. If not provided, a unique ID will be generated.
        
        optical_material (str | complex): The optical material properties of the polygon. This will be deprecated in the future.
        
        rf_eps (float | complex | int | np.ndarray | list): The relative permittivity of the polygon. Material information for RF simulations if charge transport simulation. Internally, it will be converted to a callable function of the type f(omega). If you input a float, then it is constant for all frequencies, and if you give it a ``np.array`` of shape `(2,N)` then it will interpolate to give a callable of the same signature. Note that the charge transport data will only be used to retrieve the conductivity values. You must always insert the real part of the permitivity. Otherwise it will consider it 1.
        
        electro_optic_module (None | ElectroOpticalModel): The electro-optic module associated with the polygon. This is the model that will later be used by :class:`ElectroOpticalSimulator <imodulator.ElectroOpticalSimulator.ElectroOpticalSimulator>` to calculate the electro-optic effect.
        
        electro_optic_module_kwargs (dict): The keyword arguments for the electro-optic module.
        
        calculate_current (bool): A flag to indicate whether to calculate current.
        
        d_buffer_current (bool): The buffer distance for current calculation. The current is calculated with a line integral around the object. This parameter sets the distance from the object where the line integral is calculated. The default value is 0.1 micrometers.
        
        eo_mesh_settings (dict[str, float]): The mesh settings for electro-optic simulations. It is a dictionary with the following keys:

            * **SizeMax** (float): the maximum size of the mesh at a ``distance`` from the object;
            * **resolution** (float): the maximum size of the mesh *inside* the polygon;
            * **distance** (float): the distance from the object where the mesh size is ``SizeMax``;
            * **distance_junction** (float): If there is a boundary of the type Semiconductor-Semiconductor, Semiconductor-Insulator a new polygon can be automatically inserted that captures the junction. This new polygon is inserted only in the direction of the semiconductor. ``distance_junction`` controls how deep into the semiconductor it goes;
            * **resolution_junction** (float): the maximum size of the mesh *inside* the junction;

            .. warning::
                The algorithm to search for junctions is not active at the moment. It is a feature that will be implemented in the future.
        
        rf_mesh_settings (dict[str, float]): The mesh settings for RF simulations. It is a dictionary with the same keys as ``eo_mesh_settings``.
        
        optical_mesh_settings (dict[str, float]): The mesh settings for optical simulations. It is a dictionary with the same keys as ``eo_mesh_settings``.
        
        charge_transport_simulator_kwargs (dict | None): Additional keyword arguments for the charge transport simulator.
            charge_transport_simulator_kwargs;   
                    "material_definition":None, If not defined, default is "Ga(x)In(1-x)As(y)P(1-y)" should be available in Nextnano
                    "doping_conc": None,
                    "doping_type": None,
                    "alloy_y": None,
                    "alloy_x":None,
        
    """

    polygon: Union[Polygon, MultiPolygon]
    name: str = field(default_factory=lambda: str(uuid.uuid4()))
    optical_material: Union[str, complex] = 1 + 0j
    rf_eps: Union[float, complex, int, np.ndarray, list, Callable[[float], float]] = 1
    electro_optic_module: None = None
    electro_optic_module_kwargs: Dict = field(default_factory=dict)
    calculate_current: bool = False
    d_buffer_current: float = 0.1
    charge_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
        "resolution_junction": 0,
        "distance_junction": 0,
    })
    eo_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
        "resolution_junction": 0,
        "distance_junction": 0,
    })
    rf_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
        "resolution_junction": 0,
        "distance_junction": 0,
    })
    optical_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
        "resolution_junction": 0,
        "distance_junction": 0,
    })
    charge_transport_simulator_kwargs: Union[Dict, None] = field(default_factory=lambda:{
        "material_definition":None,
        "doping_conc": None,
        "doping_type": None,
        "alloy_y": None,
        "alloy_x":None,
    })
    has_charge_transport_data: bool = False
    _rf_eps_func: Callable[[float], float] = field(repr=False, init=False)

    def __post_init__(self):
        # Generate unique ID if name is not provided
        if self.name is None:
            self.name = str(uuid.uuid4())

        self._make_rf_eps_callable(self.rf_eps)

        default_mesh_settings = {
            "resolution": 100,
            "SizeMax": 100,
            "distance": 0.0,
            "resolution_junction": 100,
            "distance_junction": 0,
        }

        for key, value in default_mesh_settings.items():
            self.eo_mesh_settings.setdefault(key, value)
            self.charge_mesh_settings.setdefault(key, value)
            self.rf_mesh_settings.setdefault(key, value)
            self.optical_mesh_settings.setdefault(key, value)


    def _make_rf_eps_callable(self, value):
        if isinstance(value, (float, complex, int)):
            super().__setattr__("_rf_eps_func", lambda omega: value)
            super().__setattr__("rf_eps", lambda omega: value)
        elif isinstance(value, np.ndarray) or isinstance(value, list):
            value = np.array(value)
            if np.shape(value)[0] != 2:
                raise ValueError("If rf_eps is given as an array or list, it must have shape (2,N)")
            call = lambda omega: interp1d(2 * np.pi * value[0].real, value[1].real, kind='cubic')(omega) + 1j*interp1d(2 * np.pi * value[0].real, value[1].imag, kind='cubic')(omega)

            super().__setattr__("_rf_eps_func", call)
            super().__setattr__("rf_eps", call)

        elif callable(value):
            super().__setattr__("_rf_eps_func", value)
            super().__setattr__("rf_eps", value)
        else:
            raise ValueError("rf_eps must be a float, complex, int, np.ndarray, or callable")

    def __setattr__(self, name, value):
        # Intercept assignments to rf_eps
        if name == "rf_eps":
            # For initialization, ensure _rf_eps_func also updated
            super().__setattr__(name, value)
            if hasattr(self, "_rf_eps_func"):  # skip during init until post_init
                self._make_rf_eps_callable(value)
        else:
            super().__setattr__(name, value)

@dataclass
class MetalPolygon:
    """
    Initialize a ``MetalPolygon`` object with the specified parameters. All the input spatial dimensions are assumed to be in micrometers, and the frequencies are in GHz.

    Attributes:
        polygon (Polygon | MultiPolygon): 
            The polygon geometry. It must be a `shapely`_ object, either a `Polygon` or a `MultiPolygon`. 

        name (str, optional): 
            The name of the polygon. If not provided, a unique ID will be generated.

        optical_material (str or complex): 
            The optical material properties of the polygon. This will be deprecated in the future.

        rf_eps (float, complex, int, or np.ndarray): 
            The relative permittivity of the polygon. Material information for RF simulations if charge transport simulation. Internally, it will be converted to a callable function of the type f(omega). If you input a float, then it is constant for all frequencies, and if you give it a ``np.array`` of shape `(2,N)` then it will interpolate to give a callable of the same signature.

        calculate_current (bool): 
            A flag to indicate whether to calculate current.

        d_buffer_current (float): 
            The buffer distance for current calculation. The current is calculated with a line integral around the object. This parameter sets the distance from the object where the line integral is calculated.

        eo_mesh_settings (dict[str, float]): 
            The mesh settings for electro-optic simulations. It is a dictionary with the following keys:

                * **SizeMax** (float): the maximum size of the mesh at a ``distance`` from the object;
                * **resolution** (float): the maximum size of the mesh *inside* the polygon;
                * **distance** (float): the distance from the object where the mesh size is ``SizeMax``;
                * **distance_junction** (float): If there is a boundary of the type Semiconductor-Semiconductor, Semiconductor-Insulator a new polygon can be automatically inserted that captures the junction. This new polygon is inserted only in the direction of the semiconductor. ``distance_junction`` controls how deep into the semiconductor it goes;
                * **resolution_junction** (float): the maximum size of the mesh *inside* the junction;

                .. warning::
                    The algorithm to search for junctions is not active at the moment. It is a feature that will be implemented in the future.

        rf_mesh_settings (dict[str, float]): 
            The mesh settings for RF simulations. It is a dictionary with the same keys as ``eo_mesh_settings``.

        optical_mesh_settings (dict[str, float]): 
            The mesh settings for optical simulations. It is a dictionary with the same keys as ``eo_mesh_settings``.

        charge_transport_simulator_kwargs (dict, optional): 
            You don't input metals into charge transport simulation; the simulation line defines the boundaries.
    """
    polygon: Union[Polygon, MultiPolygon] 
    name: Optional[str] = field(default=None)
    optical_material: Union[str, complex] = field(default=1 + 0j)
    calculate_current: bool = field(default=True)
    d_buffer_current: float = field(default=0.1)
    rf_eps: Union[float, complex, int, np.ndarray, list, Callable[[float], float]] = 1
    
    charge_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
        "resolution_junction": 0,
        "distance_junction": 0,
    })
        
    eo_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
    })
    
    rf_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
    })
    
    optical_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
        "resolution_junction": 0.0,
        "distance_junction": 0.0,
    })
    
    _rf_eps_func: Callable[[float], float] = field(repr=False, init=False)

    def __post_init__(self):
        # Generate unique ID if name is not provided
        if self.name is None:
            self.name = str(uuid.uuid4())
        
        self._make_rf_eps_callable(self.rf_eps)

        # Default mesh settings if keys are missing
        default_mesh_settings = {
            "resolution": 100,
            "SizeMax": 100,
            "distance": 0.0,
            "resolution_junction": 100,
            "distance_junction": 0,
        }

        for key, value in default_mesh_settings.items():
            self.eo_mesh_settings.setdefault(key, value)
            self.charge_mesh_settings.setdefault(key, value)
            self.rf_mesh_settings.setdefault(key, value)
            self.optical_mesh_settings.setdefault(key, value)

    def _make_rf_eps_callable(self, value):
        if isinstance(value, (float, complex, int)):
            super().__setattr__("_rf_eps_func", lambda omega: value)
            super().__setattr__("rf_eps", lambda omega: value)
        elif isinstance(value, np.ndarray) or isinstance(value, list):
            value = np.array(value)
            if np.shape(value)[0] != 2:
                raise ValueError("If rf_eps is given as an array or list, it must have shape (2,N)")
            
            call = lambda omega: interp1d(2 * np.pi * value[0].real, value[1].real, kind='cubic')(omega) + 1j*interp1d(2 * np.pi * value[0].real, value[1].imag, kind='cubic')(omega)

            super().__setattr__("_rf_eps_func", call)
            super().__setattr__("rf_eps", call)

        elif callable(value):
            super().__setattr__("_rf_eps_func", value)
            super().__setattr__("rf_eps", value)
        else:
            raise ValueError("rf_eps must be a float, complex, int, np.ndarray, or callable")

    def __setattr__(self, name, value):
        # Intercept assignments to rf_eps
        if name == "rf_eps":
            # For initialization, ensure _rf_eps_func also updated
            super().__setattr__(name, value)
            if hasattr(self, "_rf_eps_func"):  # skip during init until post_init
                self._make_rf_eps_callable(value)
        else:
            super().__setattr__(name, value)


@dataclass
class InsulatorPolygon:
    """
    Initialize an ``InsulatorPolygon`` object with the specified parameters. All the input spatial dimensions are assumed to be in micrometers, and the frequencies are in GHz.

    Attributes:
        polygon (Polygon | MultiPolygon): 
            The polygon geometry. It must be a `shapely`_ object, either a `Polygon` or a `MultiPolygon`. 

        name (str, optional): 
            The name of the polygon. If not provided, a unique ID will be generated.

        optical_material (str or complex): 
            The optical material properties of the polygon. This will be deprecated in the future.

        rf_eps (float, complex, int, or np.ndarray): 
            The relative permittivity of the polygon. Material information for RF simulations if charge transport simulation. Internally, it will be converted to a callable function of the type f(omega). If you input a float, then it is constant for all frequencies, and if you give it a ``np.array`` of shape `(2,N)` then it will interpolate to give a callable of the same signature.

        electro_optic_module (None): 
            The electro-optic module associated with the polygon.

        eo_mesh_settings (dict[str, float]): 
            The mesh settings for electro-optic simulations. It is a dictionary with the following keys:

                * **SizeMax** (float): the maximum size of the mesh at a ``distance`` from the object;
                * **resolution** (float): the maximum size of the mesh *inside* the polygon;
                * **distance** (float): the distance from the object where the mesh size is ``SizeMax``;
                * **distance_junction** (float): If there is a boundary of the type Semiconductor-Semiconductor, Semiconductor-Insulator a new polygon can be automatically inserted that captures the junction. This new polygon is inserted only in the direction of the semiconductor. ``distance_junction`` controls how deep into the semiconductor it goes;
                * **resolution_junction** (float): the maximum size of the mesh *inside* the junction;

                .. warning::
                    The algorithm to search for junctions is not active at the moment. It is a feature that will be implemented in the future.

        rf_mesh_settings (dict[str, float]): 
            The mesh settings for RF simulations. It is a dictionary with the same keys as ``eo_mesh_settings``.

        optical_mesh_settings (dict[str, float]): 
            The mesh settings for optical simulations. It is a dictionary with the same keys as ``eo_mesh_settings``.

        charge_transport_simulator_kwargs (dict, optional): 
            Additional keyword arguments for the charge transport simulator.

        calculate_current (bool): 
            Flag to indicate whether to calculate current.

        d_buffer_current (bool): 
            Buffer distance for current calculation.
    """
    polygon: Union[Polygon, MultiPolygon]
    name: str = field(default_factory=lambda: str(uuid.uuid4()))
    optical_material: Union[str, complex] = 1 + 0j
    rf_eps: Union[float, complex, int, np.ndarray, list, Callable[[float], float]] = 1 
    electro_optic_module: None = None
    charge_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
        "resolution_junction": 0,
        "distance_junction": 0,
    })
    eo_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
    })
    rf_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
    })
    optical_mesh_settings: Dict[str, float] = field(default_factory=lambda: {
        "resolution": 0.1,
        "SizeMax": 0.1,
        "distance": 0.1,
        "resolution_junction": 0,
        "distance_junction": 0,
    })
    charge_transport_simulator_kwargs: Union[Dict, None] = None
    calculate_current: bool = False
    d_buffer_current: bool = False
    
    _rf_eps_func: Callable[[float], float] = field(repr=False, init=False)
    def __post_init__(self):
        if self.name is None:
            self.name = str(uuid.uuid4())
        
        self._make_rf_eps_callable(self.rf_eps)

        default_mesh_settings = {
            "resolution": 100,
            "SizeMax": 100,
            "distance": 0.0,
            "resolution_junction": 100,
            "distance_junction": 0,
        }

        for key, value in default_mesh_settings.items():
            self.eo_mesh_settings.setdefault(key, value)
            self.charge_mesh_settings.setdefault(key, value)
            self.rf_mesh_settings.setdefault(key, value)
            self.optical_mesh_settings.setdefault(key, value)

    def _make_rf_eps_callable(self, value):
        if isinstance(value, (float, complex, int)):
            super().__setattr__("_rf_eps_func", lambda omega: value)
            super().__setattr__("rf_eps", lambda omega: value)
        elif isinstance(value, np.ndarray) or isinstance(value, list):
            value = np.array(value)
            if np.shape(value)[0] != 2:
                raise ValueError("If rf_eps is given as an array or list, it must have shape (2,N)")

            call = lambda omega: interp1d(2 * np.pi * value[0].real, value[1].real, kind='linear')(omega) + 1j*interp1d(2 * np.pi * value[0].real, value[1].imag, kind='linear')(omega)

            super().__setattr__("_rf_eps_func", call)
            super().__setattr__("rf_eps", call)

        elif callable(value):
            super().__setattr__("_rf_eps_func", value)
            super().__setattr__("rf_eps", value)
        else:
            raise ValueError("rf_eps must be a float, complex, int, np.ndarray, or callable")

    def __setattr__(self, name, value):
        # Intercept assignments to rf_eps
        if name == "rf_eps":
            # For initialization, ensure _rf_eps_func also updated
            super().__setattr__(name, value)
            if hasattr(self, "_rf_eps_func"):  # skip during init until post_init
                self._make_rf_eps_callable(value)
        else:
            super().__setattr__(name, value)