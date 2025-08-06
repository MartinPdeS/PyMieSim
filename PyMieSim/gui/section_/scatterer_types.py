#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Scatterer Types Manager

This module provides a clean interface for managing different scatterer types
and their associated parameters in the PyMieSim GUI.
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Any
from PyMieSim.gui.section_.base import length_units


class BaseScattererType(ABC):
    """
    Abstract base class for scatterer types.

    This class defines the interface that all scatterer types must implement,
    ensuring consistent behavior across different particle geometries.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the name of the scatterer type."""
        pass

    @property
    @abstractmethod
    def label(self) -> str:
        """Return the display label for the scatterer type."""
        pass

    @property
    @abstractmethod
    def description(self) -> str:
        """Return a description of the scatterer type."""
        pass

    @property
    @abstractmethod
    def use_cases(self) -> List[str]:
        """Return a list of common use cases."""
        pass

    @abstractmethod
    def get_input_definitions(self) -> Dict[str, Dict[str, Any]]:
        """Return the input field definitions for this scatterer type."""
        pass

    @abstractmethod
    def get_parameter_mapping(self) -> Dict[str, str]:
        """Return the mapping from GUI parameters to PyMieSim parameters."""
        pass

    @abstractmethod
    def validate_inputs(self, input_values: Dict[str, str]) -> tuple[bool, str]:
        """Validate the input values for this scatterer type."""
        pass


class SphereScatterer(BaseScattererType):
    """
    Sphere scatterer type implementation.

    Handles homogeneous spherical particles with uniform refractive index.
    """

    @property
    def name(self) -> str:
        return "sphere"

    @property
    def label(self) -> str:
        return "Sphere"

    @property
    def description(self) -> str:
        return "Homogeneous spherical particle with uniform refractive index throughout the volume."

    @property
    def use_cases(self) -> List[str]:
        return [
            "Simple dielectric spheres (glass, plastic)",
            "Metallic nanoparticles (gold, silver)",
            "Biological cells (simplified model)",
            "Water droplets and aerosols",
            "Polystyrene microspheres"
        ]

    def get_input_definitions(self) -> Dict[str, Dict[str, Any]]:
        """Get input field definitions for sphere scatterer."""
        return {
            "diameter": {
                "id": "scatterer-diameter",
                "label": f"Diameter [{length_units}]",
                "default": "100:2000:50",
                "tooltip": "Sphere diameter in nanometers. Use 'start:end:count' for parameter sweeps or comma-separated values for discrete points.",
                "placeholder": "e.g., 500 or 100:1000:20 or 200,400,600"
            },
            "property": {
                "id": "scatterer-property",
                "label": "Refractive Index",
                "default": "1.5",
                "tooltip": "Complex refractive index of sphere material (n+ik). Use real values for dielectrics, complex for metals.",
                "placeholder": "e.g., 1.5 or 1.5+0.1j or 0.5+2.0j"
            },
            "medium_property": {
                "id": "scatterer-medium-property",
                "label": "Medium Refractive Index",
                "default": "1.0",
                "tooltip": "Refractive index of surrounding medium (1.0 for air, 1.33 for water, 1.5 for glass).",
                "placeholder": "e.g., 1.0 (air) or 1.33 (water)"
            }
        }

    def get_parameter_mapping(self) -> Dict[str, str]:
        """Get parameter mapping for PyMieSim Sphere class."""
        return {
            'diameter': 'diameter',
            'property': 'property',
            'medium_property': 'medium_property'
        }

    def validate_inputs(self, input_values: Dict[str, str]) -> tuple[bool, str]:
        """Validate sphere input parameters."""
        if not input_values.get('diameter'):
            return False, "Diameter is required for sphere scatterer"

        if not input_values.get('property'):
            return False, "Refractive index is required for sphere scatterer"

        if not input_values.get('medium_property'):
            return False, "Medium refractive index is required"

        # Additional validation could be added here (e.g., positive diameter values)
        return True, ""


class CoreShellScatterer(BaseScattererType):
    """
    Core-shell scatterer type implementation.

    Handles spherical particles with distinct core and shell materials.
    """

    @property
    def name(self) -> str:
        return "coreshell"

    @property
    def label(self) -> str:
        return "Core-Shell"

    @property
    def description(self) -> str:
        return "Spherical particle with distinct core and shell materials, enabling complex optical responses."

    @property
    def use_cases(self) -> List[str]:
        return [
            "Metal-dielectric nanoparticles (Au core, SiO₂ shell)",
            "Coated microspheres for drug delivery",
            "Biological cells with membrane structure",
            "Composite particles for enhanced properties",
            "Plasmonic enhancement structures"
        ]

    def get_input_definitions(self) -> Dict[str, Dict[str, Any]]:
        """Get input field definitions for core-shell scatterer."""
        return {
            "core_diameter": {
                "id": "scatterer-core-diameter",
                "label": f"Core Diameter [{length_units}]",
                "default": "100:1000:30",
                "tooltip": "Core diameter in nanometers. Must be smaller than total particle diameter.",
                "placeholder": "e.g., 200 or 50:500:20"
            },
            "shell_thickness": {
                "id": "scatterer-shell-thickness",
                "label": f"Shell Thickness [{length_units}]",
                "default": "10:200:20",
                "tooltip": "Shell thickness in nanometers. Total diameter = core_diameter + 2 * shell_thickness.",
                "placeholder": "e.g., 50 or 10:100:10"
            },
            "core_property": {
                "id": "scatterer-core-property",
                "label": "Core Refractive Index",
                "default": "2.4+1.2j",
                "tooltip": "Complex refractive index of core material (typically metal: high imaginary part).",
                "placeholder": "e.g., 2.4+1.2j (Au) or 1.8 (TiO₂)"
            },
            "shell_property": {
                "id": "scatterer-shell-property",
                "label": "Shell Refractive Index",
                "default": "1.5",
                "tooltip": "Complex refractive index of shell material (typically dielectric: low imaginary part).",
                "placeholder": "e.g., 1.5 (SiO₂) or 2.0 (TiO₂)"
            },
            "medium_property": {
                "id": "scatterer-medium-property",
                "label": "Medium Refractive Index",
                "default": "1.0",
                "tooltip": "Refractive index of surrounding medium.",
                "placeholder": "e.g., 1.0 (air) or 1.33 (water)"
            }
        }

    def get_parameter_mapping(self) -> Dict[str, str]:
        """Get parameter mapping for PyMieSim CoreShell class."""
        return {
            'core_diameter': 'core_diameter',
            'shell_thickness': 'shell_thickness',
            'core_property': 'core_property',
            'shell_property': 'shell_property',
            'medium_property': 'medium_property'
        }

    def validate_inputs(self, input_values: Dict[str, str]) -> tuple[bool, str]:
        """Validate core-shell input parameters."""
        required_fields = [
            ('core_diameter', 'Core diameter'),
            ('shell_thickness', 'Shell thickness'),
            ('core_property', 'Core refractive index'),
            ('shell_property', 'Shell refractive index'),
            ('medium_property', 'Medium refractive index')
        ]

        for field, field_name in required_fields:
            if not input_values.get(field):
                return False, f"{field_name} is required for core-shell scatterer"

        # Additional validation could be added here
        # (e.g., core_diameter > 0, shell_thickness > 0)
        return True, ""


class ScattererTypeManager:
    """
    Manager class for handling different scatterer types.

    This class provides a clean interface for registering, accessing, and
    managing different scatterer types in the GUI.
    """

    def __init__(self):
        """Initialize the scatterer type manager."""
        self._types: Dict[str, BaseScattererType] = {}
        self._register_default_types()

    def _register_default_types(self):
        """Register the default scatterer types."""
        self.register_type(SphereScatterer())
        self.register_type(CoreShellScatterer())

    def register_type(self, scatterer_type: BaseScattererType):
        """
        Register a new scatterer type.

        Parameters
        ----------
        scatterer_type : BaseScattererType
            The scatterer type instance to register.
        """
        self._types[scatterer_type.name] = scatterer_type

    def get_type(self, name: str) -> BaseScattererType:
        """
        Get a scatterer type by name.

        Parameters
        ----------
        name : str
            Name of the scatterer type.

        Returns
        -------
        BaseScattererType
            The scatterer type instance.

        Raises
        ------
        KeyError
            If the scatterer type is not found.
        """
        if name not in self._types:
            available = list(self._types.keys())
            raise KeyError(f"Unknown scatterer type: {name}. Available types: {available}")
        return self._types[name]

    def get_available_types(self) -> List[Dict[str, str]]:
        """
        Get list of available scatterer types for dropdown.

        Returns
        -------
        List[Dict[str, str]]
            List of dictionaries with 'label' and 'value' keys for Dash dropdown.
        """
        return [
            {'label': scatterer_type.label, 'value': scatterer_type.name}
            for scatterer_type in self._types.values()
        ]

    def get_available_names(self) -> List[str]:
        """
        Get list of available scatterer type names.

        Returns
        -------
        List[str]
            List of scatterer type names.
        """
        return list(self._types.keys())

    def get_dropdown_options(self) -> List[Dict[str, str]]:
        """
        Get dropdown options for scatterer types (alias for get_available_types).

        Returns
        -------
        List[Dict[str, str]]
            List of dictionaries with 'label' and 'value' keys for Dash dropdown.
        """
        return self.get_available_types()

    def get_input_definitions(self, scatterer_name: str) -> Dict[str, Dict[str, Any]]:
        """
        Get input definitions for a specific scatterer type.

        Parameters
        ----------
        scatterer_name : str
            Name of the scatterer type.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            Input definitions for the scatterer type.
        """
        return self.get_type(scatterer_name).get_input_definitions()

    def validate_inputs(self, scatterer_name: str, input_values: Dict[str, str]) -> tuple[bool, str]:
        """
        Validate inputs for a specific scatterer type.

        Parameters
        ----------
        scatterer_name : str
            Name of the scatterer type.
        input_values : Dict[str, str]
            Dictionary of input values to validate.

        Returns
        -------
        tuple[bool, str]
            (is_valid, error_message) where is_valid is bool and error_message is str.
        """
        return self.get_type(scatterer_name).validate_inputs(input_values)

    def get_parameter_mapping(self, scatterer_name: str) -> Dict[str, str]:
        """
        Get parameter mapping for a specific scatterer type.

        Parameters
        ----------
        scatterer_name : str
            Name of the scatterer type.

        Returns
        -------
        Dict[str, str]
            Parameter mapping for PyMieSim classes.
        """
        return self.get_type(scatterer_name).get_parameter_mapping()

    def get_scatterer_info(self, scatterer_name: str) -> Dict[str, Any]:
        """
        Get descriptive information about a scatterer type.

        Parameters
        ----------
        scatterer_name : str
            Name of the scatterer type.

        Returns
        -------
        Dict[str, Any]
            Information dictionary with description and use cases.
        """
        scatterer_type = self.get_type(scatterer_name)
        return {
            'name': scatterer_type.name,
            'label': scatterer_type.label,
            'description': scatterer_type.description,
            'use_cases': scatterer_type.use_cases
        }
// Add a default manager instance for easy imports
scatterer_manager = ScattererTypeManager()
