# GitHub Copilot Instructions for PyMieSim

## Project Overview

PyMieSim is an open-source Python package for fast and flexible Mie scattering simulations, targeting both single-scatterer studies and large parametric experiments. The package combines high-performance C++ computational cores with an intuitive Python API for optical physics research and education.

### Core Domain Knowledge

- **Mie Scattering Theory**: Classical electromagnetic scattering by spherical, cylindrical, and core-shell particles
- **Optical Physics**: Light-matter interaction, refractive indices, scattering cross-sections, efficiencies (Qsca, Qabs, Qext)
- **Computational Electromagnetics**: Special functions (Bessel, spherical harmonics), series expansions, complex analysis
- **Scientific Computing**: Parametric studies, data analysis with pandas, visualization with matplotlib

### Key Concepts

- **Scatterers**: Sphere, Cylinder, CoreShell geometries with physical properties (diameter, refractive index)
- **Sources**: Gaussian beams, plane waves with wavelength, polarization, optical power, numerical aperture
- **Detectors**: Photodiodes, coherent mode detectors for capturing scattered light
- **Experiments**: Setup class orchestrating scatterer-source-detector combinations for parametric studies
- **Units**: Physical quantities with proper dimensional analysis using pint/PyOptik integration
- **Cross-sections**: Geometric, scattering, absorption, extinction cross-sections and efficiencies

## Architecture & Components

### Core Structure
```
PyMieSim/
├── cpp/                    # High-performance C++ computational cores
│   ├── scatterer/         # Sphere, cylinder, core-shell implementations
│   ├── source/            # Gaussian, plane wave light sources
│   ├── detector/          # Photodiode, coherent mode detectors
│   ├── bessel/            # Special functions (Bessel, spherical)
│   ├── sets/              # Parameter set management
│   ├── experiment/        # Experiment orchestration
│   └── utils/             # Utilities, numpy interface
├── experiment/            # Python experiment framework
│   ├── scatterer.py       # Python scatterer interfaces
│   ├── source.py          # Python source interfaces
│   ├── detector.py        # Python detector interfaces
│   ├── setup.py           # Main experiment orchestration
│   └── dataframe_subclass.py # Custom pandas DataFrame
├── single/                # Single particle calculations
├── gui/                   # Graphical user interface
├── units.py              # Physical units and dimensional analysis
├── polarization.py       # Polarization handling
└── validation_data/      # Reference data for validation
```

### Technology Stack
- **Core Engine**: Modern C++ (C++17) with CMake build system
- **Python Bindings**: pybind11 for seamless C++/Python integration
- **Scientific Stack**: NumPy, pandas, matplotlib for data handling and visualization
- **Units**: pint/PyOptik for physical quantities and dimensional analysis
- **Documentation**: Sphinx with gallery examples, Doxygen for C++
- **Testing**: pytest with parametric tests and validation against reference data
- **Building**: scikit-build-core for hybrid C++/Python packages

## Coding Standards & Conventions

### Python Code
- **Style**: PEP 8 compliance with meaningful variable names reflecting physical quantities
- **Docstrings**: NumPy-style with physics-aware parameter descriptions
- **Type Hints**: Comprehensive typing for scientific code clarity
- **Units**: Mandatory dimensional analysis with pint quantities
- **Error Handling**: Informative exceptions with physics context

#### Physics-Aware Naming
```python
# Preferred - reflects physical meaning
wavelength = 532 * nanometer
scattering_efficiency = sphere.compute_Qsca()
refractive_index = 1.33 + 0.01j * RIU

# Avoid - generic or unclear naming
x = 532
result = sphere.calc()
n = 1.33 + 0.01j
```

#### Docstring Standards
```python
def compute_scattering_efficiency(self, wavelength: Quantity) -> float:
    """
    Compute the scattering efficiency (Qsca) for the particle.

    The scattering efficiency is the ratio of the scattering cross-section
    to the geometric cross-section of the particle, providing a measure
    of how effectively the particle scatters incident light.

    Parameters
    ----------
    wavelength : Quantity
        Incident wavelength in length units (typically nanometers).

    Returns
    -------
    float
        Dimensionless scattering efficiency Qsca.

    Notes
    -----
    For spherical particles, this implements the exact Mie solution.
    The efficiency depends strongly on the size parameter x = πd/λ
    and the relative refractive index m = n_particle/n_medium.
    """
```

### C++ Code
- **Style**: Modern C++17 with clear class hierarchies
- **Documentation**: Doxygen comments with mathematical context
- **Memory Management**: RAII principles, smart pointers where appropriate
- **Performance**: Vectorized operations, OpenMP parallelization
- **Interface**: Clean pybind11 bindings with Python-friendly signatures

#### C++ Documentation Standards
```cpp
/**
 * @brief Computes Mie scattering coefficients for a spherical particle.
 *
 * This method calculates the complex scattering coefficients a_n and b_n
 * using the exact Mie theory formulation. The coefficients are computed
 * up to the specified maximum order using recurrence relations for
 * spherical Bessel and Hankel functions.
 *
 * @param max_order Maximum multipole order for coefficient computation.
 *                   If 0, automatically determined from size parameter.
 *
 * @note The computation uses the Wiscombe algorithm for numerical stability
 *       in the calculation of spherical Bessel functions for large arguments.
 *
 * @see Bohren & Huffman, "Absorption and Scattering of Light by Small Particles"
 */
void compute_an_bn(size_t max_order = 0);
```

## Testing Philosophy

### Test Categories
1. **Unit Tests**: Individual component validation
2. **Integration Tests**: Multi-component workflow validation
3. **Physics Validation**: Comparison against analytical solutions and literature
4. **Performance Tests**: Benchmarking computational efficiency
5. **Cross-platform Tests**: Windows, macOS, Linux compatibility

### Test Standards
```python
import pytest
import numpy as np
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import PlaneWave
from PyMieSim.units import nanometer, RIU

class TestMieScattering:
    """Test suite for Mie scattering calculations."""

    @pytest.fixture
    def rayleigh_sphere(self):
        """Small sphere in Rayleigh regime for analytical validation."""
        return Sphere(
            diameter=10 * nanometer,
            property=1.5 * RIU,
            medium_property=1.0 * RIU
        )

    @pytest.mark.parametrize("size_param,expected_regime", [
        (0.1, "rayleigh"),
        (1.0, "resonance"),
        (10.0, "geometric")
    ])
    def test_scattering_regimes(self, size_param, expected_regime):
        """Test scattering behavior across different size regimes."""
        # Test implementation
        pass

    def test_rayleigh_limit_analytical(self, rayleigh_sphere):
        """Validate against analytical Rayleigh scattering formula."""
        # Compare computed vs analytical results
        pass
```

## Documentation Standards

### Sphinx Gallery Examples
- **Cell Structure**: Use `# %%` for code cells in example scripts
- **Physics Focus**: Emphasize physical interpretation of results
- **Visualization**: Professional plots with proper units and labels
- **Self-Contained**: Each example should run independently
- **Progressive Complexity**: Simple to advanced examples

#### Example Structure
```python
"""
Sphere Scattering Efficiency vs Size Parameter
==============================================

This example demonstrates the characteristic resonance behavior
of Mie scattering efficiency as a function of particle size.
"""

# %%
# Import the necessary packages
import numpy as np
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import PlaneWave
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, RIU

# %%
# Define the experimental parameters
wavelength = 532 * nanometer
source = PlaneWave(
    wavelength=wavelength,
    polarization=0 * degree
)

# %%
# Create spheres with varying diameter
diameters = np.linspace(10, 1000, 200) * nanometer
sphere = Sphere(
    diameter=diameters,
    property=1.5 * RIU,  # Polystyrene-like
    medium_property=1.0 * RIU,  # Air
    source=source
)

# %%
# Execute the experiment
experiment = Setup(scatterer=sphere, source=source)
dataframe = experiment.get('Qsca')

# %%
# Visualize the results
dataframe.plot_data(
    x='scatterer:diameter',
    title='Mie Scattering Efficiency vs Particle Size'
)
```

### API Documentation
- **Physics Context**: Explain the underlying physics principles
- **Mathematical Notation**: Use proper mathematical symbols and equations
- **Units**: Clearly specify expected units and dimensional analysis
- **Examples**: Provide realistic physics examples
- **Cross-references**: Link to relevant theory and literature

## Development Workflow

### Feature Development
1. **Physics Research**: Understand the theoretical background
2. **C++ Implementation**: High-performance core algorithms
3. **Python Bindings**: User-friendly interface with proper units
4. **Validation**: Compare against analytical/reference solutions
5. **Documentation**: Examples with physical interpretation
6. **Testing**: Comprehensive test coverage including edge cases

### Code Review Focus
- **Physical Correctness**: Verify equations and implementations
- **Numerical Stability**: Check for precision issues, special cases
- **Performance**: Ensure computational efficiency for large datasets
- **API Design**: Intuitive interface for physicists and engineers
- **Documentation**: Clear explanations with proper physics context

### Build & Integration
- **CMake**: Cross-platform C++ builds with dependency management
- **pybind11**: Seamless Python/C++ integration with proper error handling
- **CI/CD**: Automated testing across platforms and Python versions
- **Packaging**: Proper wheel building for scientific Python ecosystem

## Domain-Specific Guidelines

### Physics Accuracy
- Always validate against known analytical solutions (Rayleigh, geometric optics limits)
- Use proper electromagnetic field formulations and boundary conditions
- Implement numerical stability measures for extreme parameter ranges
- Reference authoritative sources (Bohren & Huffman, Van de Hulst, etc.)

### Performance Considerations
- Leverage vectorization for large parameter sweeps
- Use OpenMP for computationally intensive calculations
- Optimize special function evaluations (Bessel functions, series convergence)
- Memory-efficient handling of large result datasets

### User Experience
- Provide sensible defaults for physical parameters
- Clear error messages with physics context
- Comprehensive units support with automatic conversions
- Rich visualization capabilities with proper scientific plotting standards

### Validation & Testing
- Compare against established reference data (e.g., Wiscombe's Mie scattering test cases)
- Test limiting cases (small particles, large particles, extreme refractive indices)
- Cross-validate with other established Mie scattering codes
- Include performance benchmarks for computational efficiency

## Common Patterns & Best Practices

### Experiment Setup Pattern
```python
# Standard experimental workflow
source = PlaneWave(wavelength=wavelength, polarization=polarization)
scatterer = Sphere(diameter=diameter, property=material, source=source)
experiment = Setup(scatterer=scatterer, source=source)
results = experiment.get('Qsca')
```

### Unit Handling
```python
# Always use physical units
from PyMieSim.units import nanometer, micrometer, RIU, degree

wavelength = 632.8 * nanometer  # He-Ne laser
diameter = 1.0 * micrometer    # Typical aerosol particle
n_water = 1.33 * RIU           # Water refractive index
```

### Error Handling
```python
# Physics-aware error messages
if size_parameter > 1000:
    warnings.warn(
        f"Size parameter {size_parameter:.1f} is very large. "
        "Consider using geometric optics approximation for better performance."
    )
```

This document serves as a comprehensive guide for developing high-quality, physics-accurate code within the PyMieSim ecosystem. When contributing, always prioritize physical correctness, numerical stability, and user experience appropriate for the scientific computing community.
