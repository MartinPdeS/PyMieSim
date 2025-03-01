from PyMieSim.experiment.source.gaussian import Gaussian
from PyMieSim.experiment.scatterer.sphere import Sphere
from PyMieSim.experiment.detector.photodiode import Photodiode
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, milliwatt, AU, RIU
import numpy


length_units = nanometer
power_units = milliwatt
angle_units = degree


def parse_string_to_array_or_float(input_str):
    """
    Parse a string to return either a numpy array or a float.

    If the string is in the format 'start:end:count', return np.linspace(start, end, count).
    If the string is a comma-separated list, return a numpy array.
    If the string is a single numeric value, return it as a float.
    """
    try:
        if ":" in input_str:
            start, end, count = map(float, input_str.split(":"))
            return numpy.linspace(start, end, int(count))
        elif "," in input_str:
            return numpy.array(list(map(float, input_str.split(","))))
        else:
            return float(input_str)
    except ValueError:
        raise ValueError("Invalid input string format. Expected 'start:end:count', a comma-separated list, or a single numeric value.")


def interface(source_kwargs: dict, scatterer_kwargs: dict, detector_kwargs: dict, measure: str, **kwargs):
    source = Gaussian(
        wavelength=parse_string_to_array_or_float(source_kwargs['wavelength']) * length_units,
        polarization=parse_string_to_array_or_float(source_kwargs['polarization']) * angle_units,
        NA=parse_string_to_array_or_float(source_kwargs['NA']) * AU,
        optical_power=parse_string_to_array_or_float(source_kwargs['optical_power']) * milliwatt
    )

    scatterer = Sphere(
        diameter=parse_string_to_array_or_float(scatterer_kwargs['diameter']) * length_units,
        property=parse_string_to_array_or_float(scatterer_kwargs['property']) * RIU,
        medium_property=parse_string_to_array_or_float(scatterer_kwargs['medium_property']) * RIU,
        source=source
    )

    detector = Photodiode(
        NA=parse_string_to_array_or_float(detector_kwargs['NA']) * AU,
        gamma_offset=parse_string_to_array_or_float(detector_kwargs['gamma_offset']) * angle_units,
        phi_offset=parse_string_to_array_or_float(detector_kwargs['phi_offset']) * angle_units,
        sampling=parse_string_to_array_or_float(detector_kwargs['sampling']) * AU,
    )

    setup = Setup(source=source, scatterer=scatterer, detector=detector)

    dataframe = setup.get(measure, **kwargs)

    return dataframe
