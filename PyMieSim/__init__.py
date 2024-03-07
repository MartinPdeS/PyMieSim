import numpy
from pathlib import Path
from PyMieSim.tools.directories import lp_mode_path


def load_single_lp_mode(mode_number: str, structure_type: str = 'unstructured', sampling: int = 500) -> numpy.ndarray:
    """
    Load a single LP mode from a file based on mode number, structure type, and sampling rate.

    Args:
        mode_number (str): The mode number to load.
        structure_type (str, optional): The type of structure (e.g., 'unstructured'). Defaults to 'unstructured'.
        sampling (int, optional): The sampling rate. Defaults to 500.

    Returns:
        numpy.ndarray: The loaded LP mode field as a numpy array.
    """
    structure_type = structure_type.lower()
    mode_number = mode_number.split(':')[0]  # Keep only the mode number before ':' if present
    lp_file_directory = lp_mode_path.joinpath(mode_number)

    if not lp_file_directory.exists():
        raise FileNotFoundError(f"This LP mode: {mode_number} is not available.\nAvailable are {list_available_modes()}")

    lp_sampling_file_directory = lp_file_directory.joinpath(f"{structure_type}_{sampling}.npy")

    if not lp_sampling_file_directory.exists():
        raise FileNotFoundError(f"This sampling: {sampling} is not available.\nAvailable are {list_available_sampling(lp_file_directory, structure_type)}")

    return numpy.load(lp_sampling_file_directory).astype(complex).squeeze()


def list_available_modes() -> list:
    """List all available LP modes."""
    return [path.name for path in lp_mode_path.iterdir() if path.is_dir() and path.name.startswith('LP')]


def list_available_sampling(directory: Path, structure_type: str) -> list:
    """List all available samplings for a given structure type in the specified directory."""
    return [path.stem for path in directory.iterdir() if path.is_file() and path.name.startswith(structure_type)]


def load_lp_mode(mode_numbers: list, structure_type: str = 'unstructured', sampling: int = 500) -> numpy.ndarray:
    """
    Load multiple LP modes specified by their mode numbers.

    Args:
        mode_numbers (list): A list of mode numbers to load.
        structure_type (str, optional): The type of structure. Defaults to 'unstructured'.
        sampling (int, optional): The sampling rate. Defaults to 500.

    Returns:
        numpy.ndarray: An array of loaded LP mode fields.
    """
    fields_array = numpy.asarray([
        load_single_lp_mode(mode, structure_type, sampling) for mode in mode_numbers
    ])
    return fields_array
