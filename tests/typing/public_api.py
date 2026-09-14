"""Static regression checks for public overloads; checked by mypy, not executed."""
from typing_extensions import assert_type
from pint import Quantity

from PyMieSim.single.api import Simulation
from PyMieSim.experiment.setup import Setup as Experiment
from PyMieSim.results import SimulationResult, SimulationResults
from PyMieSim.labeled_array import LabeledArray
from PyMieSim.measures import Measure
from PyMieSim import ParticleSizeDistribution


def check_result_types(simulation: Simulation, experiment: Experiment) -> None:
    assert_type(simulation.run(Measure.QSCA), Quantity)
    assert_type(simulation.get('Qsca', 'Csca'), dict[str, Quantity])
    assert_type(simulation.run('Qsca', as_result=True), SimulationResult)
    assert_type(simulation.get('Qsca', 'Csca', as_result=True), SimulationResults)
    assert_type(experiment.run(Measure.QSCA), LabeledArray)
    assert_type(experiment.get('Qsca', 'Csca'), LabeledArray)
    assert_type(experiment.run('Qsca', as_result=True), SimulationResult)
    assert_type(experiment.get('Qsca', 'Csca', as_result=True), SimulationResults)
    assert_type(experiment.run('Csca', as_result=True).to('nanometer ** 2').as_labeled_array(), LabeledArray)


def check_distribution_types(experiment: Experiment, distribution: ParticleSizeDistribution) -> None:
    assert_type(experiment.average_size_distribution(distribution, "Csca"), LabeledArray)
    assert_type(experiment.average_size_distribution(distribution, "Csca", as_result=True), SimulationResult)
    assert_type(experiment.average_size_distribution(distribution, "Csca", "g", as_result=True), SimulationResults)
