from PyMieSim.Data._Material.utils import LoadOnlineSave
from PyMieSim import Material

LoadOnlineSave(filename='Silver', url='https://refractiveindex.info/data_csv.php?datafile=data/main/Ag/Johnson.yml')

Mat = Material('Silver')

Mat.Plot()
