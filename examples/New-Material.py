from PyMieSim._Material.utils import LoadOnlineSave
from PyMieSim import Material

LoadOnlineSave(filename='BK7', url='https://refractiveindex.info/data_csv.php?datafile=data/glass/schott/N-BK7.yml')

Mat = Material('BK7')

Mat.Plot()
