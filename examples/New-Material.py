matplotlib=True
mlab=False

def run():
    from PyMieSim.Data._Material.utils import LoadOnlineSave
    from PyMieSim                      import Material

    LoadOnlineSave(filename='BK7', url='https://refractiveindex.info/data_csv.php?datafile=data/glass/schott/N-BK7.yml')

    Mat = Material('BK7')

    Mat.Plot()


if __name__ == '__main__':
    run()
