# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/experiment/dataframe\_subclass.py     |       67 |        8 |       16 |        3 |     87% |70->73, 83, 192, 202-206, 215-219 |
| PyMieSim/experiment/detector/base.py           |       43 |        0 |       16 |        1 |     98% |  125->128 |
| PyMieSim/experiment/detector/coherent\_mode.py |       19 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/detector/photodiode.py     |       22 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       44 |        2 |       14 |        2 |     93% |    78, 93 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       25 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/cylinder.py      |       22 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/sphere.py        |       22 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/setup.py                   |      104 |        3 |       26 |        4 |     95% |275-276, 331->334, 339->343, 348 |
| PyMieSim/experiment/source/base.py             |       44 |        2 |       22 |        5 |     89% |34->37, 47->50, 61, 63->exit, 78 |
| PyMieSim/experiment/source/gaussian.py         |       20 |        1 |        2 |        1 |     91% |        44 |
| PyMieSim/experiment/source/planewave.py        |       26 |        1 |        4 |        2 |     90% |38->41, 55 |
| PyMieSim/experiment/utils.py                   |       34 |        1 |       14 |        2 |     94% |52->56, 64 |
| PyMieSim/mesh.py                               |       61 |        8 |        4 |        1 |     83% |68, 76-79, 134, 146, 158, 170 |
| PyMieSim/polarization.py                       |       45 |        1 |        8 |        0 |     98% |        87 |
| PyMieSim/single/detector/base.py               |       46 |        4 |        2 |        0 |     88% | 35, 42-45 |
| PyMieSim/single/detector/coherent.py           |       23 |        2 |        4 |        2 |     85% |    63, 68 |
| PyMieSim/single/detector/uncoherent.py         |       18 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/representations.py             |      168 |        3 |       12 |        1 |     98% | 58, 62-63 |
| PyMieSim/single/scatterer/base.py              |       93 |        1 |        4 |        1 |     98% |       458 |
| PyMieSim/single/scatterer/core\_shell.py       |       36 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/scatterer/cylinder.py          |       38 |        4 |        0 |        0 |     89% |72, 76, 80, 84 |
| PyMieSim/single/scatterer/sphere.py            |       31 |        2 |        0 |        0 |     94% |    99-107 |
| PyMieSim/single/source/base.py                 |        3 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/source/gaussian.py             |       32 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/source/planewave.py            |       21 |        0 |        0 |        0 |    100% |           |
| PyMieSim/special\_functions.py                 |       18 |        0 |        2 |        0 |    100% |           |
| PyMieSim/units.py                              |       48 |        2 |       20 |        2 |     94% |   87, 128 |
|                                      **TOTAL** | **1173** |   **45** |  **170** |   **27** | **94%** |           |


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/MartinPdeS/PyMieSim/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/MartinPdeS/PyMieSim/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2FMartinPdeS%2FPyMieSim%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.