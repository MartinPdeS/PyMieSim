# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/experiment/dataframe\_subclass.py     |       67 |        8 |       16 |        3 |     87% |74->77, 87, 200, 210-214, 222-226 |
| PyMieSim/experiment/detector/base.py           |       38 |        0 |        8 |        0 |    100% |           |
| PyMieSim/experiment/detector/coherent\_mode.py |       20 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/detector/photodiode.py     |       23 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       42 |        2 |       12 |        2 |     93% |   89, 106 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       26 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/cylinder.py      |       23 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/sphere.py        |       23 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/setup.py                   |      104 |        3 |       26 |        4 |     95% |299-302, 363->366, 369->373, 378 |
| PyMieSim/experiment/source/base.py             |       35 |        2 |        8 |        2 |     91% |    47, 65 |
| PyMieSim/experiment/source/gaussian.py         |       21 |        1 |        2 |        1 |     91% |        49 |
| PyMieSim/experiment/source/planewave.py        |       27 |        1 |        4 |        2 |     90% |43->46, 60 |
| PyMieSim/experiment/utils.py                   |       32 |        1 |       14 |        2 |     93% |44->50, 58 |
| PyMieSim/mesh.py                               |       50 |        4 |        0 |        0 |     92% |100, 112, 124, 136 |
| PyMieSim/polarization.py                       |       45 |        1 |        8 |        0 |     98% |        87 |
| PyMieSim/single/detector/base.py               |       39 |        1 |        0 |        0 |     97% |        33 |
| PyMieSim/single/detector/coherent.py           |       27 |        2 |        4 |        2 |     87% |    72, 79 |
| PyMieSim/single/detector/uncoherent.py         |       22 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/representations/base.py        |       41 |        3 |        2 |        1 |     91% | 51, 55-56 |
| PyMieSim/single/representations/far\_fields.py |       33 |        0 |        6 |        0 |    100% |           |
| PyMieSim/single/representations/footprint.py   |       45 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/representations/near\_field.py |       74 |       53 |       14 |        0 |     24% |66-70, 75-84, 89-91, 112-121, 134-166, 189-237 |
| PyMieSim/single/representations/s1s2.py        |       21 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/representations/spf.py         |       28 |        0 |        2 |        0 |    100% |           |
| PyMieSim/single/representations/stokes.py      |       31 |        0 |        2 |        0 |    100% |           |
| PyMieSim/single/scatterer/base.py              |      127 |       23 |       12 |        1 |     77% |500-537, 573, 671 |
| PyMieSim/single/scatterer/core\_shell.py       |       33 |        2 |        0 |        0 |     94% |   119-124 |
| PyMieSim/single/scatterer/cylinder.py          |       37 |        6 |        0 |        0 |     84% |73, 77, 81, 85, 107-116 |
| PyMieSim/single/scatterer/sphere.py            |       30 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/source/base.py                 |       10 |        0 |        2 |        0 |    100% |           |
| PyMieSim/single/source/gaussian.py             |       36 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/source/planewave.py            |       24 |        0 |        0 |        0 |    100% |           |
| PyMieSim/special\_functions.py                 |       18 |        0 |        2 |        0 |    100% |           |
| PyMieSim/utils.py                              |        2 |        0 |        0 |        0 |    100% |           |
|                                      **TOTAL** | **1254** |  **113** |  **144** |   **20** | **89%** |           |


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