# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/experiment/dataframe\_subclass.py     |       67 |       15 |       16 |        3 |     76% |68->73, 81, 158-165, 200, 208-212, 221-225 |
| PyMieSim/experiment/detector/base.py           |       43 |        0 |       16 |        1 |     98% |  125->128 |
| PyMieSim/experiment/detector/coherent\_mode.py |       19 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/detector/photodiode.py     |       22 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       40 |        2 |       14 |        2 |     93% |    72, 87 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       26 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/cylinder.py      |       23 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/sphere.py        |       23 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/setup.py                   |      107 |        3 |       32 |        3 |     96% |269-270, 333->337, 342 |
| PyMieSim/experiment/source/base.py             |       52 |        3 |       24 |        6 |     88% |20, 48->51, 61->64, 75, 77->exit, 92 |
| PyMieSim/experiment/source/gaussian.py         |       15 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/source/planewave.py        |       21 |        0 |        2 |        1 |     96% |    37->40 |
| PyMieSim/experiment/utils.py                   |       34 |        1 |       14 |        2 |     94% |52->56, 64 |
| PyMieSim/mesh.py                               |       88 |        4 |        8 |        2 |     94% |114, 126, 138, 150, 192->exit, 211->exit |
| PyMieSim/polarization.py                       |       45 |        1 |        8 |        0 |     98% |        87 |
| PyMieSim/single/detector/base.py               |       86 |        6 |       14 |        1 |     89% |99-105, 243 |
| PyMieSim/single/detector/coherent.py           |       39 |        2 |       10 |        3 |     90% |49, 54, 66->70 |
| PyMieSim/single/detector/uncoherent.py         |       26 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/representations.py             |      168 |        3 |       12 |        1 |     98% | 58, 62-63 |
| PyMieSim/single/scatterer/base.py              |      101 |        1 |        8 |        1 |     98% |       451 |
| PyMieSim/single/scatterer/core\_shell.py       |       25 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/scatterer/cylinder.py          |       38 |        4 |        0 |        0 |     89% |157, 161, 165, 169 |
| PyMieSim/single/scatterer/sphere.py            |       26 |        2 |        0 |        0 |     92% |  135, 163 |
| PyMieSim/single/source/base.py                 |       27 |        0 |       10 |        0 |    100% |           |
| PyMieSim/single/source/gaussian.py             |       35 |        1 |        2 |        1 |     95% |        41 |
| PyMieSim/single/source/planewave.py            |       26 |        1 |        2 |        1 |     93% |        37 |
| PyMieSim/special\_functions.py                 |       20 |        0 |        2 |        0 |    100% |           |
| PyMieSim/units.py                              |       19 |        0 |        4 |        0 |    100% |           |
|                                      **TOTAL** | **1261** |   **49** |  **198** |   **28** | **94%** |           |


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