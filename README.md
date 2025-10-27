# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/experiment/dataframe\_subclass.py     |       65 |       12 |       10 |        1 |     80% |150-156, 195, 221-225, 244-248 |
| PyMieSim/experiment/detector/base.py           |       39 |        0 |        8 |        0 |    100% |           |
| PyMieSim/experiment/detector/coherent\_mode.py |       14 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/detector/photodiode.py     |       18 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       38 |        2 |       12 |        2 |     92% |    77, 94 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       20 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/cylinder.py      |       19 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/sphere.py        |       19 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/setup.py                   |      102 |        3 |       26 |        4 |     95% |300-303, 364->367, 370->374, 379 |
| PyMieSim/experiment/source/base.py             |       40 |        2 |       10 |        3 |     90% |47, 56->59, 73 |
| PyMieSim/experiment/source/gaussian.py         |       15 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/source/planewave.py        |       15 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/utils.py                   |       32 |        1 |       14 |        2 |     93% |44->50, 58 |
| PyMieSim/single/detector/base.py               |       39 |        1 |        0 |        0 |     97% |        33 |
| PyMieSim/single/detector/coherent.py           |       23 |        2 |        4 |        2 |     85% |    72, 79 |
| PyMieSim/single/detector/uncoherent.py         |       19 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/mesh.py                        |       46 |        4 |        0 |        0 |     91% |100, 112, 124, 136 |
| PyMieSim/single/polarization.py                |       43 |        1 |        8 |        0 |     98% |        87 |
| PyMieSim/single/representations/base.py        |       38 |        3 |        2 |        1 |     90% | 50, 54-55 |
| PyMieSim/single/representations/far\_fields.py |       32 |        0 |        6 |        0 |    100% |           |
| PyMieSim/single/representations/footprint.py   |       42 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/representations/near\_field.py |       68 |       53 |       14 |        0 |     18% |66-70, 75-84, 89-91, 112-121, 134-166, 189-237 |
| PyMieSim/single/representations/s1s2.py        |       18 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/representations/spf.py         |       27 |        0 |        2 |        0 |    100% |           |
| PyMieSim/single/representations/stokes.py      |       30 |        0 |        2 |        0 |    100% |           |
| PyMieSim/single/scatterer/base.py              |      127 |       23 |       12 |        1 |     77% |502-539, 575, 673 |
| PyMieSim/single/scatterer/core\_shell.py       |       27 |        2 |        0 |        0 |     93% |   119-124 |
| PyMieSim/single/scatterer/cylinder.py          |       33 |        6 |        0 |        0 |     82% |73, 77, 81, 85, 107-116 |
| PyMieSim/single/scatterer/sphere.py            |       25 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/source/base.py                 |       10 |        0 |        2 |        0 |    100% |           |
| PyMieSim/single/source/gaussian.py             |       32 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/source/planewave.py            |       21 |        0 |        0 |        0 |    100% |           |
| PyMieSim/utils.py                              |       20 |        0 |        2 |        0 |    100% |           |
|                                      **TOTAL** | **1156** |  **115** |  **134** |   **16** | **88%** |           |


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