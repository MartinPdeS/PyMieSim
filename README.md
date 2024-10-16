# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/\_\_main\_\_.py                       |        9 |        9 |        2 |        0 |      0% |      2-15 |
| PyMieSim/experiment/detector/base.py           |       57 |        9 |       16 |        5 |     75% |80-84, 103, 106, 126, 129, 148->151 |
| PyMieSim/experiment/detector/coherent\_mode.py |       11 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/detector/photodiode.py     |       11 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       46 |        8 |       14 |        5 |     75% |48-50, 60, 89-94, 109->107, 114 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       21 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/cylinder.py      |       18 |        6 |        0 |        0 |     67% |     44-55 |
| PyMieSim/experiment/scatterer/sphere.py        |       18 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/setup.py                   |       95 |       10 |       26 |        6 |     85% |106, 133, 210, 238-239, 264-269, 304->309, 314 |
| PyMieSim/experiment/source/base.py             |       37 |        5 |       12 |        2 |     78% |32->35, 47-53, 67->64 |
| PyMieSim/experiment/source/gaussian.py         |       23 |        2 |        6 |        3 |     83% |36, 39, 46->49 |
| PyMieSim/experiment/source/planewave.py        |       15 |        4 |        2 |        0 |     65% | 31-34, 44 |
| PyMieSim/mesh.py                               |       94 |        8 |        8 |        2 |     90% |116, 128, 140, 152, 194->exit, 213->exit, 274-278 |
| PyMieSim/polarization.py                       |       45 |        5 |        8 |        1 |     85% |70, 86, 89-92 |
| PyMieSim/single/detector/base.py               |       77 |        4 |       10 |        4 |     91% |91, 102, 105, 229 |
| PyMieSim/single/detector/coherent.py           |       41 |        2 |       10 |        3 |     90% |47, 52, 64->68 |
| PyMieSim/single/detector/uncoherent.py         |       28 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/representations.py             |      168 |        3 |       12 |        1 |     98% | 57, 61-62 |
| PyMieSim/single/scatterer/base.py              |      103 |        4 |       10 |        4 |     93% |47, 50, 60, 451 |
| PyMieSim/single/scatterer/core\_shell.py       |       26 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/scatterer/cylinder.py          |       39 |        4 |        0 |        0 |     90% |159, 163, 167, 171 |
| PyMieSim/single/scatterer/sphere.py            |       26 |        2 |        0 |        0 |     92% |  136, 164 |
| PyMieSim/single/source/base.py                 |       33 |        7 |       16 |        7 |     71% |22, 25, 34, 37, 46, 49, 63 |
| PyMieSim/single/source/gaussian.py             |       35 |        1 |        2 |        1 |     95% |        41 |
| PyMieSim/single/source/planewave.py            |       26 |        1 |        2 |        1 |     93% |        37 |
| PyMieSim/special\_functions.py                 |       20 |        0 |        2 |        0 |    100% |           |
| PyMieSim/units.py                              |       19 |        0 |        4 |        0 |    100% |           |
|                                      **TOTAL** | **1141** |   **94** |  **162** |   **45** | **88%** |           |


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