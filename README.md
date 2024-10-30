# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/\_\_main\_\_.py                       |        9 |        9 |        2 |        0 |      0% |      2-15 |
| PyMieSim/experiment/detector/base.py           |       58 |       27 |       16 |        0 |     42% |82-86, 104-110, 127-130, 147-153, 162-180, 189-193 |
| PyMieSim/experiment/detector/coherent\_mode.py |       11 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/detector/photodiode.py     |       11 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       49 |       24 |       14 |        0 |     40% |42-50, 57-62, 79-89, 102-110 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       21 |        7 |        0 |        0 |     67% |     56-70 |
| PyMieSim/experiment/scatterer/cylinder.py      |       18 |        6 |        0 |        0 |     67% |     44-56 |
| PyMieSim/experiment/scatterer/sphere.py        |       18 |        6 |        0 |        0 |     67% |     48-60 |
| PyMieSim/experiment/setup.py                   |       93 |       66 |       24 |        0 |     23% |41-42, 48-50, 55-62, 65-69, 86-90, 113-129, 146, 168-187, 203-210, 228-237, 257-262, 279-316 |
| PyMieSim/experiment/source/base.py             |       58 |       34 |       24 |        0 |     29% |31-42, 49-52, 59-65, 72-78, 85-92, 103-114 |
| PyMieSim/experiment/source/gaussian.py         |        9 |        1 |        0 |        0 |     89% |        36 |
| PyMieSim/experiment/source/planewave.py        |       15 |        4 |        2 |        0 |     65% | 32-35, 45 |
| PyMieSim/mesh.py                               |       94 |       94 |        8 |        0 |      0% |     4-467 |
| PyMieSim/polarization.py                       |       45 |       22 |        8 |        0 |     43% |65, 68-71, 74-84, 87, 90-93, 114, 135, 169-174 |
| PyMieSim/single/detector/base.py               |       77 |       77 |       10 |        0 |      0% |     4-394 |
| PyMieSim/single/detector/coherent.py           |       39 |       39 |       10 |        0 |      0% |     4-115 |
| PyMieSim/single/detector/uncoherent.py         |       26 |       26 |        0 |        0 |      0% |     4-104 |
| PyMieSim/single/representations.py             |      168 |      168 |       12 |        0 |      0% |     4-734 |
| PyMieSim/single/scatterer/base.py              |      101 |      101 |        8 |        0 |      0% |     4-451 |
| PyMieSim/single/scatterer/core\_shell.py       |       25 |       25 |        0 |        0 |      0% |     4-101 |
| PyMieSim/single/scatterer/cylinder.py          |       38 |       38 |        0 |        0 |      0% |     4-169 |
| PyMieSim/single/scatterer/sphere.py            |       26 |       26 |        0 |        0 |      0% |     4-163 |
| PyMieSim/single/source/base.py                 |       27 |       27 |       10 |        0 |      0% |      4-58 |
| PyMieSim/single/source/gaussian.py             |       35 |       35 |        2 |        0 |      0% |     4-118 |
| PyMieSim/single/source/planewave.py            |       26 |       26 |        2 |        0 |      0% |     4-103 |
| PyMieSim/special\_functions.py                 |       20 |       20 |        2 |        0 |      0% |      4-99 |
| PyMieSim/units.py                              |       19 |        0 |        4 |        0 |    100% |           |
|                                      **TOTAL** | **1136** |  **908** |  **158** |    **0** | **18%** |           |


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