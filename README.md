# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/\_\_main\_\_.py                       |        9 |        9 |        2 |        0 |      0% |      2-15 |
| PyMieSim/experiment/detector/base.py           |       58 |        9 |       16 |        5 |     76% |82-86, 105, 108, 128, 148, 150->153 |
| PyMieSim/experiment/detector/coherent\_mode.py |       11 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/detector/photodiode.py     |       11 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       49 |        9 |       14 |        5 |     75% |48-50, 60, 84-89, 105, 110 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       21 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/cylinder.py      |       18 |        6 |        0 |        0 |     67% |     44-56 |
| PyMieSim/experiment/scatterer/sphere.py        |       18 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/setup.py                   |       93 |        8 |       24 |        4 |     88% |127, 231-232, 257-262, 295->299, 304 |
| PyMieSim/experiment/source/base.py             |       58 |        7 |       24 |       10 |     79% |34, 50, 60, 62->65, 73, 75->78, 86, 89, 91->exit, 107 |
| PyMieSim/experiment/source/gaussian.py         |        9 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/source/planewave.py        |       15 |        4 |        2 |        0 |     65% | 32-35, 45 |
| PyMieSim/mesh.py                               |       94 |       50 |        8 |        0 |     43% |89-105, 117, 129, 141, 153, 165, 177, 193-196, 212-215, 227, 239, 251, 260-264, 275-279, 290, 301-306, 321, 341-355, 376-390, 425, 464-467 |
| PyMieSim/polarization.py                       |       45 |       15 |        8 |        1 |     58% |65, 71, 74-84, 87, 90-93, 114, 135 |
| PyMieSim/single/detector/base.py               |       77 |       22 |       10 |        3 |     69% |88, 102, 113, 176, 207-230, 251-288 |
| PyMieSim/single/detector/coherent.py           |       39 |       14 |       10 |        3 |     57% |49, 54, 63-68, 97-115 |
| PyMieSim/single/detector/uncoherent.py         |       26 |        2 |        0 |        0 |     92% |   67, 104 |
| PyMieSim/single/representations.py             |      168 |      112 |       12 |        0 |     31% |58, 61-65, 90-110, 135-155, 244-275, 302, 333-373, 401, 437-461, 485-510, 532-534, 544, 558-581, 615, 628-658, 688-703, 718-734 |
| PyMieSim/single/scatterer/base.py              |      101 |       22 |        8 |        3 |     77% |48, 58, 64-68, 73, 98, 103, 108, 118, 128, 133, 138, 143, 148, 153, 222, 325, 373, 423, 451 |
| PyMieSim/single/scatterer/core\_shell.py       |       25 |        2 |        0 |        0 |     92% |   83, 101 |
| PyMieSim/single/scatterer/cylinder.py          |       38 |       13 |        0 |        0 |     66% |36-40, 46-48, 78, 103, 128, 153, 157, 161, 165, 169 |
| PyMieSim/single/scatterer/sphere.py            |       26 |        4 |        0 |        0 |     85% |80, 107, 135, 163 |
| PyMieSim/single/source/base.py                 |       27 |        5 |       10 |        5 |     73% |23, 33, 43, 53, 58 |
| PyMieSim/single/source/gaussian.py             |       35 |       10 |        2 |        1 |     70% |41, 72-85, 106-118 |
| PyMieSim/single/source/planewave.py            |       26 |       10 |        2 |        0 |     57% |36-41, 65-74, 94-103 |
| PyMieSim/special\_functions.py                 |       20 |       15 |        2 |        0 |     23% |21, 42-45, 66-72, 94-99 |
| PyMieSim/units.py                              |       19 |        0 |        4 |        0 |    100% |           |
|                                      **TOTAL** | **1136** |  **348** |  **158** |   **40** | **66%** |           |


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