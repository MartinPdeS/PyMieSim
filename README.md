# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/\_\_main\_\_.py                       |        9 |        9 |        2 |        0 |      0% |      2-15 |
| PyMieSim/experiment/detector/base.py           |       57 |        9 |       16 |        5 |     75% |81-85, 104, 107, 127, 130, 149->152 |
| PyMieSim/experiment/detector/coherent\_mode.py |       11 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/detector/photodiode.py     |       11 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       47 |        9 |       14 |        5 |     74% |47-49, 59, 87-92, 108, 113 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       21 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/cylinder.py      |       18 |        6 |        0 |        0 |     67% |     44-55 |
| PyMieSim/experiment/scatterer/sphere.py        |       18 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/setup.py                   |       95 |       10 |       26 |        6 |     85% |106, 133, 210, 236-237, 262-267, 302->306, 311 |
| PyMieSim/experiment/source/base.py             |       38 |        6 |       12 |        2 |     76% |33->36, 48-54, 69 |
| PyMieSim/experiment/source/gaussian.py         |       23 |        2 |        6 |        3 |     83% |37, 40, 47->50 |
| PyMieSim/experiment/source/planewave.py        |       15 |        4 |        2 |        0 |     65% | 32-35, 45 |
| PyMieSim/mesh.py                               |       94 |        8 |        8 |        2 |     90% |117, 129, 141, 153, 195->exit, 214->exit, 275-279 |
| PyMieSim/polarization.py                       |       45 |        5 |        8 |        1 |     85% |71, 87, 90-93 |
| PyMieSim/single/detector/base.py               |       75 |        4 |       10 |        4 |     91% |91, 102, 105, 222 |
| PyMieSim/single/detector/coherent.py           |       41 |        2 |       10 |        3 |     90% |48, 53, 65->69 |
| PyMieSim/single/detector/uncoherent.py         |       28 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/representations.py             |      168 |        3 |       12 |        1 |     98% | 58, 62-63 |
| PyMieSim/single/scatterer/base.py              |      103 |        4 |       10 |        4 |     93% |48, 51, 61, 454 |
| PyMieSim/single/scatterer/core\_shell.py       |       25 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/scatterer/cylinder.py          |       38 |        4 |        0 |        0 |     89% |157, 161, 165, 169 |
| PyMieSim/single/scatterer/sphere.py            |       26 |        2 |        0 |        0 |     92% |  135, 163 |
| PyMieSim/single/source/base.py                 |       33 |        7 |       16 |        7 |     71% |22, 25, 34, 37, 46, 49, 63 |
| PyMieSim/single/source/gaussian.py             |       35 |        1 |        2 |        1 |     95% |        41 |
| PyMieSim/single/source/planewave.py            |       26 |        1 |        2 |        1 |     93% |        37 |
| PyMieSim/special\_functions.py                 |       20 |        0 |        2 |        0 |    100% |           |
| PyMieSim/units.py                              |       19 |        0 |        4 |        0 |    100% |           |
|                                      **TOTAL** | **1139** |   **96** |  **162** |   **45** | **88%** |           |


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