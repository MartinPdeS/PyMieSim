# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/\_\_main\_\_.py                       |        9 |        9 |        2 |        0 |      0% |      2-15 |
| PyMieSim/experiment/detector/base.py           |       57 |        8 |       16 |        4 |     78% |81-85, 107, 127, 130, 149->152 |
| PyMieSim/experiment/detector/coherent\_mode.py |       11 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/detector/photodiode.py     |       11 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       47 |        2 |       14 |        2 |     93% |   92, 108 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       21 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/cylinder.py      |       18 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/sphere.py        |       18 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/setup.py                   |       94 |        6 |       26 |        3 |     91% |106, 133, 262-267, 300->304 |
| PyMieSim/experiment/source/base.py             |       58 |        3 |       24 |        6 |     89% |34, 62->65, 75->78, 89, 91->exit, 107 |
| PyMieSim/experiment/source/gaussian.py         |        9 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/source/planewave.py        |       15 |        0 |        2 |        1 |     94% |    32->35 |
| PyMieSim/mesh.py                               |       94 |        8 |        8 |        2 |     90% |117, 129, 141, 153, 195->exit, 214->exit, 275-279 |
| PyMieSim/polarization.py                       |       45 |        1 |        8 |        0 |     98% |        87 |
| PyMieSim/single/detector/base.py               |       78 |        1 |       10 |        1 |     98% |       229 |
| PyMieSim/single/detector/coherent.py           |       41 |        2 |       10 |        3 |     90% |48, 53, 65->69 |
| PyMieSim/single/detector/uncoherent.py         |       28 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/representations.py             |      168 |        3 |       12 |        1 |     98% | 58, 62-63 |
| PyMieSim/single/scatterer/base.py              |      103 |        4 |       10 |        4 |     93% |48, 51, 61, 454 |
| PyMieSim/single/scatterer/core\_shell.py       |       25 |        0 |        0 |        0 |    100% |           |
| PyMieSim/single/scatterer/cylinder.py          |       38 |        4 |        0 |        0 |     89% |157, 161, 165, 169 |
| PyMieSim/single/scatterer/sphere.py            |       26 |        2 |        0 |        0 |     92% |  135, 163 |
| PyMieSim/single/source/base.py                 |       27 |        0 |       10 |        0 |    100% |           |
| PyMieSim/single/source/gaussian.py             |       35 |        1 |        2 |        1 |     95% |        41 |
| PyMieSim/single/source/planewave.py            |       26 |        1 |        2 |        1 |     93% |        37 |
| PyMieSim/special\_functions.py                 |       20 |        0 |        2 |        0 |    100% |           |
| PyMieSim/units.py                              |       19 |        0 |        4 |        0 |    100% |           |
|                                      **TOTAL** | **1141** |   **55** |  **162** |   **29** | **93%** |           |


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