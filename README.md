# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyMieSim/experiment/dataframe\_subclass.py     |       65 |       12 |       10 |        1 |     80.00% |150-156, 195, 221-225, 244-248 |
| PyMieSim/experiment/scatterer/base.py          |       38 |        2 |       12 |        2 |     92.00% |    77, 94 |
| PyMieSim/experiment/setup.py                   |      102 |        3 |       26 |        4 |     94.53% |302-305, 366->369, 372->376, 381 |
| PyMieSim/experiment/source/base.py             |       40 |        2 |       10 |        3 |     90.00% |47, 56->59, 73 |
| PyMieSim/experiment/utils.py                   |       32 |        1 |       14 |        2 |     93.48% |44->50, 58 |
| PyMieSim/single/detector/base.py               |       30 |        2 |        0 |        0 |     93.33% |   28, 197 |
| PyMieSim/single/detector/coherent.py           |       23 |        2 |        4 |        2 |     85.19% |    72, 79 |
| PyMieSim/single/mesh.py                        |       46 |        7 |        0 |        0 |     84.78% |100, 112, 124, 136, 160-166 |
| PyMieSim/single/polarization.py                |       43 |        5 |        8 |        1 |     84.31% |71, 87, 90-93 |
| PyMieSim/single/representations/base.py        |       38 |        3 |        2 |        1 |     90.00% | 50, 54-55 |
| PyMieSim/single/representations/near\_field.py |       68 |       53 |       14 |        0 |     18.29% |66-70, 75-84, 89-91, 112-121, 134-166, 189-237 |
| PyMieSim/single/scatterer/base.py              |      127 |       27 |       12 |        1 |     74.10% |21-25, 502-539, 575, 673 |
| PyMieSim/single/scatterer/core\_shell.py       |       27 |        2 |        0 |        0 |     92.59% |   119-124 |
| PyMieSim/single/scatterer/cylinder.py          |       33 |        7 |        0 |        0 |     78.79% |69, 73, 77, 81, 85, 107-116 |
| PyMieSim/single/scatterer/sphere.py            |       25 |        2 |        0 |        0 |     92.00% |    77, 82 |
| **TOTAL**                                      | **1147** |  **130** |  **134** |   **17** | **86.49%** |           |

18 files skipped due to complete coverage.


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