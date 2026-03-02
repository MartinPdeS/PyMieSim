# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyMieSim/experiment/dataframe\_subclass.py     |       65 |       12 |       10 |        1 |     80.00% |150-156, 195, 221-225, 244-248 |
| PyMieSim/experiment/setup.py                   |       88 |        3 |       26 |        3 |     94.74% |281-284, 351->355, 360 |
| PyMieSim/experiment/utils.py                   |       32 |        2 |       14 |        4 |     86.96% |44->50, 58, 66->63, 123 |
| PyMieSim/single/mesh.py                        |       38 |        7 |        0 |        0 |     81.58% |94, 106, 118, 130, 142-148 |
| PyMieSim/single/plottings.py                   |      101 |       22 |       18 |        8 |     74.79% |126, 132->135, 135->139, 146->149, 155->153, 157, 175-176, 194-201, 225, 269-276, 300-333, 410-420, 446 |
| PyMieSim/single/representations/base.py        |       32 |       22 |        2 |        0 |     29.41% |53-57, 83-105, 131-153 |
| PyMieSim/single/representations/far\_fields.py |       58 |        2 |        8 |        1 |     95.45% |     78-79 |
| PyMieSim/single/representations/near\_field.py |      156 |      137 |       42 |        0 |      9.60% |29-32, 49-62, 96-115, 134-143, 173-213, 231-251, 268-277, 286-324, 368-425 |
| PyMieSim/single/representations/stokes.py      |       41 |        2 |        4 |        1 |     93.33% |     77-78 |
| **TOTAL**                                      |  **809** |  **209** |  **136** |   **18** | **71.11%** |           |

9 files skipped due to complete coverage.


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