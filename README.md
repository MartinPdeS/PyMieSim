# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyMieSim/experiment/dataframe\_subclass.py     |       65 |       12 |       10 |        1 |     80.00% |150-156, 195, 221-225, 244-248 |
| PyMieSim/experiment/setup.py                   |       99 |        4 |       26 |        4 |     93.60% |24, 296-299, 360->363, 366->370, 375 |
| PyMieSim/experiment/source/base.py             |       21 |        9 |       10 |        0 |     51.61% |     20-33 |
| PyMieSim/experiment/source/gaussian.py         |       19 |        0 |        2 |        1 |     95.24% |    44->48 |
| PyMieSim/experiment/source/planewave.py        |       18 |        0 |        2 |        1 |     95.00% |    40->44 |
| PyMieSim/experiment/utils.py                   |       32 |        1 |       14 |        2 |     93.48% |44->50, 58 |
| PyMieSim/single/mesh.py                        |       38 |        7 |        0 |        0 |     81.58% |94, 106, 118, 130, 142-148 |
| PyMieSim/single/plottings.py                   |      101 |       22 |       18 |        8 |     74.79% |126, 132->135, 135->139, 146->149, 155->153, 157, 175-176, 194-201, 225, 269-276, 300-333, 410-420, 446 |
| PyMieSim/single/representations/base.py        |       32 |       22 |        2 |        0 |     29.41% |53-57, 83-105, 131-153 |
| PyMieSim/single/representations/far\_fields.py |       60 |        2 |        8 |        1 |     95.59% |     78-79 |
| PyMieSim/single/representations/near\_field.py |       71 |       53 |       14 |        0 |     21.18% |67-72, 78-87, 92-94, 115-124, 137-169, 192-240 |
| PyMieSim/single/representations/stokes.py      |       41 |        2 |        4 |        1 |     93.33% |     77-78 |
| **TOTAL**                                      |  **854** |  **134** |  **132** |   **19** | **82.05%** |           |

13 files skipped due to complete coverage.


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