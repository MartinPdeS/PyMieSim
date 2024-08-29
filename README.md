# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                               |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/\_\_main\_\_.py           |        9 |        9 |        2 |        0 |      0% |      2-15 |
| PyMieSim/experiment/detector.py    |       59 |        6 |       10 |        0 |     86% |   169-176 |
| PyMieSim/experiment/measure.py     |       37 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/parameters.py  |       23 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer.py   |      118 |        1 |       16 |        1 |     99% |        86 |
| PyMieSim/experiment/setup.py       |       59 |        1 |       26 |        1 |     98% |        93 |
| PyMieSim/experiment/source.py      |       50 |        2 |        8 |        0 |     97% |   137-143 |
| PyMieSim/mesh.py                   |      111 |       15 |       28 |       11 |     81% |66->65, 67, 70->69, 71, 74->73, 75, 78->77, 79, 82->81, 86->85, 93->exit, 100->exit, 104->103, 108->107, 112->111, 179-183, 203-214 |
| PyMieSim/physics.py                |       18 |       18 |        0 |        0 |      0% |      4-56 |
| PyMieSim/polarization.py           |       50 |        3 |       12 |        1 |     94% |28, 39, 47->46, 67 |
| PyMieSim/single/detector.py        |      112 |       14 |       20 |        3 |     84% |33-34, 83-106, 288, 293, 305->309 |
| PyMieSim/single/representations.py |      165 |       81 |       22 |        1 |     49% |8, 52, 55-59, 81-101, 123-143, 224-255, 305-345, 394-419, 464-487, 617-633 |
| PyMieSim/single/scatterer.py       |      170 |       14 |       54 |       23 |     83% |28-39, 42->41, 46->45, 50->49, 55->54, 60->59, 65->64, 70->69, 75->74, 80->79, 85->84, 90->89, 95->94, 100->99, 105->104, 110->109, 115->114, 120->119, 248, 259, 265, 375, 398, 601->600, 602, 605->604, 606, 609->608, 610, 613->612, 614 |
| PyMieSim/single/source.py          |       46 |       17 |        8 |        1 |     63% |39-44, 57-61, 83->86, 102-121 |
| PyMieSim/special\_functions.py     |       49 |       24 |        4 |        1 |     49% |18-24, 53-56, 70-75, 108->111, 130-135, 151-156 |
|                          **TOTAL** | **1076** |  **205** |  **210** |   **43** | **79%** |           |


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