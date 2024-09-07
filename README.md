# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                               |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/\_\_main\_\_.py           |        9 |        9 |        2 |        0 |      0% |      2-15 |
| PyMieSim/experiment/detector.py    |       59 |        1 |       10 |        1 |     97% |       174 |
| PyMieSim/experiment/measure.py     |       37 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/parameters.py  |       23 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer.py   |      118 |        2 |       16 |        2 |     97% |    86, 99 |
| PyMieSim/experiment/setup.py       |       59 |        1 |       26 |        1 |     98% |        93 |
| PyMieSim/experiment/source.py      |       50 |        2 |        8 |        0 |     97% |   137-143 |
| PyMieSim/mesh.py                   |      111 |       15 |       28 |       11 |     81% |66->65, 67, 70->69, 71, 74->73, 75, 78->77, 79, 82->81, 86->85, 93->exit, 100->exit, 104->103, 108->107, 112->111, 179-183, 203-214 |
| PyMieSim/polarization.py           |       50 |        3 |       12 |        1 |     94% |28, 39, 47->46, 67 |
| PyMieSim/single/detector.py        |      120 |        4 |       18 |        3 |     93% |34-35, 358, 363, 375->379 |
| PyMieSim/single/representations.py |      170 |        4 |       24 |        2 |     97% |8, 52, 56-57 |
| PyMieSim/single/scatterer.py       |      172 |        7 |       54 |       22 |     87% |37->36, 42->41, 47->46, 52->51, 57->56, 62->61, 67->66, 72->71, 77->76, 82->81, 87->86, 92->91, 97->96, 102->101, 107->106, 112->111, 117->116, 259, 381, 405, 620->619, 621, 624->623, 625, 628->627, 629, 632->631, 633 |
| PyMieSim/single/source.py          |       50 |        0 |        8 |        0 |    100% |           |
| PyMieSim/special\_functions.py     |       49 |       24 |        4 |        0 |     51% |18-24, 53-56, 70-75, 130-135, 151-156 |
|                          **TOTAL** | **1077** |   **72** |  **210** |   **43** | **91%** |           |


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