# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/\_\_main\_\_.py                       |        9 |        9 |        2 |        0 |      0% |      2-15 |
| PyMieSim/experiment/detector/base.py           |       56 |        3 |       18 |        6 |     88% |54->53, 63, 68->67, 74, 77, 82->81 |
| PyMieSim/experiment/detector/coherent\_mode.py |       20 |        1 |        8 |        2 |     89% |37->36, 42 |
| PyMieSim/experiment/detector/photodiode.py     |       13 |        0 |        4 |        0 |    100% |           |
| PyMieSim/experiment/measure.py                 |       37 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/parameters.py              |       23 |        0 |        0 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base\_class.py   |       43 |       10 |       16 |        2 |     76% |45-54, 76, 91 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       42 |        2 |       12 |        5 |     87% |51->50, 53->56, 59->58, 63, 66 |
| PyMieSim/experiment/scatterer/cylinder.py      |       38 |        2 |       12 |        5 |     86% |39->38, 41->44, 47->46, 51, 54 |
| PyMieSim/experiment/scatterer/sphere.py        |       38 |        2 |       12 |        5 |     86% |44->43, 46->49, 52->51, 56, 59 |
| PyMieSim/experiment/setup.py                   |       97 |        8 |       30 |        3 |     90% |108, 135, 272-277, 297-301, 329->334 |
| PyMieSim/experiment/source/base.py             |       35 |        6 |       12 |        3 |     72% |32->35, 44->43, 47-53, 73 |
| PyMieSim/experiment/source/gaussian.py         |       25 |        2 |       12 |        5 |     81% |30->29, 34, 37, 42->41, 44->47 |
| PyMieSim/experiment/source/planewave.py        |       10 |        1 |        2 |        0 |     92% |        33 |
| PyMieSim/mesh.py                               |      111 |       15 |       28 |       11 |     81% |66->65, 67, 70->69, 71, 74->73, 75, 78->77, 79, 82->81, 86->85, 93->exit, 100->exit, 104->103, 108->107, 112->111, 179-183, 203-214 |
| PyMieSim/polarization.py                       |       41 |        4 |       10 |        2 |     84% |29->32, 34->37, 40, 43-46 |
| PyMieSim/single/detector.py                    |      129 |      129 |       30 |        0 |      0% |     4-442 |
| PyMieSim/single/detector/base.py               |       73 |        3 |       14 |        5 |     91% |46->45, 55, 60->59, 66, 69 |
| PyMieSim/single/detector/coherent.py           |       42 |        2 |       12 |        3 |     91% |39, 44, 56->60 |
| PyMieSim/single/detector/uncoherent.py         |       29 |        0 |        4 |        0 |    100% |           |
| PyMieSim/single/representations.py             |      170 |        4 |       24 |        2 |     97% |8, 52, 56-57 |
| PyMieSim/single/scatterer/base.py              |       88 |        1 |       40 |       18 |     85% |35->34, 40->39, 45->44, 50->49, 55->54, 60->59, 65->64, 70->69, 75->74, 80->79, 85->84, 90->89, 95->94, 100->99, 105->104, 110->109, 115->114, 257 |
| PyMieSim/single/scatterer/core\_shell.py       |       49 |        5 |       16 |        7 |     82% |49->48, 54, 57, 62->61, 67, 70, 73 |
| PyMieSim/single/scatterer/cylinder.py          |       60 |        9 |       24 |       11 |     76% |44->43, 49, 52, 57->56, 62, 65, 68, 175->174, 176, 179->178, 180, 183->182, 184, 187->186, 188 |
| PyMieSim/single/scatterer/sphere.py            |       48 |        7 |       16 |        7 |     78% |42->41, 47, 50, 55->54, 60, 63, 66, 159, 183 |
| PyMieSim/single/source/base.py                 |       33 |        7 |       24 |       11 |     68% |18->17, 22, 25, 30->29, 34, 37, 42->41, 46, 49, 54->53, 63 |
| PyMieSim/single/source/gaussian.py             |       33 |        1 |        4 |        1 |     95% |        38 |
| PyMieSim/single/source/planewave.py            |       26 |        1 |        4 |        1 |     93% |        33 |
| PyMieSim/special\_functions.py                 |       49 |       24 |        4 |        0 |     51% |18-24, 53-56, 70-75, 130-135, 151-156 |
| PyMieSim/units.py                              |       24 |        0 |        4 |        0 |    100% |           |
| PyMieSim/utils.py                              |       77 |       36 |       20 |        2 |     48% |29-54, 74-120, 132-155, 226-227, 309-310 |
|                                      **TOTAL** | **1568** |  **294** |  **418** |  **117** | **76%** |           |


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