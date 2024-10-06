# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/\_\_main\_\_.py                       |        9 |        9 |        2 |        0 |      0% |      2-15 |
| PyMieSim/experiment/detector/base.py           |       56 |        4 |       18 |        8 |     84% |64->63, 70, 73, 78->77, 84, 87, 92->91, 94->97 |
| PyMieSim/experiment/detector/coherent\_mode.py |       20 |        5 |        8 |        1 |     64% |47->46, 49-53 |
| PyMieSim/experiment/detector/photodiode.py     |       11 |        0 |        4 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       37 |        5 |       16 |        4 |     77% |40->39, 69-74, 88->85, 97 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       38 |        2 |       14 |        6 |     85% |53->52, 59->58, 61->64, 67->66, 71, 74 |
| PyMieSim/experiment/scatterer/cylinder.py      |       32 |       14 |       12 |        2 |     50% |41->40, 43-46, 49->48, 52-58, 67-78 |
| PyMieSim/experiment/scatterer/sphere.py        |       32 |        2 |       12 |        5 |     84% |44->43, 46->49, 52->51, 56, 59 |
| PyMieSim/experiment/setup.py                   |       98 |       12 |       30 |        6 |     84% |108, 135, 212, 245-246, 272-277, 297-301, 329->334, 339 |
| PyMieSim/experiment/source/base.py             |       39 |        2 |       16 |        5 |     87% |33->36, 45->44, 49, 52, 68->65 |
| PyMieSim/experiment/source/gaussian.py         |       23 |        2 |       12 |        5 |     80% |32->31, 36, 39, 44->43, 46->49 |
| PyMieSim/experiment/source/planewave.py        |       21 |        0 |        6 |        2 |     93% |29->28, 31->34 |
| PyMieSim/mesh.py                               |      111 |       15 |       28 |       11 |     81% |66->65, 67, 70->69, 71, 74->73, 75, 78->77, 79, 82->81, 86->85, 93->exit, 100->exit, 104->103, 108->107, 112->111, 179-183, 203-214 |
| PyMieSim/polarization.py                       |       45 |        3 |       12 |        2 |     91% |30, 46, 52 |
| PyMieSim/single/detector/base.py               |       72 |        3 |       14 |        5 |     91% |45->44, 54, 59->58, 65, 68 |
| PyMieSim/single/detector/coherent.py           |       42 |        2 |       12 |        3 |     91% |49, 54, 66->70 |
| PyMieSim/single/detector/uncoherent.py         |       28 |        0 |        4 |        0 |    100% |           |
| PyMieSim/single/representations.py             |      170 |        4 |       24 |        2 |     97% |8, 52, 56-57 |
| PyMieSim/single/scatterer/base.py              |       89 |        1 |       40 |       18 |     85% |50->49, 55->54, 60->59, 65->64, 70->69, 75->74, 80->79, 85->84, 90->89, 95->94, 100->99, 105->104, 110->109, 115->114, 120->119, 125->124, 130->129, 275 |
| PyMieSim/single/scatterer/core\_shell.py       |       41 |        9 |       16 |        4 |     67% |53->52, 58, 61, 66->65, 70-79 |
| PyMieSim/single/scatterer/cylinder.py          |       54 |       13 |       24 |        8 |     65% |40->39, 45, 48, 53->52, 57-66, 171->170, 172, 175->174, 176, 179->178, 180, 183->182, 184 |
| PyMieSim/single/scatterer/sphere.py            |       42 |       11 |       16 |        4 |     64% |40->39, 45, 48, 53->52, 57-66, 154, 178 |
| PyMieSim/single/source/base.py                 |       33 |        7 |       24 |       11 |     68% |18->17, 22, 25, 30->29, 34, 37, 42->41, 46, 49, 54->53, 63 |
| PyMieSim/single/source/gaussian.py             |       33 |        1 |        4 |        1 |     95% |        42 |
| PyMieSim/single/source/planewave.py            |       26 |        1 |        4 |        1 |     93% |        37 |
| PyMieSim/special\_functions.py                 |       49 |       24 |        4 |        0 |     51% |18-24, 53-56, 70-75, 130-135, 151-156 |
| PyMieSim/units.py                              |       24 |        0 |        4 |        0 |    100% |           |
| PyMieSim/utils.py                              |       39 |        4 |        8 |        2 |     87% |79-80, 162-163 |
|                                      **TOTAL** | **1314** |  **155** |  **388** |  **116** | **82%** |           |


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