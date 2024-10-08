# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| PyMieSim/\_\_main\_\_.py                       |        9 |        9 |        2 |        0 |      0% |      2-15 |
| PyMieSim/experiment/detector/base.py           |       56 |        4 |       18 |        8 |     84% |64->63, 70, 73, 78->77, 84, 87, 92->91, 94->97 |
| PyMieSim/experiment/detector/coherent\_mode.py |       20 |        5 |        8 |        1 |     64% |47->46, 49-53 |
| PyMieSim/experiment/detector/photodiode.py     |       11 |        0 |        4 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/base.py          |       48 |        8 |       22 |        7 |     74% |40->39, 49-51, 54->53, 61, 90-95, 109->106, 118 |
| PyMieSim/experiment/scatterer/core\_shell.py   |       23 |        0 |        2 |        0 |    100% |           |
| PyMieSim/experiment/scatterer/cylinder.py      |       20 |        6 |        2 |        0 |     73% |     47-58 |
| PyMieSim/experiment/scatterer/sphere.py        |       20 |        0 |        2 |        0 |    100% |           |
| PyMieSim/experiment/setup.py                   |       94 |       10 |       30 |        6 |     85% |108, 135, 212, 244-245, 271-276, 304->309, 314 |
| PyMieSim/experiment/source/base.py             |       37 |        5 |       16 |        3 |     77% |32->35, 44->43, 47-53, 67->64 |
| PyMieSim/experiment/source/gaussian.py         |       23 |        2 |       12 |        5 |     80% |32->31, 36, 39, 44->43, 46->49 |
| PyMieSim/experiment/source/planewave.py        |       15 |        4 |        6 |        1 |     67% |29->28, 31-34, 44 |
| PyMieSim/mesh.py                               |       94 |        8 |       28 |       11 |     84% |107->106, 116, 119->118, 128, 131->130, 140, 143->142, 152, 155->154, 167->166, 194->exit, 213->exit, 217->216, 229->228, 241->240, 274-278 |
| PyMieSim/polarization.py                       |       45 |        5 |       12 |        1 |     86% |70, 86, 89-92 |
| PyMieSim/single/detector/base.py               |       77 |        4 |       16 |        6 |     89% |82->81, 91, 96->95, 102, 105, 229 |
| PyMieSim/single/detector/coherent.py           |       42 |        2 |       12 |        3 |     91% |49, 54, 66->70 |
| PyMieSim/single/detector/uncoherent.py         |       28 |        0 |        4 |        0 |    100% |           |
| PyMieSim/single/representations.py             |      168 |        3 |       28 |        1 |     98% | 57, 61-62 |
| PyMieSim/single/scatterer/base.py              |      103 |        4 |       50 |       23 |     82% |42->41, 47, 50, 55->54, 60, 73->72, 78->77, 83->82, 88->87, 93->92, 98->97, 103->102, 108->107, 113->112, 118->117, 123->122, 128->127, 133->132, 138->137, 143->142, 148->147, 153->152, 451 |
| PyMieSim/single/scatterer/core\_shell.py       |       26 |        0 |        2 |        0 |    100% |           |
| PyMieSim/single/scatterer/cylinder.py          |       39 |        4 |       10 |        4 |     84% |142->141, 143, 146->145, 147, 150->149, 151, 154->153, 155 |
| PyMieSim/single/scatterer/sphere.py            |       26 |        2 |        2 |        0 |     93% |  124, 148 |
| PyMieSim/single/source/base.py                 |       33 |        7 |       24 |       11 |     68% |18->17, 22, 25, 30->29, 34, 37, 42->41, 46, 49, 54->53, 63 |
| PyMieSim/single/source/gaussian.py             |       35 |        1 |        4 |        1 |     95% |        42 |
| PyMieSim/single/source/planewave.py            |       26 |        1 |        4 |        1 |     93% |        37 |
| PyMieSim/special\_functions.py                 |       20 |        0 |        2 |        0 |    100% |           |
| PyMieSim/units.py                              |       25 |        0 |        4 |        0 |    100% |           |
| PyMieSim/utils.py                              |       39 |        4 |        8 |        2 |     87% |79-80, 162-163 |
|                                      **TOTAL** | **1202** |   **98** |  **334** |   **95** | **86%** |           |


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