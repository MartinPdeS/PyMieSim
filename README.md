# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                          |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|---------------------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyMieSim/experiment/setup.py                  |       80 |        3 |       30 |        3 |     94.55% |96, 103, 220 |
| PyMieSim/experiment/utils.py                  |       32 |       32 |       14 |        0 |      0.00% |     1-134 |
| PyMieSim/materials.py                         |      211 |       44 |      100 |       33 |     72.03% |82, 94-95, 101, 103, 124, 126, 128, 149, 157, 160-163, 173, 175, 192-\>200, 198, 218-\>226, 247, 255-261, 262-\>264, 265-\>267, 271, 274-\>276, 283-287, 294, 299, 302-307, 314, 320, 350-\>358, 356, 370-\>372, 380, 384, 385-\>exit, 386-\>exit |
| PyMieSim/measures.py                          |       59 |        6 |        4 |        0 |     84.13% | 37, 47-51 |
| PyMieSim/results.py                           |       29 |        6 |        0 |        0 |     79.31% |23, 29, 40, 43, 53, 62 |
| PyMieSim/single/api.py                        |       36 |        7 |        6 |        2 |     78.57% |30, 42, 52, 57, 62, 67, 72 |
| PyMieSim/single/representations/farfields.py  |      207 |       95 |       60 |       18 |     50.19% |211-222, 286-318, 342-354, 388-400, 471-\>exit, 499-518, 543-560, 584-603, 632-671, 688-714, 794, 814, 817, 820-831, 840, 845, 851, 858, 865-875, 888-\>891, 892, 913-915, 929-\>exit, 947-953, 970-974, 1005, 1028 |
| PyMieSim/single/representations/nearfields.py |      251 |       86 |       88 |       14 |     58.70% |36, 108-130, 206, 246-\>253, 248-\>250, 250-\>253, 269, 305-333, 363-375, 392-401, 413, 426, 428, 432, 442, 494-\>464, 498, 541, 622-680 |
| PyMieSim/single/representations/s1s2.py       |       39 |        4 |       16 |        6 |     78.18% |109-110, 112-\>115, 116, 118-\>121, 160, 162-\>166 |
| PyMieSim/single/representations/spf.py        |       99 |       26 |       32 |       12 |     64.89% |185, 189, 243, 248, 253, 259, 267-284, 298, 317-319, 336-\>exit, 357, 373-383 |
| PyMieSim/single/representations/stokes.py     |       93 |       10 |       24 |       10 |     82.91% |258-\>exit, 290, 311, 316, 322, 337, 358-360, 377-\>exit, 409, 427 |
| **TOTAL**                                     | **1201** |  **319** |  **376** |   **98** | **67.34%** |           |

3 files skipped due to complete coverage.


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