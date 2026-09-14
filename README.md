# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                          |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|---------------------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyMieSim/experiment/distributions.py          |       67 |        1 |       26 |        1 |     97.85% |        36 |
| PyMieSim/experiment/setup.py                  |       93 |        3 |       30 |        3 |     95.12% |174, 181, 276 |
| PyMieSim/experiment/utils.py                  |       32 |       32 |       14 |        0 |      0.00% |     1-134 |
| PyMieSim/materials.py                         |      211 |       44 |      100 |       33 |     72.03% |82, 94-95, 101, 103, 124, 126, 128, 149, 157, 160-163, 173, 175, 192-\>200, 198, 218-\>226, 247, 255-261, 262-\>264, 265-\>267, 271, 274-\>276, 283-287, 294, 299, 302-307, 314, 320, 350-\>358, 356, 370-\>372, 380, 384, 385-\>exit, 386-\>exit |
| PyMieSim/measures.py                          |       70 |        1 |       10 |        0 |     98.75% |        37 |
| PyMieSim/results.py                           |       42 |        3 |        2 |        1 |     90.91% |32, 62, 86 |
| PyMieSim/single/api.py                        |       40 |        3 |        4 |        1 |     90.91% |41, 78, 94 |
| PyMieSim/single/representations/\_plotting.py |       63 |        2 |       26 |        3 |     94.38% |48, 50, 91-\>exit |
| PyMieSim/single/representations/farfields.py  |      152 |       77 |       28 |        6 |     46.11% |213-224, 288-320, 344-356, 390-402, 473-\>exit, 501-520, 545-562, 586-605, 634-673, 690-716, 796, 805, 809, 812-813, 834-840, 857-861 |
| PyMieSim/single/representations/nearfields.py |      251 |       86 |       88 |       14 |     58.70% |36, 108-130, 206, 246-\>253, 248-\>250, 250-\>253, 269, 305-333, 363-375, 392-401, 413, 426, 428, 432, 442, 494-\>464, 498, 541, 622-680 |
| PyMieSim/single/representations/s1s2.py       |       39 |        4 |       16 |        6 |     78.18% |109-110, 112-\>115, 116, 118-\>121, 160, 162-\>166 |
| PyMieSim/single/representations/spf.py        |       47 |        2 |        4 |        2 |     92.16% |  187, 191 |
| PyMieSim/single/representations/stokes.py     |       53 |        1 |        4 |        2 |     94.74% |260-\>exit, 292 |
| **TOTAL**                                     | **1225** |  **259** |  **354** |   **72** | **73.97%** |           |

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