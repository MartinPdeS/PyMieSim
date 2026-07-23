# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                          |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|---------------------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyMieSim/application/wsgi.py                  |        2 |        2 |        0 |        0 |      0.00% |      8-11 |
| PyMieSim/experiment/dataframe\_subclass.py    |      375 |       86 |      180 |       46 |     72.25% |112, 148, 162, 165, 168, 171, 174-175, 177-\>180, 260-261, 312, 338-343, 391, 403, 405-\>408, 409, 412, 436-\>441, 441-\>exit, 465, 484-485, 613-\>616, 640, 643, 646, 677-685, 689-690, 707, 729-732, 749-751, 757, 784-785, 800, 814-815, 822, 825, 835, 924-925, 951, 954, 957, 987-994, 998-999, 1016, 1038-1041, 1058-1060, 1066, 1079, 1122-1123, 1152-1153, 1160, 1163, 1173, 1224-1230, 1246-1248 |
| PyMieSim/experiment/setup.py                  |      102 |        3 |       40 |        5 |     94.37% |103-\>106, 219, 245, 257, 322-\>316 |
| PyMieSim/experiment/utils.py                  |       32 |       32 |       14 |        0 |      0.00% |     1-127 |
| PyMieSim/single/representations/farfields.py  |      207 |       95 |       60 |       18 |     50.19% |193-204, 268-300, 324-336, 370-382, 453-\>exit, 481-500, 525-542, 566-585, 614-653, 670-696, 776, 796, 799, 802-813, 822, 827, 833, 840, 847-857, 870-\>873, 874, 895-897, 911-\>exit, 929-935, 952-956, 987, 1010 |
| PyMieSim/single/representations/nearfields.py |      186 |      164 |       62 |        0 |      8.87% |33-36, 53-66, 107-129, 163-182, 201-210, 240-283, 301-329, 343-355, 372-381, 390-428, 498-556 |
| PyMieSim/single/representations/s1s2.py       |       39 |        4 |       16 |        6 |     78.18% |108-109, 111-\>114, 115, 117-\>120, 159, 161-\>165 |
| PyMieSim/single/representations/spf.py        |       99 |       26 |       32 |       12 |     64.89% |176, 180, 234, 239, 244, 250, 258-275, 289, 308-310, 327-\>exit, 348, 364-374 |
| PyMieSim/single/representations/stokes.py     |       93 |       10 |       24 |       10 |     82.91% |247-\>exit, 279, 300, 305, 311, 326, 347-349, 366-\>exit, 398, 416 |
| **TOTAL**                                     | **1200** |  **422** |  **430** |   **97** | **60.31%** |           |

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