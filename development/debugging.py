from pydantic.dataclasses import dataclass
from typing import List

import numpy


@dataclass
class ComplexNumbersDataClass:
    numbers: List[int]
