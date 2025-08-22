from pydantic import ConfigDict


# Configuration dictionary for the Pydantic dataclass
config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)
