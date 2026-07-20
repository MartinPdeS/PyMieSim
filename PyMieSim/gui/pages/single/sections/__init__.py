"""Sections composing the single-particle page."""

from .representations import build_representation_section
from .setup import build_scatterer_section, build_source_section

__all__ = ["build_representation_section", "build_scatterer_section", "build_source_section"]
