"""Test-only: a copy of ``rigel.config.CONSTANTS`` with one value changed, the way an experiment varies a constant
without editing a module."""

from dataclasses import replace

from rigel.config import Constants


def replaced(constants: Constants, dotted: str, value) -> Constants:
    """A copy of *constants* with one value changed, named ``"section.field"``; the section re-validates."""
    section_name, _, field_name = dotted.partition(".")
    section = getattr(constants, section_name)
    return replace(constants, **{section_name: replace(section, **{field_name: value})})
