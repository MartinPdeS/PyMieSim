"""Reusable card primitives shared by all dashboard pages."""

from __future__ import annotations

from collections.abc import Sequence

from dash import html


class Card:
    """Build a RosettaX-style card with an optional colored left accent."""

    @classmethod
    def classes(cls, *, color: str = "blue", extra: str = "") -> str:
        """Return the shared card class names for any page component."""
        return " ".join(filter(None, ["ui-card", f"ui-card--{color}", extra]))

    def __init__(self, children, *, class_name: str = "", color: str = "blue"):
        self.children = children
        self.class_name = class_name
        self.color = color

    def render(self):
        return html.Section(className=self.classes(color=self.color, extra=self.class_name), children=self.children)


class HeaderCard(Card):
    """Build a workflow header with an intro and embedded step cards."""

    def __init__(self, title: str, description: str, steps: Sequence[tuple[str, str, str, str]], *, color: str = "green"):
        super().__init__(children=[], class_name="workflow-header-card", color=color)
        self.title = title
        self.description = description
        self.steps = steps

    def render(self):
        return html.Section(
            className=self.classes(color=self.color, extra=self.class_name),
            children=[
                html.Div(
                    className="header-card-intro",
                    children=[html.H1(self.title), html.P(self.description)],
                ),
                html.Div(
                    className="header-card-steps",
                    children=[
                        html.Div(
                            className=f"header-step-card header-step-card--{color}",
                            children=[html.Span(number, className="header-step-number"), html.H3(step_title), html.P(step_description)],
                        )
                        for number, step_title, step_description, color in self.steps
                    ],
                ),
            ]
        )
