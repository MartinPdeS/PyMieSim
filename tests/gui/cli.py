#!/usr/bin/env python
"""Regression tests for the dashboard command-line interface."""

from PyMieSim.__main__ import _build_argument_parser
from PyMieSim.gui.interface import create_dash_app


def test_cli_defaults_are_stable():
    args = _build_argument_parser().parse_args([])

    assert args.host == "0.0.0.0"
    assert args.port == "8050"
    assert args.debug is False
    assert args.no_browser is False


def test_cli_accepts_server_options():
    args = _build_argument_parser().parse_args(
        ["--host", "127.0.0.1", "--port", "9000", "--debug", "--no-browser"]
    )

    assert args.host == "127.0.0.1"
    assert args.port == "9000"
    assert args.debug is True
    assert args.no_browser is True


def test_dashboard_registers_sidebar_routes():
    app = create_dash_app()

    assert any("page-content.children" in callback_id for callback_id in app.callback_map)
