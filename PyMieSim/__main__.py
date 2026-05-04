"""Entry point for launching the PyMieSim experiment dashboard."""

import argparse
import logging

from PyMieSim.gui.interface import OpticalSetupGUI


def _build_argument_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for the dashboard launcher."""
    parser = argparse.ArgumentParser(description="Launch the PyMieSim experiment dashboard.")
    parser.add_argument("--host", default="127.0.0.1", help="Host interface for the Dash server.")
    parser.add_argument("--port", default="8050", help="Port for the Dash server.")
    parser.add_argument("--debug", action="store_true", help="Enable Dash debug mode and verbose debug logging.")
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Start the server without opening a browser window.",
    )
    return parser

if __name__ == "__main__":
    args = _build_argument_parser().parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )
    _gui = OpticalSetupGUI()
    _gui.run(
        host=args.host,
        port=args.port,
        open_browser=not args.no_browser,
        debug=args.debug,
    )
