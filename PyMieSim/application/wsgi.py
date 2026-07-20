"""WSGI entry point for production servers such as Gunicorn.

Render can load the dashboard with::

    gunicorn PyMieSim.application.wsgi:server --bind 0.0.0.0:$PORT
"""

from PyMieSim.gui.interface import app


server = app.server

