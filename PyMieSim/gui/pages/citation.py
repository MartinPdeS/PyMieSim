"""Citation page composition."""

from dash import html

from PyMieSim.gui.components import Card


PYMIESIM_DOI = "10.1364/OPTCON.473102"
PYMIESIM_ARTICLE_URL = "https://opg.optica.org/optcon/fulltext.cfm?uri=optcon-2-3-520"
PYMIESIM_BIBTEX = """@article{PoinsinetdeSivry-Houle:23,
  author = {Martin Poinsinet de Sivry-Houle and Nicolas Godbout and Caroline Boudoux},
  journal = {Opt. Continuum},
  title = {PyMieSim: an open-source library for fast and flexible far-field Mie scattering simulations},
  volume = {2},
  number = {3},
  pages = {520--534},
  year = {2023},
  doi = {10.1364/OPTCON.473102},
}"""


def build_citation_page():
    """Build the citation page shown by the home-page citation button."""
    return html.Div(
        className="page-content-stack",
        children=[
            html.Section(
                className="page-hero documentation-page-hero",
                children=[
                    html.P("PyMieSim reference", className="eyebrow"),
                    html.H1("Cite PyMieSim"),
                    html.P(
                        "Use this reference when PyMieSim contributes to a scientific workflow, analysis, or publication.",
                        className="hero-text",
                    ),
                ],
            ),
            html.Section(
                className=Card.classes(color="blue", extra="panel documentation-card"),
                children=[
                    html.Div(className="card-header panel-header", children=[html.H2("BibTeX")]),
                    html.Div(className="card-body", children=[
                        html.Div("Copy this reference into your bibliography file.", className="home-info-copy"),
                        html.Pre(PYMIESIM_BIBTEX, className="citation-bibtex"),
                    ]),
                ],
            ),
            html.Section(
                className=Card.classes(color="blue", extra="panel documentation-card"),
                children=[
                    html.Div(className="card-header panel-header", children=[html.H2("Publication")]),
                    html.Div(className="card-body", children=[
                        html.Div(f"DOI: {PYMIESIM_DOI}", className="home-info-copy"),
                        html.A("View publication", href=PYMIESIM_ARTICLE_URL, target="_blank", rel="noopener noreferrer", className="inline-action"),
                    ]),
                ],
            ),
            html.Section(
                className=Card.classes(color="blue", extra="panel documentation-card"),
                children=[
                    html.Div(className="card-header panel-header", children=[html.H2("Citation guidance")]),
                    html.Div(className="card-body", children=[html.P("Please cite the publication above when using PyMieSim in academic work.")]),
                ],
            ),
        ],
    )
