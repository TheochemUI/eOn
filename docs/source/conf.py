# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from datetime import datetime

project = "eOn"
author = "the eOn developers"
copyright = f"{datetime.now().date().year}, {author}"
try:
    import tomllib

    with open("../../pyproject.toml", "rb") as f:
        release = tomllib.load(f)["project"]["version"]
except Exception:
    release = "0.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.githubpages",
    "sphinx.ext.viewcode",
    "sphinx_togglebutton",
    "sphinx_favicon",
    "sphinx_contributors",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxcontrib.spelling",
    "sphinxcontrib.autodoc_pydantic",
    "sphinxcontrib.bibtex",
    "sphinx_sitemap",
    "autodoc2",
]

# pyeonclient pure-Python modules (steps, ase_bridge, bridge) for automodule.
# Mock the nanobind extension when docs build without a compiled wheel.
autodoc_mock_imports = ["pyeonclient._core"]
import sys
from pathlib import Path as _Path

_pyeon_src = _Path(__file__).resolve().parents[2] / "client" / "python"
if _pyeon_src.is_dir():
    sys.path.insert(0, str(_pyeon_src))

sitemap_show_lastmod = True

bibtex_bibfiles = ["bibtex/eonDocs.bib"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
}

myst_enable_extensions = [
    "deflist",
    "fieldlist",
]

# CI sets NB_EXECUTION_MODE=cache (still executes on miss; faster reruns).
# Local/default remains "force" so authors always see live notebook output.
import os as _os
from pathlib import Path as _Path

nb_execution_mode = _os.environ.get("NB_EXECUTION_MODE", "force")
nb_execution_timeout = 600
nb_execution_raise_on_error = False
nb_execution_show_tb = True
# Persist myst-nb execution results under docs/source for actions/cache restore.
nb_execution_cache_path = str(_Path(__file__).resolve().parent / ".jupyter_cache")

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
}

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "shibuya"
html_title = "eOn"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_baseurl = "https://eondocs.org/"

html_context = {
    "source_type": "github",
    "source_user": "TheochemUI",
    "source_repo": "eOn",
    "source_version": "main",
    "source_docs_path": "/docs/source/",
}

html_theme_options = {
    "github_url": "https://github.com/TheochemUI/eOn",
    "accent_color": "teal",
    "dark_code": True,
    "globaltoc_expand_depth": 2,
    "nav_links": [
        {
            "title": "C++ API",
            "url": "/api-cpp/index.html",
            "resource": True,
        },
        {
            "title": "Ecosystem",
            "children": [
                {
                    "title": "rgpycrumbs",
                    "url": "https://rgpycrumbs.rgoswami.me",
                    "summary": "Python helpers for eOn workflows",
                    "external": True,
                },
                {
                    "title": "chemparseplot",
                    "url": "https://chemparseplot.rgoswami.me",
                    "summary": "Parsing and plotting for computational chemistry",
                    "external": True,
                },
            ],
        },
    ],
    "logo": {
        "light": "_static/logo/eon_v3_light.svg",
        "dark": "_static/logo/eon_v3_dark.svg",
    },
}

# Each configuration model is a field list, in class order, with its default.
autodoc_pydantic_model_show_json = False
autodoc_pydantic_model_show_config_summary = False
autodoc_pydantic_model_show_validator_summary = False
autodoc_pydantic_model_show_validator_members = False
autodoc_pydantic_model_hide_paramlist = True
autodoc_pydantic_model_show_field_summary = True
autodoc_pydantic_model_summary_list_order = "bysource"
autodoc_pydantic_model_member_order = "bysource"
autodoc_pydantic_field_list_validators = False

# The PDF is a manual: logo and title centered, no blank verso, no wrap arrows.
latex_engine = "xelatex"
latex_logo = "_static/favicons/android-chrome-512x512.png"
latex_elements = {
    "extraclassoptions": "openany",
    "fncychap": "",
    "sphinxsetup": (
        r"verbatimvisiblespace=\mbox{}, "
        r"verbatimcontinued=\mbox{}, "
        r"verbatimhintsturnover=false"
    ),
    "preamble": r"""
\makeatletter
\renewcommand{\sphinxmaketitle}{%
  \let\sphinxrestorepageanchorsetting\relax
  \ifHy@pageanchor\def\sphinxrestorepageanchorsetting{\Hy@pageanchortrue}\fi
  \hypersetup{pageanchor=false}%
  \begin{titlepage}%
    \let\footnotesize\small
    \let\footnoterule\relax
    \begingroup
      \def\endgraf{ }\def\and{\& }%
      \pdfstringdefDisableCommands{\def\\{, }}%
      \hypersetup{pdfauthor={\@author}, pdftitle={\@title}}%
    \endgroup
    \vspace*{1.6cm}%
    \begin{center}%
      \sphinxincludegraphics[height=3.2cm]{android-chrome-512x512.png}\par
      \vspace{1.2cm}%
      {\py@HeaderFamily\Huge \@title \par}
      \vspace{0.45cm}%
      {\large long timescale dynamics software \par}
      \vspace{0.7cm}%
      {\itshape\large \py@release \releaseinfo \par}
      \vspace{1.4cm}%
      {\Large
        \begin{tabular}[t]{c}
          \@author
        \end{tabular}\kern-\tabcolsep \par}
      \vfill
      {\large \@date \par}
      \py@authoraddress \par
    \end{center}%
    \@thanks
  \end{titlepage}%
  \setcounter{footnote}{0}%
  \let\thanks\relax\let\maketitle\relax
  \clearpage
  \ifdefined\sphinxbackoftitlepage\sphinxbackoftitlepage\fi
  \if@openright\cleardoublepage\else\clearpage\fi
  \sphinxrestorepageanchorsetting
}
\makeatother
""",
}

# --- Plugin options -------------------------------------------------

# --- autodoc2 options
autodoc2_render_plugin = "myst"
autodoc2_packages = [
    {
        "path": f"../../{project.lower()}",
        "exclude_dirs": [
            "__pycache__",
        ],
        "exclude_files": [
            "*schema*",
        ],
    }
]
autodoc2_hidden_objects = ["dunder", "private", "inherited"]

# ----- sphinx_favicon options

favicons = [
    "favicons/favicon-16x16.png",
    "favicons/favicon-32x32.png",
    "favicons/favicon.ico",
    "favicons/android-chrome-192x192.png",
    "favicons/android-chrome-512x512.png",
    "favicons/apple-touch-icon.png",
    "favicons/site.webmanifest",
]


# --- Generated figures (JSON SSoT → SVG at build time) ---------------------
def _generate_docs_figures(app):
    """Build metatomic backend benchmark figure from committed JSON."""
    import sys

    from sphinx.util import logging as _sphinx_logging

    log = _sphinx_logging.getLogger(__name__)
    repo = _Path(__file__).resolve().parents[2]
    scripts = repo / "scripts"
    if str(scripts) not in sys.path:
        sys.path.insert(0, str(scripts))
    try:
        from plot_metatomic_backend_bench import generate_docs_assets
    except Exception as exc:  # pragma: no cover
        log.warning("metatomic backend figure not generated: %s", exc)
        return
    json_path = repo / "docs/source/fig/data/metatomic_backend_bench.json"
    if not json_path.is_file():
        log.warning("missing benchmark JSON: %s", json_path)
        return
    try:
        svg, table = generate_docs_assets(json_path)
        log.info("generated %s", svg)
        log.info("generated %s", table)
    except Exception as exc:
        log.warning("metatomic backend figure failed: %s", exc)


def setup(app):
    app.connect("builder-inited", _generate_docs_figures)
