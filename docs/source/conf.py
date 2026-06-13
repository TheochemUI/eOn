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

sitemap_show_lastmod = True

bibtex_bibfiles = ["bibtex/eonDocs.bib"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
}

myst_enable_extensions = [
    "deflist",
    "fieldlist",
]

nb_execution_mode = "force"
nb_execution_timeout = 600
nb_execution_raise_on_error = False
nb_execution_show_tb = True

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
html_js_files = [
    ("https://antics-api.turtletech.us/antics.js", {"defer": "defer"}),
]
html_baseurl = 'https://eondocs.org/'

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
