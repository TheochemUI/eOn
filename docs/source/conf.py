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
    with open("../../CHANGELOG.md", "r") as f:
        # The first non-empty line should be like: ## [2.8.0] - 2025-09-04
        for line in f:
            if line.strip():
                # Extract version from "## [2.8.0]" -> "2.8.0"
                release = line.split("[")[1].split("]")[0]
                break
except FileNotFoundError:
    release = "0.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
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

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
}

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_title = "eOn"
html_static_path = ["_static"]
html_baseurl = 'https://eondocs.org/'

html_theme_options = {
    "repository_url": "https://github.com/TheochemUI/eOn",
    "use_repository_button": True,
    "logo": {
        "image_light": "_static/logo/ev2_trans.png",
        "image_dark": "_static/logo/ev2_trans_dark.png",
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
