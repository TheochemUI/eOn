# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "eON"
copyright = "2024, Henkelmann Group, Jonsson Group, Rohit Goswami"
author = "Rohit Goswami"
release = "2.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "autodoc2",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.githubpages",
    "sphinx_contributors",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxcontrib.spelling",
]

autodoc2_render_plugin = "myst"
autodoc2_packages = [
    f"../../src/{project}",
]

myst_enable_extensions = [
    "deflist",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
}

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_title = "eON"
html_static_path = ["_static"]

html_theme_options = {
    "repository_url": "https://github.com/TheochemUI/EONgit",
    "use_repository_button": True,
   "logo": {
      "image_light": "_static/ev2_trans.png",
      "image_dark": "_static/ev2_trans_dark.png",
   }
}
