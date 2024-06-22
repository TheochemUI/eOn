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
    "autodoc2",
]

autodoc2_render_plugin = "myst"
autodoc2_packages = [
    f"../../{project.lower()}",
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
    },
}

# --- Plugin options
autodoc2_render_plugin = "myst"

favicons = [
    "favicons/favicon-16x16.png",
    "favicons/favicon-32x32.png",
    "favicons/favicon.ico",
    "favicons/android-chrome-192x192.png",
    "favicons/android-chrome-512x512.png",
    "favicons/apple-touch-icon.png",
    "favicons/site.webmanifest",
]

GH_ORGANIZATION = "TheochemUI"
GH_PROJECT = "EONgit"
MODULE = "eon"


def linkcode_resolve(domain, info):
    import subprocess
    import sys, os

    """Generate link to GitHub.
    References:
    - https://github.com/scikit-learn/scikit-learn/blob/f0faaee45762d0a5c75dcf3d487c118b10e1a5a8/doc/conf.py
    - https://github.com/chainer/chainer/pull/2758/
    """
    if domain != "py" or not info["module"]:
        return None

    # tag
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip()
    except (subprocess.CalledProcessError, OSError):
        print("Failed to execute git to get revision")
        return None
    revision = revision.decode("utf-8")

    obj = sys.modules.get(info["module"])
    if obj is None:
        return None
    for comp in info["fullname"].split("."):
        obj = getattr(obj, comp)

    # filename
    try:
        filename = inspect.getsourcefile(obj)
    except Exception:
        return None
    if filename is None:
        return None

    # relpath
    pkg_root_dir = os.path.dirname(__import__(MODULE).__file__)
    filename = os.path.realpath(filename)
    if not filename.startswith(pkg_root_dir):
        return None
    relpath = os.path.relpath(filename, pkg_root_dir)

    # line number
    try:
        linenum = inspect.getsourcelines(obj)[1]
    except Exception:
        linenum = ""

    return "https://github.com/{}/{}/blob/{}/{}/{}#L{}".format(
        GH_ORGANIZATION, GH_PROJECT, revision, MODULE, relpath, linenum
    )
