# Sphinx configuration for statphys-viz documentation
import statphys_viz

project = "statphys-viz"
author = "Pedro Sordo Martínez"
copyright = "2026, Pedro Sordo Martínez"
release = statphys_viz.__version__
version = release
language = "es"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "nbsphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = "statphys-viz Documentation"
html_short_title = "statphys-viz"

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

napoleon_google_docstring = True
napoleon_numpy_docstring = True
autosummary_generate = True
