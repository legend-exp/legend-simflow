# Configuration file for the Sphinx documentation builder.
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, Path(__file__).parents[2].resolve().as_posix())

project = "legend-simflow"
copyright = "2025, the LEGEND Collaboration"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "sphinx_inline_tabs",
    "myst_parser",
]

source_suffix = [".rst", ".md"]
exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "Thumbs.db",
    ".DS_Store",
    ".env",
    ".venv",
]

# Furo theme
html_theme = "furo"
html_theme_options = {
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/legend-exp/legend-simflow",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
    "source_repository": "https://github.com/legend-exp/legend-simflow",
    "source_branch": "main",
    "source_directory": "docs/source",
}

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
]

# sphinx-napoleon
# enforce consistent usage of NumPy-style docstrings
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_ivar = True
napoleon_use_rtype = False

# intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "awkward": ("https://awkward-array.org/doc/stable", None),
    "numba": ("https://numba.readthedocs.io/en/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "h5py": ("https://docs.h5py.org/en/stable", None),
    "pint": ("https://pint.readthedocs.io/en/stable", None),
    "hist": ("https://hist.readthedocs.io/en/latest", None),
    "dspeed": ("https://dspeed.readthedocs.io/en/stable", None),
    "daq2lh5": ("https://legend-daq2lh5.readthedocs.io/en/stable", None),
    "lgdo": ("https://legend-pydataobj.readthedocs.io/en/stable", None),
    "dbetto": ("https://dbetto.readthedocs.io/en/stable", None),
    "pylegendmeta": ("https://pylegendmeta.readthedocs.io/en/stable", None),
    "reboost": ("https://reboost.readthedocs.io/en/stable", None),
}  # add new intersphinx mappings here

# sphinx-autodoc
autodoc_default_options = {"ignore-module-all": True}
# Include __init__() docstring in class docstring
autoclass_content = "both"
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented_params"
autodoc_typehints_format = "short"
