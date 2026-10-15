# Configuration file for the Sphinx documentation builder.
#
# The package is installed (``uv sync --group docs``), so no sys.path tweaks
# are needed. Build with:
#
#   uv run --group docs sphinx-build -W --keep-going -b html docs docs/_build/html

import datetime
from importlib.metadata import version as _pkg_version

# -- Project information -----------------------------------------------

project = "orgdisord"
author = "Kane Shenton"
copyright = f"{datetime.date.today().year}, {author}"

release = _pkg_version("orgdisord")
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx_click",
    "sphinx_copybutton",
]

root_doc = "index"
source_suffix = ".rst"
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
language = "en"
pygments_style = "sphinx"

# autodoc / autosummary
autodoc_class_signature = "separated"
autodoc_member_order = "bysource"
autodoc_default_options = {
    "show-inheritance": True,
}
autoclass_content = "both"
autosummary_generate = True

# napoleon (Google and NumPy docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_include_special_with_doc = True


# -- Options for HTML output -------------------------------------------

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/CCP-NC/orgdisord",
    "repository_branch": "main",
    "path_to_docs": "docs",
    "use_repository_button": True,
    "use_edit_page_button": True,
}

# CopyButton configuration
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True
copybutton_line_continuation_character = "\\"
copybutton_here_doc_delimiter = "EOT"
copybutton_selector = "div:not(.no-copybutton) > div.highlight > pre"
