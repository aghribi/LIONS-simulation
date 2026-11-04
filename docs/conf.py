project = "LIONS"
author = "Omar Kamalou, Adnan Ghribi"
copyright = "2026, Omar Kamalou, Adnan Ghribi"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build"]

myst_heading_anchors = 3
autodoc_member_order = "bysource"
autodoc_default_options = {"members": True, "undoc-members": True}

html_theme = "furo"
html_title = "LIONS"


def _skip_qt_signals(app, what, name, obj, skip, options):
    # PyQt signal docstrings are bare signatures ("*args") that docutils
    # rejects.
    from PyQt5.QtCore import pyqtSignal

    return True if isinstance(obj, pyqtSignal) else skip


def setup(app):
    app.connect("autodoc-skip-member", _skip_qt_signals)
