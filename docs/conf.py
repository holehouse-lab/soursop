# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/stable/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
#import os
#import sys
#sys.path.insert(0, os.path.abspath('.'))

import os
import sys
sys.path.insert(0, os.path.abspath('..'))  # Adjust the path as needed

# -- Project information -----------------------------------------------------

project = 'soursop'
copyright = ("2015-2026, Alex Holehouse, Jared Lalmansingh [a MOLSSI sponsored project]")
author = 'Alex Holehouse'

# The version is read automatically so the docs always build with the current
# soursop version rather than a hardcoded string. We try, in order:
#   1. versioningit computed straight from the git tags. This works on Read the
#      Docs (versioningit is a docs dependency) even though RTD never installs
#      the soursop package, and it ignores any stale installed distribution.
#   2. the versioningit-written soursop/_version.py, for built/installed trees
#      (e.g. an sdist) that have no .git directory. NOTE this file is gitignored,
#      so it is absent on a fresh clone - hence versioningit is tried first.
#   3. the installed package metadata.
#   4. "unknown".
import re

_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")


def _get_soursop_version():
    # 1. versioningit from git
    try:
        import versioningit

        _v = versioningit.get_version(project_dir=_REPO_ROOT)
        # reject the pyproject default-version ("1+unknown") used when git/tags
        # cannot be resolved (e.g. a too-shallow clone with no reachable tag).
        if _v and "unknown" not in _v:
            return _v
    except Exception:
        pass

    # 2. versioningit-written _version.py
    try:
        with open(os.path.join(_REPO_ROOT, "soursop", "_version.py")) as _fh:
            _m = re.search(r"""__version__\s*=\s*['"]([^'"]+)['"]""", _fh.read())
            if _m:
                return _m.group(1)
    except OSError:
        pass

    # 3. installed package metadata
    try:
        from importlib.metadata import version as _dist_version, PackageNotFoundError

        try:
            return _dist_version("soursop")
        except PackageNotFoundError:
            pass
    except Exception:
        pass

    return "unknown"


def _get_release_date(rel):
    """Month/Year the given version was released, from CHANGELOG.md.

    The changelog headers carry the release month (e.g. ``## 2.0.3 (July
    2026)``). We look up the entry matching ``rel`` and, failing that, fall
    back to the most recent versioned entry. Returns "" if nothing is found.
    """
    if rel == "unknown":
        return ""
    changelog = os.path.join(_REPO_ROOT, "CHANGELOG.md")
    try:
        with open(changelog) as _fh:
            _text = _fh.read()
    except OSError:
        return ""
    # exact match: "## <rel> (<Month Year>)"
    _m = re.search(
        r"^##\s+" + re.escape(rel) + r"\s+\(([^)]+)\)", _text, re.MULTILINE
    )
    if _m:
        return _m.group(1).strip()
    # fall back to the first "## X.Y.Z (<Month Year>)" header
    _m = re.search(r"^##\s+\d[\w.]*\s+\(([^)]+)\)", _text, re.MULTILINE)
    return _m.group(1).strip() if _m else ""


# The full version, including alpha/beta/rc tags (PEP 440 local segment, e.g.
# "+1.gabc123", is dropped for display); the short X.Y version is derived from it.
release = _get_soursop_version().split("+")[0]
version = ".".join(release.split(".")[:2]) if release != "unknown" else release

# The release Month/Year (from CHANGELOG.md). Exposed to .rst as the
# |version_info| substitution below (version, optionally with the date).
release_date = _get_release_date(release)
if release_date:
    _version_info = f"{release} (released {release_date})"
else:
    _version_info = release
rst_prolog = f".. |version_info| replace:: {_version_info}\n"


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# See: https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html
autodoc_default_options = {
    # 'members': True,
    #'undoc-members': True,
    'private-members': True,
    #'special-members': True,
    #'inherited-members': True,
    #'show-inheritance': True
}

# added july 2024 - stubs for imports
autodoc_mock_imports = ["mdtraj", "threadpoolctl", "natsort", "matplotlib", "pandas"]


autosummary_generate = True

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
#    'numpydoc', # automatically includes `sphinx.ext.autosummary`
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',    
]


#mathjax_config = {
#    'extensions': ['tex2jax.js'],
#    'jax': ['input/TeX', 'output/HTML-CSS'],
#}



# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'soursopdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'soursop.tex', 'soursop Documentation',
     'soursop', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'soursop', 'soursop Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'soursop', 'soursop Documentation',
     author, 'soursop', 'Analysis package for all-atom simulations of proteins, with a specific focus on intrinsically disordered proteins.',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------
