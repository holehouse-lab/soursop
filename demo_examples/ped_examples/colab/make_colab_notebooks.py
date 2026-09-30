"""
Build the Google Colab versions of the PED example notebooks.

The notebooks in ``demo_examples/ped_examples/`` are the originals. This
script copies each one into this folder and makes the changes needed to run
it in Colab:

* adds an "Open in Colab" badge, and a link back to the original notebook
  (which is saved with its output);
* adds a cell that installs SOURSOP;
* lets the figures fall back to Liberation Sans or DejaVu Sans, since Colab
  does not have Arial;
* points references to other notebooks at their Colab versions;
* removes all output.

Edit the original notebooks, then re-run this script to update the Colab
versions::

    python make_colab_notebooks.py
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell

REPO = "holehouse-lab/soursop"
BRANCH = "master"
FOLDER = "demo_examples/ped_examples"
BADGE = "https://colab.research.google.com/assets/colab-badge.svg"

HERE = Path(__file__).resolve().parent
ORIGINALS = HERE.parent

INSTALL = '''\
# Colab does not come with SOURSOP, so install it first. This uses uv, which
# is much faster than pip (and is itself installed first if it is missing);
# --system installs into Colab's own Python rather than a virtual environment.
# The PED module is new in SOURSOP 2.0.8; if that version is not on PyPI yet,
# this installs the latest version from GitHub instead.
!command -v uv >/dev/null || pip install -q uv
!uv pip install --system -q "soursop>=2.0.8" 2>/dev/null || uv pip install --system -q "git+https://github.com/holehouse-lab/soursop.git"'''

# (original text, Colab text) pairs applied to the set-up cell of every notebook
SETUP_EDITS = [
    (
        """\
# Figure style: Arial throughout, text kept editable in Illustrator, no grid
# lines, and the xmgrace colour order (black, red, blue, green, purple)
plt.rcParams.update({
    "font.family": "Arial",""",
        """\
# Figure style: Arial throughout, text kept editable in Illustrator, no grid
# lines, and the xmgrace colour order (black, red, blue, green, purple).
# Colab does not have Arial, so there figures fall back to Liberation Sans
# (metrically identical to Arial) or matplotlib's default, DejaVu Sans
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],""",
    ),
    (
        """\
# Downloaded ensembles are kept here and reused, so each is only fetched once
# (shared by all the notebooks in this folder)""",
        """\
# Downloaded ensembles are kept here and reused, so each is only fetched once
# per session (Colab deletes them when the session ends)""",
    ),
]

# notebook 7 saves files to disk, which in Colab means the Colab machine
SAVE_TEXT = (
    "The file can be read back in with `SSTrajectory` like any other PDB "
    "file, or opened in PyMOL, VMD and so on."
)
SAVE_NOTE = (
    "\n\nIn Colab, files are saved on the Colab machine rather than your "
    "computer. Download them from the file browser (the folder icon on the "
    "left), or with `google.colab.files.download(path)`."
)


def colab_url(name: str) -> str:
    """Return the URL that opens a Colab notebook from GitHub.

    Parameters
    ----------
    name : str
        Notebook file name, e.g. ``"01_getting_started.ipynb"``.

    Returns
    -------
    str
        The ``colab.research.google.com`` URL for the notebook in this folder
        on the ``master`` branch.
    """
    return f"https://colab.research.google.com/github/{REPO}/blob/{BRANCH}/{FOLDER}/colab/{name}"


def github_url(name: str) -> str:
    """Return the GitHub URL of an original (local) notebook.

    Parameters
    ----------
    name : str
        Notebook file name, e.g. ``"01_getting_started.ipynb"``.

    Returns
    -------
    str
        The ``github.com`` URL for the original notebook, which is saved with
        its output.
    """
    return f"https://github.com/{REPO}/blob/{BRANCH}/{FOLDER}/{name}"


def replace_once(text: str, old: str, new: str, label: str) -> str:
    """Replace ``old`` with ``new`` in ``text``, insisting it appears once.

    This makes the script fail loudly if an original notebook changes in a
    way the Colab edits no longer match, rather than silently skipping it.

    Parameters
    ----------
    text : str
        Text to edit.
    old : str
        Text to replace. Must appear exactly once in ``text``.
    new : str
        Replacement text.
    label : str
        Name used in the error message (e.g. the notebook name).

    Returns
    -------
    str
        The edited text.

    Raises
    ------
    ValueError
        If ``old`` does not appear exactly once.
    """
    count = text.count(old)
    if count != 1:
        raise ValueError(
            f"{label}: expected to find this text once, found it {count} times:\n{old}"
        )
    return text.replace(old, new)


def make_colab_notebook(path: Path, numbered: Dict[str, str]) -> Any:
    """Build the Colab version of one original notebook.

    Parameters
    ----------
    path : Path
        The original notebook.
    numbered : dict
        Maps each notebook's number (``"3"``) to its file name, so that
        references such as "notebook 3" can be linked to the Colab version.

    Returns
    -------
    nbformat.NotebookNode
        The Colab notebook, ready to write.
    """
    name = path.name
    nb = nbformat.read(path, as_version=4)  # type: ignore[no-untyped-call]

    for cell in nb.cells:
        cell.metadata.pop("execution", None)
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
        else:
            cell.source = re.sub(
                r"\bnotebook (\d)\b",
                lambda m: f"[notebook {m.group(1)}]({colab_url(numbered[m.group(1)])})",
                cell.source,
            )

    # cell 0 is the title and introduction, cell 1 the imports and set-up
    setup = nb.cells[1]
    for old, new in SETUP_EDITS:
        setup.source = replace_once(setup.source, old, new, name)

    for cell in nb.cells:
        if cell.cell_type == "markdown" and SAVE_TEXT in cell.source:
            cell.source = cell.source.replace(SAVE_TEXT, SAVE_TEXT + SAVE_NOTE)

    # fixed ids keep the output identical from one run of this script to the next
    badge = new_markdown_cell(  # type: ignore[no-untyped-call]
        f'<a href="{colab_url(name)}" target="_parent"><img src="{BADGE}" alt="Open In Colab"/></a>\n\n'
        f"This is the Google Colab version of this notebook. A copy with all its output "
        f"is [in the SOURSOP repository]({github_url(name)}).",
        id="colab-badge",
    )
    install = new_code_cell(INSTALL, id="colab-install")  # type: ignore[no-untyped-call]

    nb.cells = [badge, nb.cells[0], install] + nb.cells[1:]
    nb.metadata = {
        "colab": {"provenance": []},
        "kernelspec": {"name": "python3", "display_name": "Python 3"},
        "language_info": {"name": "python"},
    }
    nbformat.validate(nb)
    return nb


def main() -> None:
    """Write a Colab version of every original notebook into this folder."""
    originals: List[Path] = sorted(ORIGINALS.glob("0*.ipynb"))
    numbered = {str(int(path.name[:2])): path.name for path in originals}
    for path in originals:
        nb = make_colab_notebook(path, numbered)
        nbformat.write(nb, HERE / path.name)  # type: ignore[no-untyped-call]
        print(f"wrote colab/{path.name}")


if __name__ == "__main__":
    main()
