#!/usr/bin/env python3
"""Build .ipynb notebooks from jupytext-style .py sources in notebooks_src/.

Convention
----------
    # %%              starts a code cell
    # %% [markdown]   starts a markdown cell, whose lines each start with "# "

Keeping the sources as plain Python makes them easy to review and diff;
this script only converts them to the JSON notebook format.
"""

from __future__ import annotations

import re
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "notebooks_src"
DST = ROOT / "notebooks"

CELL_RE = re.compile(r"^# %%(?:\s*\[(markdown)\]|\s*code)?\s*$")


def parse(text: str) -> list[tuple[str, str]]:
    cells: list[tuple[str, str]] = []
    kind: str | None = None
    lines: list[str] = []
    for line in text.splitlines():
        match = CELL_RE.match(line)
        if match:
            if kind is not None:
                cells.append((kind, "\n".join(lines).strip("\n")))
            kind = match.group(1) or "code"
            lines = []
            continue
        if kind is None:
            continue
        lines.append(re.sub(r"^# ?", "", line) if kind == "markdown" else line)
    if kind is not None:
        cells.append((kind, "\n".join(lines).strip("\n")))
    return cells


def build(path: Path) -> Path:
    notebook = nbf.v4.new_notebook()
    for kind, text in parse(path.read_text()):
        if kind == "markdown":
            notebook.cells.append(nbf.v4.new_markdown_cell(text))
        else:
            notebook.cells.append(nbf.v4.new_code_cell(text))
    notebook.metadata.update(
        {
            "kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
            "language_info": {"name": "python"},
        }
    )
    out = DST / (path.stem + ".ipynb")
    nbf.write(notebook, str(out))
    return out


def main() -> None:
    DST.mkdir(exist_ok=True)
    sources = sorted(SRC.glob("*.py"))
    if not sources:
        raise SystemExit(f"no notebook sources in {SRC}")
    for source in sources:
        print("wrote", build(source))


if __name__ == "__main__":
    main()
