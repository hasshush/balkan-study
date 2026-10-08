#!/usr/bin/env python3
"""Replace machine-specific paths in executed notebook outputs.

Notebook sources print absolute paths (project root, plink command lines).
This keeps committed notebooks portable by rewriting them to repo-relative
paths after execution.
"""

from pathlib import Path

ROOT = str(Path.cwd().resolve())
REPLACEMENTS = [
    (ROOT, "."),
    (str(Path.home()), "~"),
]


def main() -> None:
    for path in sorted(Path("notebooks").glob("*.ipynb")):
        text = path.read_text()
        sanitized = text
        for old, new in REPLACEMENTS:
            sanitized = sanitized.replace(old, new)
        if sanitized != text:
            path.write_text(sanitized)
            print(f"sanitized {path}")
        else:
            print(f"clean     {path}")


if __name__ == "__main__":
    main()
