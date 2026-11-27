"""Execute the notebook headlessly and rebuild all figures."""
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
from IPython.display import Image, Markdown


def print_output(item):
    if isinstance(item, Image):
        return
    print(item.data if isinstance(item, Markdown) else item)


def run_notebook(path):
    notebook = json.loads(path.read_text(encoding="utf-8"))
    scope = {"__name__": "__main__"}
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        # Replace display in the first cell after it imports IPython's renderer.
        exec(compile(source, f"{path.name}:{cell['id']}", "exec"), scope)
        scope["display"] = print_output
    return scope


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    os.chdir(root)
    run_notebook(root / "flow_battery_benchmark.ipynb")
    print("\nFigures saved in figures/.")
