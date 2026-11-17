"""Shared loading for the ledger audit scripts."""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
SOURCES = {
    "LG": ROOT / "manuscripts/lg/LG_Master_Unified.tex",
    "TAUM": ROOT / "manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex",
}
THEOREM_LEVEL = ("theorem", "proposition", "corollary", "lemma")


def read_csv(name):
    with (INVENTORY / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def source_lines():
    return {name: path.read_text(encoding="utf-8").splitlines() for name, path in SOURCES.items()}


def entry_key(entry):
    """Key used by the ledger: the label, or '@line' for an unlabelled environment."""
    return entry["manuscript"], entry["label"] or "@" + entry["line"]


def statement(lines, entry):
    """Lines of the environment that begins at the entry's line."""
    start, kind = int(entry["line"]), entry["kind"]
    out = []
    for i in range(start - 1, min(start + 150, len(lines))):
        out.append(lines[i])
        if i > start - 1 and re.search(r"\\end\{" + re.escape(kind) + r"\}", lines[i]):
            break
    return out


def rows_by_entry(ledger):
    grouped = {}
    for row in ledger:
        grouped.setdefault((row["source"], row["source_key"]), []).append(row)
    return grouped
