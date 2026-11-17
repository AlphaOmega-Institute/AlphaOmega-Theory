"""Check the editorial ledger against the source inventory and the consolidated TeX."""
import csv
import re
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MANUSCRIPTS = ROOT / "manuscripts"
SOURCES = {
    "LG": MANUSCRIPTS / "lg/LG_Master_Unified.tex",
    "TAUM": MANUSCRIPTS / "taum/Toward_a_Universal_Model_Self_Contained.tex",
}
LABEL = re.compile(r"\\label\{([^}]+)\}")
FIELDS = ["source", "source_line", "source_kind", "source_key", "clause", "coverage",
          "draft_labels", "planned_homes", "evidence_or_owed"]


def check():
    with (HERE / "statements.csv").open(encoding="utf-8", newline="") as stream:
        inventory = {(r["manuscript"], r["label"] or "@" + r["line"]): r
                     for r in csv.DictReader(stream)}
    labels = {}
    for source, path in SOURCES.items():
        labels[source] = {label: line_num
                          for line_num, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
                          for label in LABEL.findall(line)}
    draft = (MANUSCRIPTS / "consolidated/Universal_Model.tex").read_text(encoding="utf-8")
    destinations = set(LABEL.findall(draft))
    represented = set()
    supplements = set()
    clauses = set()
    counts = Counter()
    with (HERE / "migration.csv").open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != FIELDS:
            raise ValueError("Unexpected ledger columns")
        for row in reader:
            if set(row) != set(FIELDS) or any(row[k] is None or not row[k].strip()
                                              for k in FIELDS if k != "draft_labels"):
                raise ValueError(f"Incomplete ledger row at line {reader.line_num}")
            source, key = row["source"], row["source_key"]
            identity = (source, key)
            location = f"{source} {key} at ledger line {reader.line_num}"
            if source not in SOURCES:
                raise ValueError(f"Unknown source: {location}")
            source_line = int(row["source_line"])
            if row["source_kind"] == "supplement":
                if identity in inventory or labels[source].get(key) != source_line:
                    raise ValueError(f"Incorrect supplemental source anchor: {location}")
                supplements.add(identity)
            else:
                entry = inventory.get(identity)
                if not entry or (entry["kind"], int(entry["line"])) != (row["source_kind"], source_line):
                    raise ValueError(f"Source inventory mismatch: {location}")
                represented.add(identity)
            clause = (*identity, row["clause"])
            if clause in clauses:
                raise ValueError(f"Duplicate clause: {location}")
            clauses.add(clause)
            coverage = row["coverage"]
            targets = set(row["draft_labels"].split())
            if coverage not in ("explicit", "recoverable", "owed"):
                raise ValueError(f"Unknown coverage: {location}")
            if (coverage == "owed") != (not targets) or targets - destinations:
                raise ValueError(f"Missing or invalid draft anchor: {location}")
            counts[coverage] += 1
    missing = sorted(inventory.keys() - represented)
    if missing:
        raise ValueError(f"Missing inventory entries: {missing}")
    print(f"{len(represented)} formal source entries, {len(supplements)} supplemental anchors, "
          f"{len(clauses)} clause rows")
    print(", ".join(f"{k}: {counts[k]}" for k in ("explicit", "recoverable", "owed")))


if __name__ == "__main__":
    check()
