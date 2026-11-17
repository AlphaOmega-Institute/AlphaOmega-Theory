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
STATUS = {"absent", "definition", "hypothesis", "model choice", "identification",
          "proved result", "conjecture", "open question", "worked example"}
COVERAGE = ("explicit", "recoverable", "owed")
FIELDS = ["source", "source_line", "source_kind", "source_key", "clause", "coverage",
          "draft_status", "draft_labels", "draft_quote", "planned_homes", "evidence_or_owed"]
OPTIONAL = {"draft_labels", "draft_quote"}


def labels_by_line(text):
    return {label: line_num for line_num, line in enumerate(text.splitlines(), 1)
            for label in LABEL.findall(line)}


def check():
    with (HERE / "statements.csv").open(encoding="utf-8", newline="") as stream:
        inventory = {(r["manuscript"], r["label"] or "@" + r["line"]): r
                     for r in csv.DictReader(stream)}
    labels = {source: labels_by_line(path.read_text(encoding="utf-8"))
              for source, path in SOURCES.items()}
    draft = (MANUSCRIPTS / "consolidated/Universal_Model.tex").read_text(encoding="utf-8")
    destinations = labels_by_line(draft)
    draft_text = " ".join(draft.split())
    represented = set()
    supplements = set()
    extensions = set()
    clauses = set()
    counts = Counter()
    with (HERE / "migration.csv").open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != FIELDS:
            raise ValueError("Unexpected ledger columns")
        for row in reader:
            if (set(row) != set(FIELDS) or None in row.values()
                    or any(not row[k].strip() for k in FIELDS if k not in OPTIONAL)):
                raise ValueError(f"Incomplete ledger row at line {reader.line_num}")
            source, key = row["source"], row["source_key"]
            identity = (source, key)
            location = f"{source} {key} at ledger line {reader.line_num}"
            source_line = int(row["source_line"])
            targets = set(row["draft_labels"].split())
            coverage = row["coverage"]
            if source == "UM":
                if (row["source_kind"] != "draft extension" or coverage != "explicit"
                        or destinations.get(key) != source_line or key not in targets):
                    raise ValueError(f"Incorrect draft-only provenance: {location}")
                extensions.add(identity)
            elif source not in SOURCES:
                raise ValueError(f"Unknown source: {location}")
            elif row["source_kind"] == "supplement":
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
            if coverage not in COVERAGE:
                raise ValueError(f"Unknown coverage: {location}")
            roles = [role.strip() for role in row["draft_status"].split("|")]
            if (set(roles) - STATUS or len(set(roles)) != len(roles)
                    or ("absent" in roles and roles != ["absent"])
                    or (coverage == "owed") != (roles == ["absent"])):
                raise ValueError(f"Invalid draft status: {location}")
            if (coverage == "owed") != (not targets) or targets - destinations.keys():
                raise ValueError(f"Missing or invalid draft anchor: {location}")
            quote = " ".join(row["draft_quote"].split())
            if coverage == "owed" and quote:
                raise ValueError(f"Owed clause has a draft quote: {location}")
            if (coverage == "explicit" and not quote
                    and all(target.startswith(("part:", "sec:", "subsec:")) for target in targets)):
                raise ValueError(f"Broad draft anchor needs a unique quote: {location}")
            if quote:
                first = draft_text.find(quote)
                if first < 0 or first != draft_text.rfind(quote):
                    raise ValueError(f"Draft quote must occur exactly once: {location}")
            counts[coverage] += 1
    missing = sorted(inventory.keys() - represented)
    if missing:
        raise ValueError(f"Missing inventory entries: {missing}")
    print(f"{len(represented)} formal source entries, {len(supplements)} supplemental source anchors, "
          f"draft-only anchors: {len(extensions)}, {len(clauses)} clause rows")
    print(", ".join(f"{k}: {counts[k]}" for k in COVERAGE))


if __name__ == "__main__":
    check()
