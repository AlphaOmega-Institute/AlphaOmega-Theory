#!/usr/bin/env python3
"""Write and check the proposed ledger rows for the symmetry-averaging fragment.

    python consolidation/inventory/proposed/build_symmetry_average_rows.py

The rows have the columns of migration.csv. They are proposals: `draft_labels` and
`draft_quote` refer to manuscripts/consolidated/fragments/symmetry_average.tex, which is
not part of Universal_Model.tex, so they cannot enter migration.csv until the text has a
permanent home. The script takes the current rows for the clauses named below, wherever
they sit in the source, applies the updates, appends a row for a source passage that had
none, and applies the validator's rules with the fragment as the draft passage.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
NAME = "symmetry_average"
FRAGMENT = ROOT / f"manuscripts/consolidated/fragments/{NAME}.tex"
DRAFT = ROOT / "manuscripts/consolidated/Universal_Model.tex"
TAUM = ROOT / "manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex"
LG = ROOT / "manuscripts/lg/LG_Master_Unified.tex"
STATUS = {"absent", "definition", "hypothesis", "model choice", "identification",
          "proved result", "conjecture", "open question", "worked example"}
FIELDS = ["source", "source_line", "source_kind", "source_key", "clause", "coverage",
          "draft_status", "draft_labels", "draft_quote", "planned_homes", "evidence_or_owed"]

# (source_key, clause) -> (draft_status, draft_labels, draft_quote, evidence)
UPDATES = {
    ("thm:Haar-projector", "Haar average and its invariant range"):
        ("proved result", "prop:symmetry-projection eq:symmetry-average eq:symmetry-absorption", "",
         "The average is defined by its matrix elements, which gives the operator of TAUM's vector integral. Boundedness, idempotence, self-adjointness, the invariant range and absorption of the representation on both sides are proved."),
    ("thm:Haar-projector", "Compatibility with the grading"):
        ("proved result", "prop:symmetry-compatibility eq:symmetry-sector-sum", "",
         "Part (a): when the representation commutes with the sector symmetry, the invariant projection commutes with it, is the direct sum of the invariant projections of the two sectors, and equals its adjoint for the indefinite form."),
    ("cor:finite-group-average", "Finite average and trivial cubic representation"):
        ("proved result", "eq:symmetry-finite-average subsec:symmetry-projection", "For a finite group of rotations of a cube",
         "The normalized finite sum is part of the proposition. The sentence after its proof identifies the average over a group of cube rotations with the projection onto the vectors on which that group acts trivially."),
    ("@1177", "Distinguishing the four projections and the two symmetry groups"):
        ("proved result|worked example", "subsec:symmetry-compatibility subsec:symmetry-pair eq:symmetry-pair-exchange",
         "Its construction does not identify it with any of them.",
         "The invariant projection is separated from the sector projection, the support projections and a prepared or moving projection. In the pair model it equals the prepared projection at one setting, commutes with it at another and does neither in between. The difference between a finite and a continuous group is shown by the cube rotations on quartic polynomials and by the cyclic subgroups of the orbital rotations."),
    ("thm:architecture", "Conditional compact symmetry and the reducing Haar projection"):
        ("proved result", "prop:symmetry-projection prop:symmetry-compatibility", "",
         "For a compact group acting by unitaries that commute with the sector symmetry, the average is the projection onto the invariant vectors and it commutes with that symmetry."),
}

# Source passages that had no row: (source_key, source_line, clause, draft_status, draft_labels, draft_quote, planned_homes, evidence)
ADDITIONS = [
    ("sec:symmetry-averaging", 1079, "The symmetry average is a separate object from the sector projections and a moving projection, with its own notation",
     "definition", "sec:symmetry-average eq:symmetry-average",
     "", "III.2",
     "The opening of TAUM's section. The fragment defines the projection with its own symbol and states the distinction after the compatibility proposition."),
]


def labels(path):
    return set(re.findall(r"\\label\{([^}]+)\}", path.read_text(encoding="utf-8")))


def label_lines(path):
    found = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        for key in re.findall(r"\\label\{([^}]+)\}", line):
            found[key] = number
    return found


def main():
    with (INVENTORY / "migration.csv").open(encoding="utf-8", newline="") as stream:
        ledger = list(csv.DictReader(stream))
    current = [r for r in ledger if r["source"] == "TAUM" and (r["source_key"], r["clause"]) in UPDATES]
    assert len(current) == len(UPDATES), sorted(set(UPDATES) - {(r["source_key"], r["clause"]) for r in current})
    # every clause of the ported section is among them
    section = [r for r in ledger if r["source"] == "TAUM" and 1078 <= int(r["source_line"]) <= 1187]
    assert all((r["source_key"], r["clause"]) in UPDATES for r in section), "a clause of the section has no update"
    fragment_labels, draft_labels = labels(FRAGMENT), labels(DRAFT)
    assert not fragment_labels & draft_labels, fragment_labels & draft_labels
    for source in (LG, TAUM):
        assert not fragment_labels & labels(source), fragment_labels & labels(source)
    fragment_text = " ".join(FRAGMENT.read_text(encoding="utf-8").split())
    source_lines = label_lines(TAUM)
    rows = []
    for row in current:
        status, targets, quote, evidence = UPDATES[(row["source_key"], row["clause"])]
        rows.append({**row, "coverage": "explicit", "draft_status": status, "draft_labels": targets,
                     "draft_quote": quote, "evidence_or_owed": evidence})
    existing_keys = {r["source_key"] for r in ledger if r["source"] == "TAUM"}
    for key, line, clause, status, targets, quote, homes, evidence in ADDITIONS:
        assert key not in existing_keys and source_lines.get(key) == line, (key, source_lines.get(key))
        rows.append({"source": "TAUM", "source_line": str(line), "source_kind": "supplement", "source_key": key,
                     "clause": clause, "coverage": "explicit", "draft_status": status, "draft_labels": targets,
                     "draft_quote": quote, "planned_homes": homes, "evidence_or_owed": evidence})
    rows.sort(key=lambda r: int(r["source_line"]))
    seen = set()
    for row in rows:
        clause = (row["source_key"], row["clause"])
        assert clause not in seen, clause
        seen.add(clause)
        roles = row["draft_status"].split("|")
        assert not set(roles) - STATUS and "absent" not in roles, row["draft_status"]
        targets = row["draft_labels"].split()
        assert targets and all(t in fragment_labels or t in draft_labels for t in targets), targets
        quote = " ".join(row["draft_quote"].split())
        if all(t.startswith(("part:", "sec:", "subsec:")) for t in targets):
            assert quote, clause
        if quote:
            assert fragment_text.count(quote) == 1, (fragment_text.count(quote), quote)
    with (HERE / f"{NAME}_rows.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    before = {c: sum(r["coverage"] == c for r in current) for c in ("explicit", "recoverable", "owed")}
    total = {c: sum(r["coverage"] == c for r in ledger) for c in ("explicit", "recoverable", "owed")}
    after = {"explicit": total["explicit"] + before["recoverable"] + before["owed"] + len(ADDITIONS),
             "recoverable": total["recoverable"] - before["recoverable"], "owed": total["owed"] - before["owed"]}
    print(f"{len(rows)} proposed rows: {len(UPDATES)} updated, {len(ADDITIONS)} for a passage without a row")
    print(f"updated clauses before: {before}; proposed: all explicit")
    print(f"ledger now: {len(ledger)} rows {total}; if accepted: {len(ledger) + len(ADDITIONS)} rows {after}")
    print(f"fragment labels: {len(fragment_labels)}, none shared with the manuscript, LG or TAUM")


if __name__ == "__main__":
    main()
