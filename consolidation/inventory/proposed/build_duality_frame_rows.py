#!/usr/bin/env python3
"""Write and check the proposed ledger rows for the duality-frame fragment.

    python consolidation/inventory/proposed/build_duality_frame_rows.py

The rows have the columns of ../migration.csv. They are proposals: `draft_labels` and
`draft_quote` refer to manuscripts/consolidated/fragments/duality_frame.tex, which is
not part of Universal_Model.tex, so ../check_migration.py would reject them. They wait
here until the text is integrated and has its final labels.

The script takes the current row for each clause named below, checks that it is not yet
explicit, replaces its coverage, status, labels, quote and evidence, writes
duality_frame_rows.csv, and applies the validator's rules with the fragment as the
draft passage. No explicit row is touched and no row is added.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
NAME = "duality_frame"
FRAGMENT = ROOT / f"manuscripts/consolidated/fragments/{NAME}.tex"
DRAFT = ROOT / "manuscripts/consolidated/Universal_Model.tex"
SOURCES = {"TAUM": ROOT / "manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex",
           "LG": ROOT / "manuscripts/lg/LG_Master_Unified.tex"}
STATUS = {"absent", "definition", "hypothesis", "model choice", "identification",
          "proved result", "conjecture", "open question", "worked example"}
FIELDS = ["source", "source_line", "source_kind", "source_key", "clause", "coverage",
          "draft_status", "draft_labels", "draft_quote", "planned_homes", "evidence_or_owed"]

# (source, source_key, clause) -> (draft_status, draft_labels, draft_quote, evidence)
UPDATES = {
    ("LG", "def:derived-duality-frame", "State-indexed frame for a differentiated normalized state"):
        ("definition", "eq:duality-frame",
         r"A sharp frame is a pair of orthogonal projections $P_+$ and $P_-$ on $\mathcal H_\Phi$ with",
         "The two orthogonal projections, their sum and the support bound in the complement of the ground are displayed, with the grading and its two properties. LG's reading of the projections as aspects and the dependence of the frame on profile and state are attributed to the cognitive model."),
    ("LG", "thm:duality-frame-vanishes-ground", "Vanishing of both aspect projections and their sum on the luminous ground"):
        ("proved result", "prop:duality-frame-ground eq:duality-frame-norms", "",
         "Proved from the support bound. The text adds the splitting of the differentiated part into the two frame components and a part the frame misses, the remark that a zero mean of the grading does not characterize the ground, and the same vanishing for unsharp effects."),
    ("LG", "choice:duality-generator", "Optional spectral construction of a duality frame"):
        ("definition|proved result", "subsec:duality-frame",
         r"One way to obtain a frame is spectral.",
         "A bounded self-adjoint operator that annihilates the ground gives a frame through its positive and negative spectral projections, with the grading as its sign. Every sharp frame arises in this way from its own grading. Assigning such an operator to each profile is attributed to LG as a choice for a realization."),
    ("LG", "q:projection-or-effect-duality-frame", "Choice of representation for the frame"):
        ("open question", "subsec:duality-frame prop:duality-frame-ground",
         r"Whether a realization should use sharp projections, unsharp effects, a state-dependent tensor factorization or a more general operator-valued decomposition is open.",
         "The four alternatives are stated as open. The vanishing at the ground is shown to hold for effects as well, so it does not depend on the choice of sharp projections."),
}

# (source, source_key, line, clause, planned_homes, draft_status, draft_labels, draft_quote, evidence)
ADDITIONS = []


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
    fragment_labels, draft_labels = labels(FRAGMENT), labels(DRAFT)
    assert not fragment_labels & draft_labels, fragment_labels & draft_labels
    for path in SOURCES.values():
        assert not fragment_labels & labels(path), fragment_labels & labels(path)
    fragment_text = " ".join(FRAGMENT.read_text(encoding="utf-8").split())
    rows, before = [], {"explicit": 0, "recoverable": 0, "owed": 0}
    for row in ledger:
        key = (row["source"], row["source_key"], row["clause"])
        if key not in UPDATES:
            continue
        assert row["coverage"] != "explicit", key
        before[row["coverage"]] += 1
        status, targets, quote, evidence = UPDATES[key]
        rows.append({**row, "coverage": "explicit", "draft_status": status, "draft_labels": targets,
                     "draft_quote": quote, "evidence_or_owed": evidence})
    assert len(rows) == len(UPDATES), set(UPDATES) - {(r["source"], r["source_key"], r["clause"]) for r in rows}
    for source, key, line, clause, homes, status, targets, quote, evidence in ADDITIONS:
        assert all(not (r["source"] == source and r["source_key"] == key) for r in ledger), key
        assert label_lines(SOURCES[source]).get(key) == line, (key, label_lines(SOURCES[source]).get(key))
        rows.append({"source": source, "source_line": str(line), "source_kind": "supplement", "source_key": key,
                     "clause": clause, "coverage": "explicit", "draft_status": status, "draft_labels": targets,
                     "draft_quote": quote, "planned_homes": homes, "evidence_or_owed": evidence})
    for row in rows:
        roles = row["draft_status"].split("|")
        assert not set(roles) - STATUS and "absent" not in roles, row["draft_status"]
        targets = row["draft_labels"].split()
        assert targets and all(t in fragment_labels or t in draft_labels for t in targets), targets
        quote = " ".join(row["draft_quote"].split())
        if all(t.startswith(("part:", "sec:", "subsec:")) for t in targets):
            assert quote, row["clause"]
        if quote:
            assert fragment_text.count(quote) == 1, (fragment_text.count(quote), quote)
    with (HERE / f"{NAME}_rows.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    def counts(selection):
        return {c: sum(r["coverage"] == c for r in selection) for c in ("explicit", "recoverable", "owed")}

    def after(selection, source=None):
        now = counts(selection)
        moved = [r for r in rows[:len(UPDATES)] if source in (None, r["source"])]
        was = {c: sum(1 for r in ledger if (r["source"], r["source_key"], r["clause"]) in UPDATES
                      and r["coverage"] == c and source in (None, r["source"])) for c in ("recoverable", "owed")}
        added = sum(1 for a in ADDITIONS if source in (None, a[0]))
        return {"explicit": now["explicit"] + len(moved) + added, "recoverable": now["recoverable"] - was["recoverable"],
                "owed": now["owed"] - was["owed"]}

    print(f"{len(rows)} proposed rows: {len(UPDATES)} made explicit (before: {before}), {len(ADDITIONS)} for a passage without a row")
    print(f"ledger now: {len(ledger)} rows {counts(ledger)}; if accepted: {len(ledger) + len(ADDITIONS)} rows {after(ledger)}")
    for source in SOURCES:
        selection = [r for r in ledger if r["source"] == source]
        print(f"{source} now: {counts(selection)}; if accepted: {after(selection, source)}")
    print(f"fragment labels: {len(fragment_labels)}, none shared with the manuscript, LG or TAUM")


if __name__ == "__main__":
    main()
