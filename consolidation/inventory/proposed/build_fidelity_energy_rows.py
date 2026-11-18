#!/usr/bin/env python3
"""Write and check the proposed ledger rows for the fidelity-energy fragment.

    python consolidation/inventory/proposed/build_fidelity_energy_rows.py

The rows have the columns of ../migration.csv. They are proposals: `draft_labels` and
`draft_quote` refer to manuscripts/consolidated/fragments/fidelity_energy.tex, which is
not part of Universal_Model.tex, so ../check_migration.py would reject them. They wait
here until the text is integrated and has its final labels.

The script takes the current row for each clause named below, checks that it is not yet
explicit, replaces its coverage, status, labels, quote and evidence, appends one row for
a source passage that had none, writes fidelity_energy_rows.csv, and applies the
validator's rules with the fragment as the draft passage. No explicit row is touched.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
NAME = "fidelity_energy"
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
    ("TAUM", "thm:reference-fidelity", "Arbitrary unit test state in Dom H"):
        ("proved result", "prop:fidelity-energy eq:fidelity-energy-factor", "",
         "Part (a) holds for every unit vector of the form domain, which contains Dom H, with the case of fidelity one stated separately and without reference to a preparation."),
    ("TAUM", "def:conditional-excitation-scale", "Nonnegative conditional excitation scale for any eligible state"):
        ("definition|proved result", "eq:fidelity-energy-factor", "",
         "The scale Z is defined as the excitation energy of the normalized orthogonal component and shown nonnegative. The text says that it depends on the state, when it is constant, and bounds it by the gap and by the width of a bounded spectrum."),
    ("TAUM", "cor:fidelity-gap", "General inequality on the full quadratic-form domain"):
        ("proved result", "eq:fidelity-energy-form eq:fidelity-energy-gap", "",
         "The excitation energy is defined on the form domain by the square root, and part (b) proves the inequality for every unit vector there under the form gap. A state of the number operator with finite mean and infinite second moment shows that the form domain is strictly larger."),
    ("TAUM", "eq:FS-gap-bound", "Fubini–Study form of the arbitrary input-state energy-gap estimate"):
        ("proved result", "eq:fidelity-energy-gap eq:fidelity-energy-distance", "",
         "The distance is defined with the modulus of the overlap, measured from the ray of the ground, and the bound is displayed in that form."),
    ("LG", "eq:fidelity-state-decomposition", "Reconstruction of a unit vector from ground fidelity, relative phase and normalized orthogonal component"):
        ("proved result", "eq:fidelity-energy-decomposition", "",
         "Displayed for fidelity strictly between zero and one, with the two boundary cases: no phase at zero, no conditional vector at one."),
    ("LG", "prop:ground-fidelity-factorization", "Arbitrary normalized form-domain states on either side"):
        ("proved result", "prop:fidelity-energy eq:fidelity-energy-factor eq:fidelity-energy-transport", "",
         "Part (a) is proved on the square-root form domain from the fact that the square root annihilates the ground. The text applies it to the input and to the output separately and says what the polar map does and does not transport."),
    ("LG", "cor:ground-gap-fidelity", "Arbitrary normalized finite-form-energy states on either side"):
        ("proved result", "eq:fidelity-energy-gap", "",
         "Part (b), with the gap stated as an inequality of forms on the form domain."),
    ("LG", "cor:BRLC-FS-isometry", "Fubini–Study distance on the aligned graph ray space"):
        ("proved result", "eq:fidelity-energy-transport subsec:fidelity-energy",
         "The graph maps $G_\\varphi$ of Section~\\ref{sec:neutral-graphs-phase} are isometries as well",
         "The graph maps, the aligned one among them, are isometries, so the rays of two vectors, of their images and of their graph vectors have the same distances."),
    ("LG", "thm:projective-ground-recovery", "Recovery of fidelity and Fubini–Study convergence from the ratio of components"):
        ("proved result", "eq:fidelity-energy-ratio", "",
         "The fidelity is displayed as a function of the ratio of the two components, convergence to the ground ray is equivalent to that ratio tending to zero, and the statement is specialized to determinations acting on the ground."),
    ("LG", "prop:no-state-below-ground", "Spectral and form lower bounds; unique ground ray when assumed"):
        ("proved result", "eq:fidelity-energy-form prop:fidelity-energy", "",
         "After the definition the text states that the excitation energy is nonnegative on the whole form domain, using only that the level is the bottom of the spectrum, and that uniqueness of the ground ray is a further assumption. Part (b) shows that a gap implies it."),
}

# (source, source_key, line, clause, planned_homes, draft_status, draft_labels, draft_quote, evidence)
ADDITIONS = [
    ("TAUM", "eq:cyclic-reference", 1343,
     "The fidelity results apply to any self-adjoint realization whose ground and domain hypotheses hold, and use no cyclicity",
     "III.1", "proved result", "subsec:fidelity-energy", "No cyclicity is used.",
     "TAUM lines 1327 to 1349. The last paragraph names the three operators to which the proposition can be applied, says that the hypotheses must be checked in each case, that a reference vector with nonzero image need not be a ground of the return and a ground in the kernel has zero image, and that cyclicity is not used."),
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
