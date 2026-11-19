#!/usr/bin/env python3
"""Write and check the proposed ledger rows for the ground-preservation fragment.

    python consolidation/inventory/proposed/build_ground_preservation_rows.py

The rows have the columns of ../migration.csv. They are proposals: `draft_labels` and
`draft_quote` refer to manuscripts/consolidated/fragments/ground_preservation.tex, which is
not part of Universal_Model.tex, so ../check_migration.py would reject them. They wait
here until the text is integrated and has its final labels.

The script takes the current row for each clause named below, checks that it is not yet
explicit, replaces its coverage, status, labels, quote and evidence, writes
ground_preservation_rows.csv, and applies the validator's rules with the fragment as the
draft passage. No explicit row is touched and no row is added.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
NAME = "ground_preservation"
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
    ("LG", "prop:strong-polar-consequences", "Output ground has positive eigenvalue a² for AA*"):
        ("proved result", "eq:ground-output-return", "",
         "Displayed with its domain statement, derived from condition (v) and the bridge. The text adds that the weak bridge alone does not give it, on the unequal-cell path."),
    ("LG", "prop:vacuum-not-return-zero-mode", "Ground with any nonzero image in the product-domain statement"):
        ("proved result", "eq:return-kernels", "",
         "The quadratic form of the input return is the squared norm of the image, so the kernels of the two returns are the kernels of A and of its adjoint. No bridge is used."),
    ("LG", "eq:shifted-physical-return-model", "Constraints on taking structural returns as proposed ground generators"):
        ("proved result|worked example", "eq:return-kernels eq:ground-output-return",
         r"It constrains any attempt to take the returns as the sector Hamiltonians.",
         "The unshifted input return cannot have the ground at energy zero under the bridge, and the ground is an eigenvector of the return only under strong compatibility. LG's shifted return on the active input has the ground at its bottom exactly when the bridge coefficient squared is the bottom of the return there; an example with a kernel shows that on the whole input space the kernel lies below. R3 and R5 alone exclude the output return as a Hamiltonian selecting the output ground, which is orthogonal to its whole ground eigenspace. On the unequal-cell path and in the pair preparation the returns select other states."),
    ("LG", "prop:no-iterated-awareness-tower", "Forward image, vanishing square, and annihilation of every output vector"):
        ("proved result", "eq:index-no-second-step", "",
         "Displayed for the differential on its whole domain, with the same statement for its negative and the remark that the model generates no iterated sequence."),
    ("LG", "id:rangrig-terminal-law", "No second act is needed to make the occurrence manifest"):
        ("identification", "eq:index-no-second-step part:problem-origin",
         r"On that reading no second act is needed to make it manifest.",
         "The consequence of the opening's reading, that a cognition is manifest in its occurrence, is stated next to the vanishing of the differential on the output summand. Whether that vanishing expresses reflexivity is left as a further identification."),
    ("LG", "def:ground-line-preserving-operator", "Preservation of the ground line and projective survival"):
        ("definition|proved result", "eq:ground-line-kept", "",
         "For a bounded operator the block condition, the vector condition and the vanishing of the variance are equated. Keeping the ground ray requires a nonzero mean, and an annihilating operator keeps the line without keeping the ray."),
    ("LG", "prop:ground-preservation-reduction", "One-sided preservation and two-sided reduction"):
        ("proved result", "eq:ground-line-kept",
         r"Both $X$ and $X^*$ map the line into itself exactly when $X$ commutes with $P_y$",
         "One-sided preservation is shown not to imply the two-sided condition, by an explicit operator. The two-sided condition is commutation with the ground projection and makes the line reducing. For a unitary one side gives the other, and for a self-adjoint operator the condition is an eigenvector."),
    ("LG", "thm:unitary-ground-preservation", "Unitary preservation, scalar character and generator eigenvector for all real times"):
        ("proved result", "eq:relaxation-unitary-group", "",
         "The three conditions are displayed and proved through the character of the group and the spectral theorem, including the domain of the generator. The return times form a closed subgroup of the real line, so a single return is not enough, as an example shows, and two returns with irrational ratio are."),
    ("LG", "prop:terminality-nonalteration-distinct", "Typed terminality places no restriction on independent ground dynamics"):
        ("proved result|worked example", "eq:index-no-second-step eq:index-channel-coherence",
         r"It places no condition on how a channel treats the grounds.",
         "The vanishing of the differential on the output summand holds for every map. In the shift example one preparation is shared by a channel that carries the input ground to the output ground and by channels that send it to the missing line."),
    ("LG", "rem:annihilation-not-destruction", "Meaning of annihilation by the typed forward operator"):
        ("proved result", "eq:index-no-second-step",
         r"It does not describe an evolution that removes $y$. It records that the map is not applied to output vectors.",
         "The zero on the output summand is attributed to the typing of the map and separated from any evolution of the ground, which the next sentences show can keep or move it."),
    ("LG", "def:momentary-phenomenal-system", "Rank-one requirement on each moment's ground projection"):
        ("definition|hypothesis", "subsec:relaxation-moments eq:ground-spectral-projection",
         r"This is the uniqueness part of R2 imposed at that moment.",
         "Each moment's ground projection is required to be the whole spectral projection of that moment's Hamiltonian at the bottom of its spectrum, which is the uniqueness part of R2. A projection onto a chosen vector has rank one in any case."),
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
