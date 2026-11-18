#!/usr/bin/env python3
"""Write and check the proposed ledger rows for the relaxation fragment.

    python consolidation/inventory/proposed/build_relaxation_rows.py

The rows have the columns of ../migration.csv. They are proposals: `draft_labels` and
`draft_quote` refer to manuscripts/consolidated/fragments/relaxation.tex, which is not
part of Universal_Model.tex, so ../check_migration.py would reject them. They wait here
until the text is integrated and has its final labels.

The script takes the current row for each clause named below, checks that it is not yet
explicit, replaces its coverage, status, labels, quote and evidence, writes
relaxation_rows.csv, and applies the validator's rules with the fragment as the draft
passage. No row is added and no explicit row is touched.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
NAME = "relaxation"
FRAGMENT = ROOT / f"manuscripts/consolidated/fragments/{NAME}.tex"
DRAFT = ROOT / "manuscripts/consolidated/Universal_Model.tex"
TAUM = ROOT / "manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex"
LG = ROOT / "manuscripts/lg/LG_Master_Unified.tex"
STATUS = {"absent", "definition", "hypothesis", "model choice", "identification",
          "proved result", "conjecture", "open question", "worked example"}
FIELDS = ["source", "source_line", "source_kind", "source_key", "clause", "coverage",
          "draft_status", "draft_labels", "draft_quote", "planned_homes", "evidence_or_owed"]
THEOREM = "thm:canonical-deconditioning-semigroup"

# (source_key, clause) -> (draft_status, draft_labels, draft_quote, evidence). All LG.
UPDATES = {
    ("thm:exact-rigpa-density", "Equivalence of unit ground fidelity, zero deviation, and the pure ground state for normal states"):
        ("proved result", "prop:relaxation-pure-distance eq:relaxation-pure-distance", "",
         "The proposition ends with zero weight outside the ground exactly for the ground projection, and the sentence after its proof lists the three equivalent conditions. Here the equivalence follows from the upper bound on the trace distance."),
    ("def:ground-nonaltering-channel", "Ground-fixing CPTP channel"):
        ("definition", "subsec:relaxation-unitary", "The family \\emph{fixes the ground} if",
         "Stated for any family of channels on the appearance space. A single channel is the case of a family with one member."),
    ("thm:no-unitary-deconditioning", "Ground fidelity and deviation under every unitary preserving the ground ray"):
        ("proved result", "eq:relaxation-unitary", "",
         "The weight outside the ground is invariant under every unitary that fixes the ground projection, for every state. The text draws the consequence that such an evolution attracts only the ground itself, and that a unitary can move a chosen vector to the ground only by moving the ground."),
    ("def:ground-attractive-family", "Ground fixed and deviation tending to zero throughout a specified basin"):
        ("definition", "subsec:relaxation-unitary", "It \\emph{attracts} a set $\\mathsf B$ of states to the ground if",
         "Defined for a specified set of states, called the basin, with the remark that the condition is trace-norm convergence. Canonical deconditioning is the instance whose basin is every state."),
    (THEOREM, "Strong continuity on the trace class"):
        ("proved result", "prop:relaxation-semigroup", "",
         "Part (a) bounds the trace norm of the change of any trace-class operator by 4(1-exp(-γt)) times its norm. This is continuity in the operator norm on the trace class, which implies the strong continuity of the source."),
    (THEOREM, "(i) Ground projection fixed"):
        ("proved result", "prop:relaxation-semigroup", "",
         "Part (c) states that the ground projection is fixed."),
    (THEOREM, "(ii) Exact decay of ground deviation"):
        ("proved result", "eq:relaxation-decay", "",
         "Part (d), for every state, by tracing the block formula against the complement of the ground."),
    (THEOREM, "(iii) Unique fixed normal state at every positive time"):
        ("proved result", "prop:relaxation-semigroup", "",
         "Part (c): at a positive time the fixed trace-class operators are the multiples of the ground projection, so it is the only fixed state."),
    (THEOREM, "(iv) Trace-norm estimate with constant 2"):
        ("proved result", "eq:relaxation-distance eq:relaxation-pure-distance eq:relaxation-pure-input", "",
         "Part (e), from the bound for a pure target proved in the preceding proposition. The text adds a lower bound and shows from pure inputs that neither the constant 2 nor the rate γ/2 can be improved."),
    (THEOREM, "Trace-class generator at zero"):
        ("proved result", "eq:relaxation-generator eq:relaxation-jumps", "",
         "Part (b) gives the generator as a bounded operator on the trace class, with the relaxation as its exponential. The form with one jump operator for each direction orthogonal to the ground is displayed."),
    ("eq:ground-deviation-lyapunov", "Lyapunov law and finite-time ground recovery"):
        ("proved result", "eq:relaxation-decay subsec:relaxation-rate", "",
         "The paragraph after the jump form gives the differential law, strict decrease while the weight is positive, and that the ground is reached at a finite time only from the ground."),
    ("rem:reduced-irreversibility", "Environmental realization of irreversible reduced evolution"):
        ("proved result|open question", "eq:relaxation-dilation eq:relaxation-dilation-trace", "",
         "An isometry into the appearance space tensored with an environment gives each relaxation step as a reduced unitary evolution. The text says that this does not identify an environment in an observer, and that one unitary group for all times needs an infinite-dimensional larger space, which is not constructed."),
    ("def:momentary-ground-covariance", "Momentary ground covariance and contraction factors"):
        ("definition", "eq:relaxation-ground-carried eq:relaxation-contraction", "",
         "Both conditions are displayed, with the remark that the first concerns the channels and the second the particular sequence of states."),
    ("prop:exact-realization-absorbing", "Exact realization persists under momentary ground covariance"):
        ("proved result", "prop:relaxation-moments", "",
         "Part (a)."),
    ("thm:momentary-ground-stabilization", "Variable contraction product and asymptotic stabilization"):
        ("proved result", "prop:relaxation-moments eq:relaxation-product", "",
         "Part (b): the product bound, the trace-norm bound in each space, and convergence when the products tend to zero. The text notes that the bound does not use the first condition, and part (c) shows that a step contracting every state carries the ground."),
    ("def:realization-event-stabilization", "Exact event and asymptotic stabilization"):
        ("definition", "subsec:relaxation-moments", "We say that the sequence \\emph{is at the ground} at moment $n$ when",
         "Both notions are defined and distinguished by the deconditioning steps, under which a sequence off the ground is never at the ground and stabilizes when the rates times durations have a divergent sum. The cognitive names are given in one sentence."),
    ("def:conditional-projective-actualization", "Projective actualization and its exact domain"):
        ("definition", "eq:relaxation-return subsec:relaxation-return", "",
         "The returned state is defined for density operators. For a pure state it is the ray of V*φ, defined exactly when φ is not in the orthogonal complement of D."),
    ("prop:return-operation-trace", "Trace of the return operation, trace nonincrease and normalized successful branch"):
        ("proved result", "eq:relaxation-return", "",
         "The compression has trace equal to the weight in D, the map is completely positive and does not increase the trace, and the normalized state is displayed with its pairing. The text says when the trace is a probability."),
    ("thm:deconditioning-contracts-excess", "Complementary probability of physical presentation"):
        ("proved result", "eq:relaxation-paired-weight", "",
         "Displayed for every appearance state when the target lies in D."),
    ("eq:restricted-return-after-deconditioning", "Restricted conditional return gives the vacuum at positive time"):
        ("proved result", "eq:relaxation-restricted-return", "",
         "The restricted return of the relaxed state is the physical reference projection times a coefficient that is positive at every positive time."),
    ("eq:global-return-after-deconditioning", "Global conditional return after deconditioning"):
        ("proved result", "eq:relaxation-return-bound", "",
         "For every initial appearance state the normalized return converges in trace norm to the partner of the ground, with an explicit bound. The partner is the vacuum exactly under strong polar compatibility, and the two returns are contrasted on the unequal-cell path."),
}


def labels(path):
    return set(re.findall(r"\\label\{([^}]+)\}", path.read_text(encoding="utf-8")))


def main():
    with (INVENTORY / "migration.csv").open(encoding="utf-8", newline="") as stream:
        ledger = list(csv.DictReader(stream))
    fragment_labels, draft_labels = labels(FRAGMENT), labels(DRAFT)
    assert not fragment_labels & draft_labels, fragment_labels & draft_labels
    for source in (LG, TAUM):
        assert not fragment_labels & labels(source), fragment_labels & labels(source)
    fragment_text = " ".join(FRAGMENT.read_text(encoding="utf-8").split())
    rows, before = [], {"explicit": 0, "recoverable": 0, "owed": 0}
    for row in ledger:
        key = (row["source_key"], row["clause"])
        if row["source"] != "LG" or key not in UPDATES:
            continue
        assert row["coverage"] != "explicit", key
        before[row["coverage"]] += 1
        status, targets, quote, evidence = UPDATES[key]
        rows.append({**row, "coverage": "explicit", "draft_status": status, "draft_labels": targets,
                     "draft_quote": quote, "evidence_or_owed": evidence})
    assert len(rows) == len(UPDATES), set(UPDATES) - {(r["source_key"], r["clause"]) for r in rows}
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
    # every clause of the theorem is explicit once these rows are accepted
    theorem = [r for r in ledger if r["source_key"] == THEOREM]
    updated = {r["clause"] for r in rows if r["source_key"] == THEOREM}
    assert all(r["coverage"] == "explicit" or r["clause"] in updated for r in theorem)
    with (HERE / f"{NAME}_rows.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    total = {c: sum(r["coverage"] == c for r in ledger) for c in ("explicit", "recoverable", "owed")}
    after = {"explicit": total["explicit"] + len(rows), "recoverable": total["recoverable"] - before["recoverable"],
             "owed": total["owed"] - before["owed"]}
    lg = {c: sum(r["coverage"] == c for r in ledger if r["source"] == "LG") for c in ("explicit", "recoverable", "owed")}
    lg_after = {"explicit": lg["explicit"] + len(rows), "recoverable": lg["recoverable"] - before["recoverable"],
                "owed": lg["owed"] - before["owed"]}
    print(f"{len(rows)} proposed rows, all made explicit; before: {before}")
    print(f"the theorem: {len(theorem)} clauses, {len(theorem) - len(updated)} already explicit, {len(updated)} made explicit")
    print(f"ledger now: {len(ledger)} rows {total}; if accepted: {len(ledger)} rows {after}")
    print(f"LG now: {lg}; if accepted: {lg_after}")
    print(f"fragment labels: {len(fragment_labels)}, none shared with the manuscript, LG or TAUM")


if __name__ == "__main__":
    main()
