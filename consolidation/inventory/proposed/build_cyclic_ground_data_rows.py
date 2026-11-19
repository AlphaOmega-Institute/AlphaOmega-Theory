#!/usr/bin/env python3
"""Write and check the proposed ledger rows for the cyclic-ground-data fragment.

    python consolidation/inventory/proposed/build_cyclic_ground_data_rows.py

The rows have the columns of ../migration.csv. They are proposals: `draft_labels` and
`draft_quote` refer to manuscripts/consolidated/fragments/cyclic_ground_data.tex, which is
not part of Universal_Model.tex, so ../check_migration.py would reject them. They wait
here until the text is integrated and has its final labels.

The script takes the current row for each clause named below, checks that it is not yet
explicit, replaces its coverage, status, labels, quote and evidence, writes
cyclic_ground_data_rows.csv, and applies the validator's rules with the fragment as the
draft passage. No explicit row is touched and no row is added.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
NAME = "cyclic_ground_data"
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
    ("LG", "def:pointed-representation", "General abstract pointed quadruple"):
        ("definition", "eq:pointed-vector-state sec:determination",
         r"Call $(\mathfrak M,\pi,\mathcal H,y)$ a pointed representation.",
         "An abstract unital C*-algebra with a nondegenerate representation and a unit vector is named, without a Hamiltonian or cyclicity. The closed image serves as the concrete algebra of the section."),
    ("LG", "def:vector-state", "Vector functional for any represented algebra"):
        ("definition", "eq:pointed-vector-state", "",
         "The vector state is defined on the abstract algebra and identified with the concrete functional composed with the representation."),
    ("LG", "prop:vector-state-positive", "Linearity, positivity and normalization of the vector functional"):
        ("proved result", "eq:pointed-vector-state",
         r"It is linear, positive and normalized, and it is the functional $\omega$ on $\mathfrak D$ composed with $\pi$.",
         "The three properties are displayed together, with the phase invariance of the state."),
    ("LG", "def:ground-datum", "Equivalent rank-one spectral projection criterion"):
        ("proved result", "eq:ground-spectral-projection",
         r"the ground ray is unique exactly when $P_0$ has rank one",
         "The spectral projection at the bottom of the spectrum projects onto the eigenspace. It is nonzero exactly when the bottom is an eigenvalue and has rank one exactly when the ground ray is unique."),
    ("LG", "rem:ground-properties-independent", "Failure of the reverse deductions"):
        ("worked example|proved result", "subsec:ground-properties eq:ground-gap-condition",
         r"Each of the three properties can fail while the other two hold.",
         "Three examples: the five-site path with the cell algebra is unique and gapped but not cyclic; two separate paths with the site algebra are cyclic and gapped but not unique; a Hamiltonian with eigenvalues 1/n and the algebra of the unilateral shift is unique and cyclic but not gapped, with a variant in continuous spectrum."),
    ("LG", "prop:cyclic-is-gns", "Unitary from the vector-state GNS representation to any cyclic pointed representation"):
        ("proved result", "eq:gns-intertwiner", "",
         "The map from GNS classes to vectors is displayed with its two intertwining properties. It is an isometry onto the generated space, a unitary in the cyclic case, and the only unitary with those properties."),
    ("LG", "def:pure-faithful-state", "Purity, irreducibility, abstract and concrete faithfulness"):
        ("definition|proved result", "subsec:ground-properties prop:pure-faithful-scalars",
         r"Its pullback to $\mathfrak M$ is faithful exactly when, in addition, $\pi$ is injective.",
         "Purity as extremality and its equivalence with irreducibility of the GNS representation are stated, and in the cyclic case with triviality of the commutant. Faithfulness is defined, identified with the separating property on the concrete algebra, and distinguished for the abstract algebra by an example."),
    ("LG", "def:physical-cyclic-subspace", "Closed orbit of any physical vector"):
        ("definition", "sec:determination",
         r"It depends only on the ray of $\psi$ when $\psi\ne0$, and it is $\{0\}$ when $\psi=0$.",
         "The closed subspace generated from any vector is defined for the specified algebra, and the text applies the construction to the physical vacuum with the algebra of R1."),
    ("LG", "prop:physical-centered-cyclic", "Centered physical orbit spans the vacuum-orthogonal subspace"):
        ("proved result", "prop:cyclic-determination-content sec:determination",
         r"the centered vectors $X^\circ x$ have closed span $\mathcal H_\Psi\ominus\C x$",
         "The proposition is applied on the input side with the physical vacuum and the algebra of R1."),
    ("LG", "def:dark-complement", "Orthogonal complement of the closure of the luminous orbit"):
        ("definition", "sec:determination",
         r"LG calls $K^\perp$ the dark complement.",
         "The complement of the generated space is shown to be invariant, so that its projection commutes with the algebra, and is named as LG does."),
    ("LG", "prop:no-dark-sector", "Vanishing of the orthogonal complement under R2"):
        ("proved result", "sec:determination subsec:ground-properties",
         r"It is $\{0\}$ exactly when $y$ is cyclic, so R2 excludes it for the algebra that R2 specifies.",
         "Stated with the remark that a vector orthogonal to the ground can lie in the generated space. On the unequal-cell path the complement for the cell algebra is the defect, while for an algebra with a cyclic ground the defect lies in the generated space."),
    ("LG", "q:luminous-purity", "Open status of purity and faithfulness of luminous vector state"):
        ("open question|proved result", "prop:pure-faithful-scalars subsec:ground-properties",
         r"It does not say whether the vector state of the luminous ground on the algebra of R2 is pure or faithful, and LG leaves both open.",
         "The two questions are posed separately. The proposition shows that a state is pure and faithful only on the scalars, so at most one holds for a ground in a space of dimension at least two. Examples on the five-site path show either one, and the text distinguishes technical purity from the expression pure luminosity."),
    ("LG", "choice:determination-algebra", "Cyclicity of the selected algebra"):
        ("hypothesis", "sec:determination",
         r"LG chooses a subalgebra $\mathfrak D$ of the algebra in R2 to represent the determinations and assumes, unless stated otherwise, that $y$ is cyclic for $\mathfrak D$.",
         "Stated as LG's default, with the remark that it is R2 for the full specified algebra and a further hypothesis for a smaller one."),
    ("LG", "prop:covariance-content-inner-product", "Variance equals squared norm of the projected component"):
        ("proved result", "eq:determination-variance", "",
         "Displayed, with the statement that the variance vanishes exactly when the image of the ground lies on its ray."),
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
