#!/usr/bin/env python3
"""Write and check the proposed ledger rows for the polar-graph-clauses fragment.

    python consolidation/inventory/proposed/build_polar_graph_clauses_rows.py

The rows have the columns of ../migration.csv. They are proposals: `draft_labels` and
`draft_quote` refer to manuscripts/consolidated/fragments/polar_graph_clauses.tex, which is
not part of Universal_Model.tex, so ../check_migration.py would reject them. They wait
here until the text is integrated and has its final labels.

The script takes the current row for each clause named below, checks that it is not yet
explicit, replaces its coverage, status, labels, quote and evidence, writes
polar_graph_clauses_rows.csv, and applies the validator's rules with the fragment as the
draft passage. No explicit row is touched and no row is added.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
NAME = "polar_graph_clauses"
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
    # Section 2: transport of observables by the polar map
    ("LG", "prop:pure-projector-correspondence", "Conjugation of rank-one projections and equality of their trace overlaps"):
        ("proved result", "eq:matched-pure-projections", "",
         "Displayed: the transport of a rank-one projection is the projection of the transported vector, and the two trace overlaps equal the squared modulus of the inner product."),
    ("LG", "thm:alphaA-star-isomorphism", "Spectra and bounded functional calculus"):
        ("proved result", "prop:matched-born-probabilities eq:matched-spectral-measure",
         r"Hence $X$ and $VXV^*$ have the same spectrum, the second taken as an operator on $D$.",
         "The resolvent is transported, so the spectra agree on D. Polynomials and bounded Borel functions of a normal operator are transported, and positive operators, effects and projections keep their kind."),
    ("LG", "cor:paired-spectral-measures", "Spectral measure of the transported observable"):
        ("proved result", "eq:matched-spectral-measure", "",
         "Displayed, with the argument from uniqueness of the spectral measure. The text adds the unbounded self-adjoint case and the spectral measure of the extension by zero to the full output space."),
    # Section 3: graph maps in blocks, paired representation, graph states
    ("LG", "thm:WA-isometry", "Closed graph range and full adjoint formula"):
        ("proved result", "eq:graph-adjoint-blocks", "",
         "The adjoint is displayed for the active carrier and described for the full output space. The range is closed because the graph map is an isometry on a complete space, also when the range of A is not closed."),
    ("LG", "prop:phase-graph-projection", "Special cases for the aligned and anti-aligned projections"):
        ("proved result", "eq:graph-aligned-opposite", "",
         "Both projections are displayed with their product, their sum and the resulting orthogonal decomposition of the active carrier. The replacement of the identity on the full direct sum is stated."),
    ("LG", "def:paired-observable-embedding", "Unit of the selected pair algebra"):
        ("proved result", "sec:neutral-graphs-phase",
         r"$\iota(I_S)=I_S\oplus I_D$ is the identity of $\mathcal K_A$",
         "The image of the input identity is the identity of the active carrier and the active projection after extension by zero."),
    ("LG", "thm:paired-algebra-representation", "Normal faithful unital maximal representation, faithful selected restriction, and isomorphisms onto both sectoral algebras"):
        ("proved result", "sec:neutral-graphs-phase prop:matched-born-probabilities",
         r"The paired algebra $\iota(\mathcal M)$ is therefore isomorphic to $\mathcal M$ through the upper block and to $V\mathcal MV^*$ through the lower one.",
         "The map is linear, multiplicative, adjoint-preserving, unital, injective through its upper block and normal through both blocks. The restriction to a selected algebra and the two isomorphisms are stated."),
    ("LG", "thm:paired-algebra-graph-intertwining", "Graph compression identity"):
        ("proved result", "eq:graph-paired-compression", "",
         "Displayed, with the invariance of every graph. The text adds that a block-diagonal operator maps a graph into itself exactly when it is paired."),
    ("LG", "prop:mixed-state-correspondence", "Four blocks of the graph density operator"):
        ("proved result", "eq:graph-density-blocks", "",
         "All four blocks are displayed, with the opposite phase factors in the two off-diagonal blocks."),
    ("LG", "eq:mixed-sector-compressions", "Diagonal blocks of a mixed graph state"):
        ("proved result", "eq:graph-density-blocks",
         r"Its compressions to the two sectors are $\tfrac12\rho$ and $\tfrac12V\rho V^*$.",
         "Both compressions, their traces and their normalizations are stated. Deleting the off-diagonal blocks gives the equal mixture of the graph states at two opposite phases."),
    ("LG", "rem:direct-sum-not-partial-trace", "Meaning of the sector compressions of a graph density"):
        ("proved result", "eq:graph-density-blocks sec:neutral-graphs-phase",
         r"the lower block $\tfrac12V\rho V^*$ is a compression and is not the reduced state of a subsystem",
         "The partial trace over the sector coordinate returns the input density with trace one and is what the paired observables read. A sector compression keeps one diagonal block with trace one half. The carrier is a direct sum of sectors and not a tensor product of two systems."),
    ("LG", "prop:phase-graph-factorization", "Factorized ranges and projectors"):
        ("proved result", "sec:neutral-graphs-phase eq:graph-density-blocks",
         r"the graph $\mathcal N_\varphi$ is $\ket{\varphi}\otimes S$, and $\Pi_\varphi$ is $\ket{\varphi}\!\bra{\varphi}\otimes I_S$",
         "In the sector coordinates the graph, its projection and the graph state of a density operator are written as products with the sector vector at the given phase."),
    # Section 4: the restricted map and the plane of the two reference vectors
    ("LG", "prop:ground-bridge-equivalences", "Compression, vanishing orthogonal part, rank-one products, normalized projector and converse"):
        ("proved result", "eq:ground-plane-restriction subsec:ground-plane",
         r"the identity $(I-P_y)AP_x=0$ places $Ax$ on the line of $y$",
         "The compression identity, the vanishing orthogonal part, the rank-one product and the formula for the output projection are stated, with the converse for a nonzero image."),
    ("LG", "lem:ground-channel-bounded", "Norm and two positive products of the restriction"):
        ("proved result", "eq:ground-plane-restriction", "",
         "Displayed. The text adds that the bounded adjoint of the restriction extends the product of the input projection with the adjoint of A, with an example in which the output reference vector lies outside the domain of that adjoint."),
    ("LG", "thm:ground-channel-polar", "Unique factor, modulus, forward and adjoint vector actions"):
        ("proved result", "eq:ground-plane-restriction subsec:ground-plane",
         r"It acts by $u_0x=y$ and $u_0^*y=x$.",
         "The modulus is displayed. The line map has the kernel of the restriction, the polar decomposition with that property is unique, and both vector actions are stated."),
    ("LG", "prop:u0-uniqueness", "The bridge and its positive coefficient fix the line orientation"):
        ("proved result", "subsec:ground-plane",
         r"among them $AP_x=aw$ with $a>0$ holds for $w=u_0$ only",
         "The partial isometries between the two lines are the phase multiples of the line map, and the bridge with a positive coefficient holds for one of them. A common phase of the two representatives leaves it unchanged."),
    ("LG", "def:ground-graph-lines", "Anti-aligned companion of the restricted reference map"):
        ("definition", "eq:ground-plane-vectors", "",
         "The graph lines of the line map are defined at every phase under the weak bridge alone, and the aligned and anti-aligned lines are named."),
    ("LG", "thm:ground-graphs-neutral", "Total and maximal neutrality of both lines inside the ground plane"):
        ("proved result", "eq:ground-plane-vectors subsec:ground-plane",
         r"Each line is its own Krein-orthogonal complement in $\mathcal K_0$ and so is a maximal neutral subspace of the plane.",
         "The Krein product of a graph vector with a general vector of the plane is computed. It vanishes exactly on the line. The plane itself is not neutral."),
    ("LG", "prop:normalized-ground-pair", "Both vectors have Hilbert norm one, zero Krein square, and half weight in each sector"):
        ("proved result", "eq:ground-plane-vectors",
         r"Every $z^0_\varphi$ has Hilbert norm one, zero Krein square, and squared norm $1/2$ in each sector.",
         "Stated for every phase, with the normalized sector components."),
    ("LG", "prop:nonduality-without-identity", "Trivial sector intersections and bijective coordinate projections of the aligned ground graph"):
        ("proved result", "eq:ground-plane-vectors subsec:ground-plane",
         r"A line $\mathcal N^0_\varphi$ meets each sector in zero only, and each coordinate projection maps it bijectively onto the corresponding reference line.",
         "Stated for every phase. A vector of the line is fixed by either component, and neither component is a vector of the line."),
    ("LG", "def:ground-polar-operators", "Ground reflection, complex structure, and phase operator supported on the reference plane"):
        ("definition", "eq:ground-plane-operators", "",
         "The reflection and the complex structure of the line map are displayed as operators that vanish off the plane. The phase operator is written as minus i times the complex structure, as in Section 3."),
    ("LG", "thm:ground-polar-algebra", "Adjoints, grading, squares, and Pauli products of the three reference-plane operators"):
        ("proved result", "eq:ground-plane-operators subsec:ground-plane",
         r"They satisfy $\mathcal R_0^*=\mathcal R_0$, $\mathcal R_0^2=P_{xy}$, $J\mathcal R_0=-\mathcal R_0J$, $\mathcal I_0^*=-\mathcal I_0$ and $\mathcal I_0^2=-P_{xy}$.",
         "Adjoints, squares and the anticommutation with the grading are stated. The three operators are the Pauli matrices in the basis of the two reference vectors, with the three products written out."),
    ("LG", "prop:ground-graph-eigenspaces", "Reflection eigenlines; action of J and of the complex structure on both signs"):
        ("proved result", "eq:ground-plane-eigenlines", "",
         "Displayed for both vectors, with the projection onto the line at a general phase and the three means."),
    ("LG", "eq:Xi-Krein-cross-pairing", "Ground-line Witt pair"):
        ("proved result", "eq:ground-plane-witt", "",
         "The Hilbert product, the two Krein squares and the Krein cross product are displayed, and the pair is named."),
    ("LG", "prop:projectivized-neutral-locus", "Entire projective neutral cone of the reference plane as a circle, including its two reflection eigenlines"):
        ("proved result", "eq:ground-plane-vectors subsec:ground-plane",
         r"The neutral rays of the plane form a circle with coordinate $\varphi$.",
         "Proved from the two conditions on the coefficients. Equal sector weights define the circle and the bridge selects one point on it."),
    ("LG", "id:companion-branch-status", "Aligned and anti-aligned rays have identical sector compressions"):
        ("proved result", "eq:ground-plane-vectors subsec:ground-plane",
         r"The sector components of $z_-$ lie on the rays of $x$ and of $y$, as do those of $z_+$ and of every $z^0_\varphi$.",
         "Stated for every point of the circle. The points differ in relative phase only, have equal means for every operator commuting with the grading, and equal mean energy under sector Hamiltonians with the reference vectors as eigenvectors."),
    ("LG", "thm:ground-compression-full-polar", "Two-sided compressions of reflection and complex structure and inclusion of the anti-aligned line"):
        ("proved result", "prop:ground-plane-reduction eq:ground-plane-compression", "",
         "Under strong compatibility the plane reduces both polar operators to the plane operators and every graph line of the line map lies in the corresponding full graph. The proposition proves these as equivalences. Under the weak bridge alone the compressions are the plane operators times the square root of the ground fidelity."),
    ("LG", "prop:three-part-phenomenal-decomposition", "Ground, differentiated paired space, and adjoint kernel under the weak bridge"):
        ("proved result", "eq:ground-plane-three-parts", "",
         "Displayed with the corresponding split of the complement of the output reference vector. The middle part is paired by the polar map with the complement of the partner in the initial support."),
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
