#!/usr/bin/env python3
"""Write and check the proposed ledger rows for the moving-subspace fragment.

    python consolidation/inventory/proposed/build_moving_subspace_rows.py

The rows have the columns of migration.csv. They are proposals: `draft_labels` and
`draft_quote` refer to manuscripts/consolidated/fragments/moving_subspace.tex, which is
not part of Universal_Model.tex, so they cannot enter migration.csv until the text has a
permanent home. The script starts from the current rows for TAUM lines 2008-2508, applies
the updates below, appends rows for source passages that had none, and applies the
validator's rules with the fragment as the draft passage.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
NAME = "moving_subspace"
FRAGMENT = ROOT / f"manuscripts/consolidated/fragments/{NAME}.tex"
DRAFT = ROOT / "manuscripts/consolidated/Universal_Model.tex"
TAUM = ROOT / "manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex"
LG = ROOT / "manuscripts/lg/LG_Master_Unified.tex"
FIRST, LAST = 2008, 2508
STATUS = {"absent", "definition", "hypothesis", "model choice", "identification",
          "proved result", "conjecture", "open question", "worked example"}
FIELDS = ["source", "source_line", "source_kind", "source_key", "clause", "coverage",
          "draft_status", "draft_labels", "draft_quote", "planned_homes", "evidence_or_owed"]

# (source_key, clause) -> (draft_status, draft_labels, draft_quote, evidence)
UPDATES = {
    ("lem:projection-tangent", "Self-adjoint derivative and its two zero diagonal blocks"):
        ("proved result", "prop:moving-tangent eq:moving-diagonal-zero eq:moving-block-forms", "",
         "Stated for a vector field of class C1 and proved by differentiating the projection identity. The decomposition into the two off-diagonal blocks is the first block form."),
    ("def:tangent-channel", "The channel of a specified tangent direction"):
        ("definition", "eq:moving-tangent-map", "",
         "The tangent map, its adjoint and the coordinate abbreviation are defined. The text adds that the map depends only on the value of the vector field at the point."),
    ("thm:tangent-realization", "Two block forms and their oriented decompositions"):
        ("proved result", "eq:moving-block-forms eq:moving-two-odd", "",
         "Both block forms are proved. The two sums are written with the differential of Section 7, as d + d* and d* - d. TAUM's letters for the oriented pair are not used, because the manuscript has not introduced them."),
    ("@2155", "Grassmannian tangent and its Krein branch"):
        ("proved result", "eq:moving-two-odd subsec:moving-derivative", "It is the tangent vector to the family of subspaces.",
         "The derivative is identified as the Hilbert-self-adjoint tangent, and its product with the moving symmetry as Hilbert-skew-adjoint and self-adjoint for the form of that point."),
    ("def:Grassmann-connection", "Projected connection on the range bundle"):
        ("definition", "eq:moving-connection", "",
         "Defined for sections of the moving subspaces, with its compatibility with the inner product."),
    ("thm:projector-curvature", "Commutator, compression and curvature as a mixed return"):
        ("proved result", "prop:moving-curvature eq:moving-curvature eq:moving-curvature-components", "",
         "Proved for vector fields of class C1, including the bracket term, and specialized to commuting coordinate fields as in TAUM."),
    ("def:mixed-return", "Ordered mixed product of two tangent channels"):
        ("definition", "eq:moving-mixed-return subsec:moving-connection", "Call it the \\emph{mixed return} of the two directions.",
         "Defined as the product of the adjoint of one tangent map with the other, on the selected subspace. In the notation of Section 7 it is d_Y* d_Z, which equals TAUM's signed product of the oriented pair."),
    ("cor:curvature-AV", "Antisymmetrization of oriented excursions"):
        ("proved result", "eq:moving-mixed-return", "",
         "The curvature is displayed as the difference of the two mixed returns, with their order, in the differential notation of Section 7."),
    ("@2238", "Directional mixed returns and nonzero curvature"):
        ("proved result", "subsec:moving-connection", "Curvature is thus not the return of one map.",
         "The contrast between one positive return and the antisymmetric exchange of two directions is stated, with the consequence that a one-parameter family has no curvature. The worked example shows two tangent maps into the same complex line whose imaginary ratio makes the mixed return anti-Hermitian."),
    ("def:local-gauge-data", "Frame connection, adjoint sign and curvature matrix"):
        ("definition", "eq:moving-frame eq:moving-connection-matrix eq:moving-frame-curvature prop:moving-frame", "",
         "Frames, the connection matrix and its curvature are defined. Part (a) of the proposition gives the anti-Hermitian sign. Local existence of frames is shown by a formula."),
    ("thm:frame-covariance", "Fiber compression and change of orthonormal frame"):
        ("proved result", "prop:moving-frame eq:moving-frame-covariance eq:moving-change-of-frame", "",
         "Parts (b) and (c). The proof uses the action of the connection on coefficient vectors in place of TAUM's differential-form identity."),
    ("def:holonomy", "Loop transport in a local frame"):
        ("definition", "eq:moving-frame-transport eq:moving-transport", "",
         "The transport equation in a frame and the unitary matrix of a loop are given. Transport is also defined without a frame, for subspaces of any dimension and loops that leave a frame domain."),
    ("@2326", "Distinction from geometric holonomy"):
        ("proved result", "subsec:moving-frames", "Holonomy in this sense should be kept apart from the sign",
         "The algebraic square of the polar rotation is distinguished from transport around a loop."),
    ("eq:fiberwise-Laplacian", "Fiberwise fixed-channel identities for tangent channels"):
        ("proved result", "eq:moving-squares eq:moving-two-odd", "",
         "Both squares are proved. The paragraph after the proposition states that each tangent map is an instance of the single-map structure of Sections 2, 3 and 7 and that the new content is the relation between directions."),
    ("thm:chern-character", "Closed form, invariant class and first Chern form"):
        ("proved result", "prop:moving-chern eq:moving-chern-character eq:moving-first-chern", "",
         "Closedness is proved directly from the oddness of the derivative. Independence of the connection and the identification with the Chern character are cited to the Chern-Weil theorem and not reproved. The first Chern form and its components are displayed."),
    ("eq:higher-chern-character", "Higher-degree Chern character"):
        ("definition|proved result", "eq:moving-chern-character", "",
         "Displayed for every degree with its factorial and normalization, and shown to equal the trace of the projection times a power of its derivative."),
    ("cor:line-flux", "Integral first Chern flux"):
        ("proved result", "cor:moving-flux eq:moving-flux", "",
         "Proved for every finite rank by a triangulation and the determinants of edge transports. TAUM's statement is the case of rank one. The sentence after the proof covers smooth singular cycles. The identification of the integer with the Chern pairing is cited."),
    ("@2477", "Holonomy and characteristic number"):
        ("proved result", "eq:moving-holonomy-flux subsec:moving-characteristic", "The holonomy of a loop is not quantized.",
         "The determinant of the holonomy of a loop bounding a disc is the exponential of minus the flux of the trace through it, which need not be trivial for a contractible loop. The integer belongs to a closed surface. The paragraph also separates curvature at a point from holonomy of a loop."),
    ("eq:ch2-return", "Second Chern character and interpretation of characteristic numbers"):
        ("proved result", "eq:moving-second-chern subsec:moving-characteristic", "Identifying a characteristic number with a physical charge requires additional data",
         "The degree-four form is displayed. The caution about physical charge is carried."),
}

# Source passages that had no row: (source_key, source_line, clause, draft_status, draft_labels, draft_quote, evidence)
ADDITIONS = [
    ("sec:projector-geometry", 2009, "Tangent operators of a moving subspace are odd and fall under the single-map results",
     "proved result", "prop:moving-tangent subsec:moving-derivative", "A moving subspace adds no new algebra at one point and in one direction.",
     "The opening claim of TAUM's section. Anticommutation with the moving symmetry is part of the proposition."),
    ("eq:projection-field", 2022, "Norm-C2 field of orthogonal projections, locally constant dimensions and the bundle of ranges",
     "definition|proved result", "eq:moving-projection-field prop:moving-unitary", "",
     "TAUM asserts local constancy of the dimensions and the bundle structure. The fragment proves them from the direct rotation of two nearby projections."),
    ("eq:moving-J", 2047, "Symmetry and splitting induced at each point, with a fixed splitting as the constant case",
     "definition", "eq:moving-symmetry", "",
     "Defined, with the remark that the Hilbert inner product is fixed while the indefinite form varies with the point."),
    ("eq:projector-summary", 2360, "Summary: the symmetry turns the tangent into the form-self-adjoint operator, and curvature is the antisymmetrized mixed return",
     "proved result", "eq:moving-two-odd eq:moving-mixed-return", "",
     "Both halves of TAUM's boxed summary are displayed equations."),
    ("eq:curvature-form", 2381, "Curvature two-form of a smooth finite-rank field and its components",
     "definition", "eq:moving-curvature-form", "",
     "Defined with the convention that relates the two-form to its coordinate components."),
    ("eq:physics-flux-integral", 2474, "Hermitian convention for the flux",
     "proved result", "subsec:moving-characteristic", "The convention here is anti-Hermitian.",
     "The corollary is restated for the Hermitian two-form."),
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
    current = [r for r in ledger if r["source"] == "TAUM" and FIRST <= int(r["source_line"]) <= LAST]
    fragment_labels, draft_labels = labels(FRAGMENT), labels(DRAFT)
    assert not fragment_labels & draft_labels, fragment_labels & draft_labels
    for source in (LG, TAUM):
        assert not fragment_labels & labels(source), fragment_labels & labels(source)
    fragment_text = " ".join(FRAGMENT.read_text(encoding="utf-8").split())
    source_lines = label_lines(TAUM)
    rows, used, unchanged = [], set(), 0
    for row in current:
        key = (row["source_key"], row["clause"])
        if key not in UPDATES:
            assert row["coverage"] == "explicit", key          # already stated in the manuscript: kept as it is
            rows.append(dict(row))
            unchanged += 1
            continue
        status, targets, quote, evidence = UPDATES[key]
        used.add(key)
        rows.append({**row, "coverage": "explicit", "draft_status": status, "draft_labels": targets,
                     "draft_quote": quote, "evidence_or_owed": evidence, "_new": True})
    assert used == set(UPDATES), set(UPDATES) - used
    existing_keys = {r["source_key"] for r in ledger if r["source"] == "TAUM"}
    for key, line, clause, status, targets, quote, evidence in ADDITIONS:
        assert key not in existing_keys and source_lines.get(key) == line, (key, source_lines.get(key))
        assert FIRST <= line <= LAST
        rows.append({"source": "TAUM", "source_line": str(line), "source_kind": "supplement", "source_key": key,
                     "clause": clause, "coverage": "explicit", "draft_status": status, "draft_labels": targets,
                     "draft_quote": quote, "planned_homes": "III.2", "evidence_or_owed": evidence, "_new": True})
    rows.sort(key=lambda r: int(r["source_line"]))
    seen = set()
    for row in rows:
        clause = (row["source_key"], row["clause"])
        assert clause not in seen, clause
        seen.add(clause)
        if not row.pop("_new", False):
            continue
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
    print(f"{len(rows)} proposed rows: {len(UPDATES)} updated, {unchanged} unchanged, {len(ADDITIONS)} for passages without a row")
    print(f"clauses of this section before: {before}; proposed: all explicit")
    print(f"ledger now: {len(ledger)} rows {total}; if accepted: {len(ledger) + len(ADDITIONS)} rows {after}")
    print(f"fragment labels: {len(fragment_labels)}, none shared with the manuscript, LG or TAUM")


if __name__ == "__main__":
    main()
