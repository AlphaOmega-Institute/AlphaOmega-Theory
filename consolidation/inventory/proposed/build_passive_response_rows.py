#!/usr/bin/env python3
"""Write and check the proposed ledger rows for the passive-response fragment.

    python consolidation/inventory/proposed/build_passive_response_rows.py

The rows have the columns of migration.csv. They are proposals: `draft_labels` and
`draft_quote` refer to manuscripts/consolidated/fragments/passive_response.tex, which is
not part of Universal_Model.tex, so they do not belong in migration.csv until the text
has a permanent home. The script starts from the current rows for TAUM lines 2864-3610,
applies the updates below, appends the clauses the port separated out, and applies the
validator's rules with the fragment as the draft passage.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
FRAGMENT = ROOT / "manuscripts/consolidated/fragments/passive_response.tex"
DRAFT = ROOT / "manuscripts/consolidated/Universal_Model.tex"
STATUS = {"absent", "definition", "hypothesis", "model choice", "identification",
          "proved result", "conjecture", "open question", "worked example"}
FIELDS = ["source", "source_line", "source_kind", "source_key", "clause", "coverage",
          "draft_status", "draft_labels", "draft_quote", "planned_homes", "evidence_or_owed"]

# (source_key, clause) -> (draft_status, draft_labels, draft_quote, evidence)
UPDATES = {
    ("ass:selfadjoint-response", "Hilbert-self-adjoint generator, bounded port, nonreal resolvent"):
        ("hypothesis", "sec:passive-response eq:passive-response prop:preparation-split",
         "Throughout this section $\\mathcal X$ is a complex Hilbert space",
         "The opening fixes a self-adjoint H in a Hilbert space, a bounded port from a port space of any dimension, and the nonreal resolvent. Section 1 states the same hypotheses on the output space."),
    ("def:boundary-response", "Compressed resolvent of a bounded port"):
        ("definition", "eq:passive-response", "",
         "The displayed equation defines the response on the complement of the real axis for a port space of any dimension."),
    ("thm:response-Herglotz", "Holomorphy off the real axis and adjoint symmetry"):
        ("proved result", "prop:passive-herglotz", "",
         "Stated and proved from norm analyticity and the adjoint symmetry of the resolvent."),
    ("thm:response-Herglotz", "Positive imaginary part and its exact factorization on the upper half-plane"):
        ("proved result", "prop:passive-herglotz eq:passive-imaginary-part", "",
         "The proof compresses the one-point resolvent identity and exhibits the imaginary part as a positive operator."),
    ("def:Nevanlinna-kernel", "Two-point Nevanlinna kernel"):
        ("definition", "eq:passive-nevanlinna-kernel", "", "Defined on pairs of points of the upper half-plane."),
    ("thm:positive-response-kernel", "Two-point resolvent factorization and positive kernel"):
        ("proved result", "prop:passive-kernel eq:passive-kernel-factor eq:passive-kernel-positive", "",
         "The two-point resolvent identity gives the factorization. The quadratic sum is displayed as a squared norm for every finite choice of points and port vectors."),
    ("cor:canonical-transition-response", "Positive self-adjoint return on the input sector"):
        ("proved result", "cor:passive-sources eq:index-cohomology", "",
         "Part (a) uses the positive self-adjoint return stated after the cohomology formula of Section 7."),
    ("cor:canonical-transition-response", "Canonical compressed resolvent for the A*A choice"):
        ("definition", "eq:passive-canonical-response", "", "The canonical response is displayed for every bounded input into the physical space."),
    ("cor:canonical-transition-response", "Herglotz and positive-kernel conclusions for this choice"):
        ("proved result", "cor:passive-sources", "", "Part (a) applies both propositions with H equal to the return."),
    ("cor:Hilbert-block-response", "Port identity on the retained sector"):
        ("proved result", "cor:passive-sources eq:passive-block-response eq:prepared-schur", "",
         "Part (b) compresses the Schur formula of Section 8 by a port that maps into the retained space."),
    ("cor:Hilbert-block-response", "Passive response of the full Hilbert block"):
        ("proved result", "cor:passive-sources eq:passive-block-response", "",
         "The block operator is self-adjoint as a bounded symmetric perturbation of the direct sum of its diagonal blocks, so its response through the embedded port satisfies both propositions."),
    ("rem:positivity-boundary", "Limit of positivity for the Krein branch"):
        ("proved result", "eq:passive-signed-witness subsec:passive-positivity", "No such positivity is assumed here.",
         "Established by a witness: a two-by-two generator self-adjoint for the indefinite form, with real eigenvalues 1 and 2, whose response through the first coordinate has imaginary part -6/5 at z = 2 + i/2. Which further assumptions restore positivity is a separate question and is not part of this clause."),
    ("def:Schur-response", "Operator Cayley transform of the passive response"):
        ("definition", "eq:passive-schur-response", "", "Defined on the upper half-plane."),
    ("lem:Cayley-denominator", "Bounded invertibility of the Cayley denominator"):
        ("proved result", "prop:passive-cayley", "",
         "First statement of the proposition, with the lower bounds for the operator and its adjoint, closed dense range, and an inverse of norm at most one for port spaces of any dimension."),
    ("thm:Schur-contractivity", "Analyticity, defect identity and contractivity"):
        ("proved result", "prop:passive-cayley eq:passive-contractive", "", "The congruence with four times the imaginary part is derived and the norm bound follows."),
    ("def:dBR-kernel", "Operator-valued de Branges–Rovnyak kernel"):
        ("definition", "eq:passive-dbr-kernel", "", "Defined with the denominator and sign convention of TAUM."),
    ("thm:positive-Schur-kernel", "Kernel factorization, positivity and reproducing space"):
        ("proved result", "prop:passive-schur-kernel eq:passive-schur-kernel-factor", "",
         "The factorization with the two Cayley factors and the factor two is proved, and positivity follows by congruence. The reproducing-kernel consequence is stated with a citation and is not reproved."),
    ("def:compressed-spectral-measure", "Positive operator weights at the response lines"):
        ("definition", "eq:passive-spectral-measure cor:finite-response-defect", "",
         "The paragraph after the definition identifies the atoms of the measure with the weights of Corollary 1.3."),
    ("def:compressed-spectral-measure", "Borel measure and total mass for arbitrary spectra"):
        ("definition", "eq:passive-spectral-measure", "",
         "Defined on all Borel sets with total mass, weak countable additivity, and the scalar complex measures used for integration."),
    ("thm:spectral-temporal-response", "Finite-line resolvent and Fourier sums"):
        ("proved result", "eq:passive-finite-lines prop:passive-spectral-temporal", "",
         "Both finite sums are displayed after the proposition, for weights of any rank."),
    ("thm:spectral-temporal-response", "General resolvent and time-domain representations"):
        ("proved result", "prop:passive-spectral-temporal eq:passive-spectral-temporal", "",
         "Both operator-valued integrals are stated for an arbitrary compressed spectral measure, with the factor of Planck's constant used in the draft."),
    ("thm:spectral-temporal-response", "Time-domain positive definiteness"):
        ("proved result", "eq:passive-time-positive", "", "The quadratic sum is displayed as the squared norm of a sum of evolved port vectors."),
    ("def:response-generated-subspace", "Closed span generated from the port in time"):
        ("definition", "prop:passive-two-spans prop:preparation-split", "", "Section 1 defines the space by the time span. The proposition restates it for a port."),
    ("def:response-generated-subspace", "Nonreal resolvent-span definition, equivalence, and minimality terminology"):
        ("definition|proved result", "prop:passive-two-spans", "is called \\emph{minimal} when",
         "The proposition proves that the resolvent span and the time span coincide for arbitrary spectra, by the Laplace representation of the resolvent in one direction and invariance under all resolvents in the other. The sentence after it defines a minimal pair."),
    ("thm:canonical-minimal-reduction", "Contains the direct range and is invariant under the entire unitary group"):
        ("proved result", "prop:passive-two-spans prop:preparation-split", "", "Stated in the proposition and proved from the strong limit of the scaled resolvent."),
    ("thm:canonical-minimal-reduction", "Reducing self-adjoint restriction and exact response invariance"):
        ("proved result", "prop:passive-minimal-reduction eq:passive-minimal-invariance", "",
         "The restriction is self-adjoint, the two responses agree at every nonreal point, and the restricted pair is minimal."),
    ("thm:minimal-uniqueness", "Unique unitary intertwining generators and ports for equal responses"):
        ("proved result", "prop:passive-uniqueness eq:passive-unitary-equivalence eq:passive-gram", "",
         "Full proof: equal Gram matrices of resolvent-port vectors, extension by minimality, intertwining of resolvents and hence of generators with their domains, recovery of the ports, and uniqueness. Two continuity steps are written out: the conjugate pair of points in the Gram formula, which TAUM mentions, and equal points in the resolvent identity, which it does not."),
    ("cor:finite-response-defect", "Rank count for the indirect space and decomposition of the full defect"):
        ("proved result", "cor:finite-response-defect prop:preparation-split", "", "Unchanged. Section 1 states and proves the count and the decomposition of the adjoint kernel. The fragment repeats the hypothesis of finite-rank weights where it uses the count."),
    ("cor:finite-response-defect", "Port-coordinate invariance and the rank condition for injectivity (following discussion)"):
        ("proved result", "subsec:passive-minimal", "an invertible change of coordinates $T$ on $\\mathcal E$ preserves it",
         "The first remark gives invariance under a change of port coordinates. The last sentences give the kernel condition for an injective port."),
    ("cor:finite-response-defect", "Relations between the weight ranges (following discussion)"):
        ("proved result", "subsec:passive-minimal", "The count is the dimension of the space of relations",
         "The second remark proves that the range of the summed weight is the sum of the ranges, identifies the count with the dimension of the space of relations, and gives the direct-sum criterion for a zero count."),
    ("prop:permanent-port-invisibility", "Invisible complement and observability of indirect directions"):
        ("proved result", "prop:preparation-split", "", "Unchanged. Proposition 1.1 identifies the permanently invisible vectors with the orthogonal complement of the reached space."),
    ("prop:opposite-sector-response", "Orbital realization with and without the filter"):
        ("worked example", "eq:rotated-orbital-readout subsec:orbital-apparatus subsec:passive-sectors",
         "The trap of Section~\\ref{subsec:orbital-apparatus} is the instance",
         "Section 8 works the trap out. The fragment names it as the instance of the general proposition, with and without the filter."),
    ("prop:opposite-sector-response", "General matched-measure factorization and departure"):
        ("proved result", "prop:passive-opposite-sectors eq:passive-sector-factorization eq:passive-sector-departure", "",
         "Stated for strongly commuting H0 and L and equal scalar measures in opposite charge sectors: the orthonormal triple, the isometric port, the cosine factor, the weight on the indirect space and the time of complete departure."),
    ("prop:opposite-sector-response", "Fixed analysis line and weight outside the three-dimensional carrier"):
        ("proved result", "eq:passive-sector-factorization eq:passive-sector-departure subsec:passive-sectors",
         "If the two shifted spectral branches overlap",
         "The amplitude on the fixed line and the weight outside the three vectors are in the proposition. The caveat that the line need not lie in the reached space is kept, with an explicit overlapping example that TAUM does not give."),
    ("ass:distributional-trace", "Bounded trace between a Hamiltonian Hilbert scale and a port space"):
        ("definition|hypothesis", "subsec:passive-heat", "A \\emph{distributional trace} is a bounded operator",
         "The scale, its anti-dual with the carrier as pivot, the bounded trace and its adjoint are defined before the proposition."),
    ("thm:heat-smoothed-distributional-port", "Bounded smoothed port, Herglotz response and finite spectral mass"):
        ("proved result", "prop:passive-heat-distributional eq:passive-heat-port", "",
         "The port is defined as the adjoint of a bounded operator on the carrier, which gives TAUM's formula without separately extending the Gaussian to the negative scale. The statement names the eight propositions that apply to a bounded port. The paragraph after it says that the finite count and the sector proposition keep their own hypotheses."),
    ("ass:discrete-response-spectrum", "Complete eigenbasis and formal port coefficients for heat smoothing"):
        ("hypothesis", "subsec:passive-heat", "Suppose $H$ has a complete orthonormal basis of eigenvectors",
         "The eigenbasis, real eigenvalues and coefficient vectors in the port space are specified with no assumption of finite spectral support."),
    ("thm:heat-smoothed-channel", "Weighted summability, Hilbert–Schmidt port, response series and mass"):
        ("proved result", "prop:passive-heat-discrete eq:passive-heat-discrete-port eq:passive-heat-discrete-response", "",
         "The port, its Hilbert-Schmidt norm, the response series with locally uniform norm convergence and the total mass with trace-norm convergence are stated and proved."),
    ("rem:heat-channel-normalization", "Support normalization of a bounded port by its polar part"):
        ("proved result", "subsec:passive-heat eq:polar-active-supports", "is the port normalized on its support",
         "The polar part restricted to the support is identified as the normalized port and tied to the polar port of Section 1."),
    ("rem:heat-channel-normalization", "Application to heat-smoothed scalar and operator-valued ports"):
        ("proved result", "subsec:passive-heat", "A smoothed port has no preferred scale.",
         "The scalar normalization to unit mass and the operator-valued case are given, including a mass spectrum that accumulates at zero."),
}

# Clauses the port separates out: (source_key, clause, draft_status, draft_labels, draft_quote, evidence)
ADDITIONS = [
    ("def:boundary-response", "Scope of the definition: interior, finite-rank and regularized singular inputs (following discussion)",
     "definition", "sec:passive-response", "It covers a port inside the system, a port of finite rank, and a singular input",
     "TAUM's sentence after the definition is carried in the opening."),
    ("cor:finite-response-defect", "Enlarging the port can change the count (following discussion)",
     "worked example", "subsec:passive-minimal eq:pair-two-setting-ranges", "Third, enlarging the port can change the count.",
     "TAUM asserts this without an example. The fragment gives one from the pair model: an interior preparation has count one, and adding the endpoint preparation makes it zero. The check also gives a case where the count rises."),
    ("prop:permanent-port-invisibility", "Scope of the invisibility statement: one port and all real times (following discussion)",
     "proved result", "subsec:passive-minimal prop:semibounded-waiting", "These statements concern one port and all real times.",
     "TAUM's caution about finite observation windows and other observables is restated, with a pointer to the semibounded waiting proposition and the translation example of Section 1."),
    ("prop:opposite-sector-response", "Matched measures from a unitary reflection or time reversal (following discussion)",
     "proved result", "subsec:passive-sectors", "follow from a unitary reflection, or from an antiunitary symmetry such as time reversal",
     "Carried with the reason the antiunitary case still gives equal measures."),
    ("prop:opposite-sector-response", "Ground needs a no-crossing condition; complete departure without complete return (following discussion)",
     "proved result", "subsec:passive-sectors", "Three limits of the proposition should be kept.",
     "Both cautions are carried: Lg = 0 keeps g an eigenvector only, and the return weight at the stated times is the squared zero-field amplitude."),
]

NEW_KEYS_LINES = {"def:boundary-response": ("2883", "definition"), "cor:finite-response-defect": ("3294", "corollary"),
                  "prop:permanent-port-invisibility": ("3364", "proposition"), "prop:opposite-sector-response": ("3409", "proposition")}


def labels(path):
    return set(re.findall(r"\\label\{([^}]+)\}", path.read_text(encoding="utf-8")))


def main():
    with (INVENTORY / "migration.csv").open(encoding="utf-8", newline="") as stream:
        current = [r for r in csv.DictReader(stream) if r["source"] == "TAUM" and 2864 <= int(r["source_line"]) <= 3610]
    fragment_labels, draft_labels = labels(FRAGMENT), labels(DRAFT)
    assert not fragment_labels & draft_labels, fragment_labels & draft_labels
    for source in (ROOT / "manuscripts/lg/LG_Master_Unified.tex", ROOT / "manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex"):
        assert not fragment_labels & labels(source), fragment_labels & labels(source)
    fragment_text = " ".join(FRAGMENT.read_text(encoding="utf-8").split())
    rows, used = [], set()
    for row in current:
        key = (row["source_key"], row["clause"])
        status, targets, quote, evidence = UPDATES[key]
        used.add(key)
        rows.append({**row, "coverage": "explicit", "draft_status": status, "draft_labels": targets,
                     "draft_quote": quote, "evidence_or_owed": evidence})
    assert used == set(UPDATES), set(UPDATES) - used
    for key, clause, status, targets, quote, evidence in ADDITIONS:
        line, kind = NEW_KEYS_LINES[key]
        rows.append({"source": "TAUM", "source_line": line, "source_kind": kind, "source_key": key, "clause": clause,
                     "coverage": "explicit", "draft_status": status, "draft_labels": targets, "draft_quote": quote,
                     "planned_homes": "III.4", "evidence_or_owed": evidence})
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
    with (HERE / "passive_response_rows.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    before = {c: sum(r["coverage"] == c for r in current) for c in ("explicit", "recoverable", "owed")}
    print(f"{len(rows)} proposed rows for {len({r['source_key'] for r in rows})} source entries "
          f"({len(current)} existing clauses updated, {len(ADDITIONS)} separated out)")
    print(f"existing clauses before: {before}; proposed: all explicit")
    print(f"fragment labels: {len(fragment_labels)}, none shared with the draft, LG or TAUM")


if __name__ == "__main__":
    main()
