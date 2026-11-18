#!/usr/bin/env python3
"""Write and check the proposed ledger rows for the four-maps fragment.

    python consolidation/inventory/proposed/build_four_maps_rows.py

The rows have the columns of migration.csv. They are proposals: `draft_labels` and
`draft_quote` refer to manuscripts/consolidated/fragments/four_maps.tex, which is not
part of Universal_Model.tex, so they cannot enter migration.csv until the text has a
permanent home. The script starts from the current rows for TAUM lines 3628-4048, applies
the updates below, keeps the row that is already explicit, appends rows for two source
passages that had none, and applies the validator's rules with the fragment as the draft
passage. It also writes one owed row agreed in review for a sentence of LG that the
ledger did not track.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
INVENTORY = HERE.parent
ROOT = INVENTORY.parents[1]
NAME = "four_maps"
FRAGMENT = ROOT / f"manuscripts/consolidated/fragments/{NAME}.tex"
DRAFT = ROOT / "manuscripts/consolidated/Universal_Model.tex"
TAUM = ROOT / "manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex"
LG = ROOT / "manuscripts/lg/LG_Master_Unified.tex"
FIRST, LAST = 3628, 4048
STATUS = {"absent", "definition", "hypothesis", "model choice", "identification",
          "proved result", "conjecture", "open question", "worked example"}
FIELDS = ["source", "source_line", "source_kind", "source_key", "clause", "coverage",
          "draft_status", "draft_labels", "draft_quote", "planned_homes", "evidence_or_owed"]

# (source_key, clause) -> (draft_status, draft_labels, draft_quote, evidence)
UPDATES = {
    ("ex:scalar-channel", "One-dimensional polar phase, support, index and diagonal Laplacian"):
        ("worked example", "eq:four-maps-scalar subsec:four-maps-scalar", "the polar part is the phase $u=b/|b|$",
         "The Laplacian is displayed. For a nonzero number the supports, the polar phase and the index zero are stated, and the case of the number zero is added."),
    ("ex:scalar-channel", "Minimal scalar resolvent"):
        ("worked example", "eq:four-maps-scalar-response", "",
         "The return on one line with the identity port gives one line of weight one. The pair is minimal and is the smallest case of the canonical response."),
    ("ex:scalar-channel", "Scalar example as the two signed matrices"):
        ("worked example", "eq:four-maps-scalar", "",
         "Both odd matrices are displayed for a complex number, with the sign of the lower entry."),
    ("ex:split-excitation-response", "Bias inequality, full resolvent and time dependence for the named branches"):
        ("worked example", "eq:four-maps-three-state-data eq:four-maps-three-state-response eq:four-maps-three-state-transfer", "",
         "The example is restated in TAUM's parameters, with the inequality that keeps the ground unique, the full response matrix, and the cosine and sine of the transfer with Planck's constant explicit."),
    ("ex:translation-channel", "Translation as an unbounded channel with continuous partner spectrum"):
        ("worked example", "eq:four-maps-translation eq:four-maps-translation-returns", "",
         "The derivative on its Sobolev domain, its self-adjointness, the two equal returns on their domain, the Fourier multipliers and the absolutely continuous spectrum are stated."),
    ("eq:no-translation-gap", "Absence of support gap and nonclosed range for the derivative"):
        ("proved result", "eq:four-maps-translation-no-gap", "",
         "Low-frequency unit vectors show that no lower bound holds on the support. The closed-range criterion of Section 7 gives a range that is not closed. The text adds that it is dense and that the map is not Fredholm."),
    ("eq:free-response", "Continuous-spectrum response of a specified translation-channel port"):
        ("proved result|worked example", "eq:four-maps-translation-response", "",
         "The Fourier integral is displayed with its domain. A second form gives the absolutely continuous measure on the spectrum of the return, and the text shows that a one-dimensional port is never minimal here."),
    ("thm:oscillator-channel", "Creation-channel domains, partner operators and index"):
        ("proved result", "prop:four-maps-oscillator eq:four-maps-ladder eq:four-maps-oscillator-returns eq:four-maps-oscillator-index", "",
         "The ladder operators and their domain, the two returns, the lower bound that gives a closed range, the kernels and the index are stated and proved."),
    ("thm:oscillator-channel", "Positive integer spectrum"):
        ("proved result", "prop:four-maps-oscillator", "",
         "The spectra of the two returns are stated in the proposition, with the extra zero of the second. The proposition also identifies the polar isometry with the shift of Section 7."),
    ("cor:oscillator-heat-index", "Linear and squared oscillator heat traces"):
        ("proved result", "eq:four-maps-oscillator-heat", "",
         "Both differences are displayed and proved from one identity for summable functions of the number operator. The first is identified with the heat-trace formula of Section 7."),
    ("thm:oscillator-zeta", "Mellin integral of the positive integer heat trace"):
        ("proved result", "prop:four-maps-mellin eq:four-maps-mellin", "",
         "Stated for real part of s greater than one and proved by integrating term by term, which the convergent sum of the integrals of the absolute values permits."),
    ("@3901", "Scope of the completed zeta factor"):
        ("proved result", "subsec:four-maps-oscillator", "It makes no statement about the zeros of the zeta function.",
         "The factor is distinguished from the entire xi function and its pole at one is noted. The calculation uses the positive integers as eigenvalues of multiplicity one, and it says nothing about the zeros."),
    ("eq:oscillator-response", "Cyclic geometric oscillator port and its minimal response"):
        ("proved result|worked example", "eq:four-maps-oscillator-port eq:four-maps-oscillator-response", "",
         "The unit vector with geometric coefficients, its response with poles at the positive integers, minimality from the simple spectrum, and uniqueness by the proposition of Section 9."),
    ("ex:quasiperiodic-channel", "Boundary domain, spectrum, partners, index and gap"):
        ("worked example", "eq:four-maps-interval-operator eq:four-maps-interval-spectrum", "",
         "The domain with the twisted condition, the eigenvalues and normalized eigenfunctions, the equal returns, the index zero, the discrete spectrum and the gap are stated."),
    ("eq:unsmoothed-point-mass", "Uniform point weights and divergent unsmoothed mass"):
        ("proved result", "eq:four-maps-interval-point", "",
         "Evaluation at a point is identified as a distributional trace on the scale of Section 9 for exponent greater than one half. Its coefficients have equal squared modulus and their sum diverges."),
    ("eq:theta-mass", "Heat-smoothed point port and its finite mass"):
        ("proved result", "eq:four-maps-interval-port eq:four-maps-interval-mass", "",
         "The smoothed vector and its squared norm are displayed, independent of the point. A second form by Poisson summation is added."),
    ("eq:quasiperiodic-smoothed-response", "Signed-energy first-order spectrum and its scalar response"):
        ("proved result", "eq:four-maps-interval-response", "",
         "The response is displayed, with the remark that its energies are the eigenvalues of the first-order operator, of both signs."),
    ("eq:normalized-smoothed-point", "Normalization of the smoothed point port"):
        ("proved result", "subsec:four-maps-interval", "gives a spectral measure of total mass one",
         "Normalizing the vector gives unit mass and leaves the positions of the lines unchanged."),
    ("@4015", "Size of the quasiperiodic gap"):
        ("proved result", "subsec:four-maps-interval", "the smallest value of $|\\lambda_n|$ is $\\pi/\\Lambda$",
         "Both cases of the twist are given, with the square of the gap for the returns."),
}

# An owed clause agreed in review of the previous block: (source, source_key, clause, planned_homes, owed text)
OWED = [
    ("LG", "eq:vacuum-symmetry-phase", "The phase freedom of a unique ground ray is not a model of individual differences among sentient beings",
     "V", "State that the phase a symmetry gives a unique ground is the freedom of its ray and explains no difference between sentient beings. LG says so in the sentence after the equation, and Section 11 carries only the mathematics."),
]

# Source passages that had no row: (source_key, source_line, clause, draft_status, draft_labels, draft_quote, evidence)
ADDITIONS = [
    ("sec:examples", 3629, "The common relations do not fix the spectrum, and the realizations differ in spectral type",
     "worked example", "sec:four-maps subsec:four-maps-table", "They do not follow from those relations.",
     "The opening of TAUM's section and its closing comparison. The fragment opens with the same statement and ends with the table of the four maps."),
    ("sec:three-state-response", 3671, "The port read as a preparation: unit return, defect line, ground correspondence, and a response Hamiltonian that is neither return",
     "worked example", "subsec:four-maps-three-state", "It is neither $A^*A$ nor $AA^*$.",
     "The paragraph after TAUM's example. The fragment also carries its last remark, that the prepared plane is not invariant although the Hamiltonian commutes with the conserved generator."),
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
                     "draft_quote": quote, "planned_homes": "IV", "evidence_or_owed": evidence, "_new": True})
    owed_rows = []
    for source, key, clause, homes, text in OWED:
        anchor = [r for r in ledger if r["source"] == source and r["source_key"] == key]
        assert anchor and all(r["clause"] != clause for r in anchor), key
        owed_rows.append({"source": source, "source_line": anchor[0]["source_line"], "source_kind": anchor[0]["source_kind"],
                          "source_key": key, "clause": clause, "coverage": "owed", "draft_status": "absent", "draft_labels": "",
                          "draft_quote": "", "planned_homes": homes, "evidence_or_owed": text})
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
        writer.writerows(rows + owed_rows)
    before = {c: sum(r["coverage"] == c for r in current) for c in ("explicit", "recoverable", "owed")}
    total = {c: sum(r["coverage"] == c for r in ledger) for c in ("explicit", "recoverable", "owed")}
    after = {"explicit": total["explicit"] + before["recoverable"] + before["owed"] + len(ADDITIONS),
             "recoverable": total["recoverable"] - before["recoverable"], "owed": total["owed"] - before["owed"] + len(OWED)}
    print(f"{len(rows) + len(owed_rows)} proposed rows: {len(UPDATES)} updated, {unchanged} unchanged, "
          f"{len(ADDITIONS)} for passages without a row, {len(OWED)} newly recorded as owed")
    print(f"clauses of this section before: {before}; proposed: all explicit")
    print(f"ledger now: {len(ledger)} rows {total}; if accepted: {len(ledger) + len(ADDITIONS) + len(OWED)} rows {after}")
    print(f"fragment labels: {len(fragment_labels)}, none shared with the manuscript, LG or TAUM")


if __name__ == "__main__":
    main()
