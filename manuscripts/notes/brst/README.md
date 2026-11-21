# BRST working notes

Two working notes that explore a BRST-type reading of the polar correspondence. They are
separate from the consolidated manuscript and are not part of it.

| File | Contents |
|---|---|
| `brst_shift_note.tex` | *The shift example as a BRST model* (5 pp.). The charge of a preparation, its cohomology, gauge-invariant operators, a spectral obstruction, and the shift example of manuscript Section 7. |
| `brst_families_note.tex` | *Several presentations: redundancy and the reduced space* (4 pp.). The charge of a chain of maps and two treatments of redundancy among presentations. |
| `test_brst_shift.py`, `test_brst_families.py` | Standalone checks of every finite-dimensional claim (15 and 7 checks). |

Section and equation numbers refer to `manuscripts/consolidated/Universal_Model.tex` at
commit `67acdd8`.

## Status

- One assumption departs from the manuscript: the Krein form serves as the inner product of
  the reduced theory. The manuscript uses the Hilbert form for probabilities.
- Section 5 of the first note records the reading adopted by the author ("neither one nor
  two", "empty but not nonexistent") and marks it as interpretation. The mathematics does not
  depend on it.
- The first note was revised after an earlier review. Identifying two descriptions of an
  event is distinguished from removing it. The ground event is exact only under strong
  compatibility. A ghost grading is built from the polar reflection.

## Build and check

    pdflatex brst_shift_note.tex       # run three times; no bibliography
    pdflatex brst_families_note.tex
    python test_brst_shift.py
    python test_brst_families.py

The checks need NumPy and SciPy. They run on their own and are not collected by
`checks/run_all.py`.

## Points for review, most important first

1. Proposition 1 of the first note and Theorem 1 of the second: nilpotence, Krein
   self-adjointness, and the cohomology with signs alternating by degree.
2. Proposition 5 and Corollary 6 of the first note: the claim that gauge consistency implies
   the spectral exclusion that manuscript Section 7.3 assumes.
3. Proposition 4 of the second note: collective redundancy for three presentations and its
   inclusion–exclusion formula.
4. Section 5 of the first note: the revisions listed above.
5. Whether every interpretive statement is marked as such.
