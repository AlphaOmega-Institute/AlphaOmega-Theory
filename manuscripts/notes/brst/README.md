# BRST working notes

Two working notes that explore a BRST-type reading of the polar correspondence. They are
separate from the consolidated manuscript and are not part of it.

| File | Contents |
|---|---|
| `brst_shift_note.tex` | *The shift example as a BRST model*. The charge of a preparation, its cohomology, evolution on the quotient, gauge-invariant operators, a spectral obstruction, and the shift example of manuscript Section 7. |
| `brst_families_note.tex` | *Several presentations: redundancy and the reduced space*. The charge of a chain of maps and two treatments of redundancy among presentations. |
| `test_brst_shift.py`, `test_brst_families.py` | Standalone checks of every finite-dimensional claim (19 and 10 checks). |

Section and equation numbers refer to `manuscripts/consolidated/Universal_Model.tex` at
commit `67acdd8`.

## Status

- One assumption departs from the manuscript: the Krein form serves as the inner product of
  the reduced theory. The manuscript uses the Hilbert form for probabilities.
- One physical premise is needed for the results on dynamics and admissibility: the
  observer's evolution acts on the quotient defined by the preparation. For a
  sector-preserving unitary evolution this is equivalent to commuting with the charge, and
  hence to intertwining by the polar map. The choice of the Krein form does not impose it.
- Reading. The reduced space records what remains unpaired by the specified presentation.
  The stronger reading, that it inventories real appearances, is noted and not adopted.
  Section 5 of the first note gives the author's terms for the relational event ("neither one
  nor two", "empty but not nonexistent") and marks them as interpretation. The mathematics
  does not depend on them.

## Build and check

    pdflatex brst_shift_note.tex       # run three times; no bibliography
    pdflatex brst_families_note.tex
    python test_brst_shift.py
    python test_brst_families.py

The checks need NumPy and SciPy. They run on their own and are not collected by
`checks/run_all.py`.

## Revisions after the second review

1. Spectral obstruction: shorter proof through the adjoint of the polar map, covering
   degenerate eigenvalues; the premise that the evolution acts on the quotient is stated, with
   a lemma; the condition is stated for unitary groups; only the polar map must intertwine.
2. Qualifiers: at resonance the single joint charge keeps the redundancies, and the reduced
   space is zero only after they are quotiented; "no excess gives a zero reduced space" needs
   R4; a second intertwiner leaves exactly the missing line only off resonance.
3. Three presentations: collective redundancy is a cancellation not assembled from pairwise
   cancellations, with an example in which every pair shares the ground; finite dimensions
   are stated; the concluding alternative is limited to three presentations.
4. Ground and grading: strong compatibility is stated where the anti-aligned ground is the
   partner of the aligned one; a definite reduced space contains no nonzero null class; every
   surviving class has ghost number zero; the use of Eq. (95) needs a specified preparation
   and readout.
5. Reading: as stated above, used consistently, with the relativity of a class illustrated by
   the shift rotated through pi/2.
6. Checks: the intersection helper uses an absolute tolerance, with tests for identical and
   full subspaces; new checks cover each revised claim.

## Points for review, most important first

1. Lemma 4 and Proposition 6 of the first note: evolution on the quotient, and the spectral
   obstruction through the adjoint of the polar map.
2. Corollary 7 of the first note and Section 3 of the second: the resonant cases.
3. Proposition 4 of the second note and its shared-ground example.
4. Section 5 of the first note: the adopted reading and its two consequences.
