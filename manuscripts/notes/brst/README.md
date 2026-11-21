# BRST working notes

Two working notes that explore a BRST-type reading of the polar correspondence. They are
separate from the consolidated manuscript and are not part of it.

| File | Contents |
|---|---|
| `brst_shift_note.tex` | *The shift example as a BRST model*. The charge of a preparation, its cohomology, evolution on the quotient, gauge-invariant operators, a spectral obstruction, and the shift example of manuscript Section 7. |
| `brst_families_note.tex` | *Several presentations: redundancy and the reduced space*. The charge of a chain of maps and two treatments of redundancy among presentations. |
| `test_brst_shift.py`, `test_brst_families.py` | Standalone checks of every finite-dimensional claim (23 and 12 checks). |

Section and equation numbers refer to `manuscripts/consolidated/Universal_Model.tex` at
commit `67acdd8`.

## Status

- One assumption departs from the manuscript: the Krein form serves as the inner product of
  the reduced theory. The manuscript uses the Hilbert form for probabilities.
- Two physical premises are used, and the choice of the Krein form imposes neither.
  - For the results on dynamics and admissibility, the observer's evolution acts on the
    quotient defined by the preparation. For a sector-preserving unitary evolution this is
    equivalent to commuting with the charge, and hence to intertwining by the polar map. For a
    family it can be imposed on each member's quotient or on the joint quotient. For general
    dynamics the two conditions are independent, and the spectral exclusion holds under either.
    In the shift, for presentations satisfying R3, both are equivalent to intertwining by
    every presentation, which is the manuscript's admissibility condition (Proposition 9).
  - To exclude transfers by interventions, we choose a sufficient premise: every admitted
    operation also acts on the quotient, so its Kraus operators map the range and the kernel
    of the charge into themselves. The first premise does not imply this one. It protects
    exact events, and for sector-preserving interventions the whole matched part.
- Reading. The reduced space records what remains unpaired by the specified presentation.
  *Empty* is represented by exactness with respect to a named charge.
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

## Revisions after the third review

Lemma 4, Proposition 6, Corollary 7 and the three-presentation result were accepted. Their
numbering is unchanged.

1. Interventions. The evolution premise constrains free evolution only. In the shift the
   dynamics preserves the charge, and the covariant operator taking the exact ground event to
   the missing line remains available. Excluding it uses a second, sufficient premise: every admitted
   operation acts on the quotient. This premise is stated in Section 4 of the first note and
   recorded in its abstract, Section 1 and Section 5.
2. Members and families. New Remark 8 of the first note: conservation of each member's charge
   and conservation of the joint charge are independent, with an example each way. Corollary 7
   now states the exclusion under either condition. Only the joint condition makes the
   evolution act on the family's reduced space.
3. Efficacy. The sentence calling the excess "less efficacious" is replaced. The excess has
   zero direct return through this presentation (the quarter-turn returns the same
   direction with probability one). Its efficacy in another operation is left open.
4. Exactness names a charge. On the unequal-cell path the relational ground is exact for the
   charge of the restricted ground map and has no class for the full charge.
5. Kinds of zero. Section 5 of the second note distinguishes a bare redundant input, which is
   not closed and so has no class, from its aligned pair, which is exact and has zero class.
   Commitment (iii) of the first note now says that the resolved reduced space equals the
   common mathematical defect. It also says that the conjecture asks for more, and that R4 on
   one preparation establishes neither the defect of a family nor that its directions are
   genuine.
6. Checks. The intersection helper stacks the complementary projectors, which separates two
   lines 1e-5 radians apart. New checks cover nearby lines, the bare redundancy and its
   aligned pair, both independence examples, the Kraus-representation independence of the
   second premise, exactness for the restricted ground charge, and the direct returns.

## Revisions after the fourth review

Proposition 6 and Corollary 7 were kept, and the revisions on emptiness and efficacy were
accepted.

1. Shift equivalence. New Proposition 9 of the first note: in the shift, for presentations
   with A f0 in C g1, conservation of every member's charge, conservation of the joint charge
   and intertwining by every presentation are equivalent, at resonance and away from it. The
   premise therefore recovers exactly the class of presentations the manuscript assumes. The
   ground condition cannot be dropped, as an example shows. The abstract and Section 5 state
   the result.
2. Remark 8. The independence is stated for general dynamics. The reverse example now has a
   shared lower ground, so that R3 holds for both members. Both examples use three input
   energies, whereas the shift has two.
3. Interventions, now a paragraph of their own.
   - The second premise is sufficient for excluding transfer from an exact event, not
     necessary: the sector sign J does not act on the quotient but preserves the
     missing-line weight.
   - For sector-preserving Kraus operators the premise is equivalent to commuting with the
     charge, and it protects the whole matched part.
   - For sector-mixing operators it protects exact events only. In an example, an operation
     acts on the quotient, leaves the exact ground event fixed, and moves half the weight of
     the bare input into the missing line.
   - The channels of Eq. (95), extended by zero, fail the premise under R4, so they enter as
     presentations.
   - Using Eq. (95) as a test requires saying how the preparation and intervention realize
     the graph event being tested.
4. Wording. The README and test comments call the rotation through pi/2 a quarter-turn.
5. Checks. New checks cover Proposition 9 on random families at three excitation energies,
   the role of the ground condition, the reverse example with R3, and each claim about the
   scope of the intervention premise.

## Revisions after the fifth review

Proposition 9, the revised Remark 8 and the claim about sector-preserving interventions
were accepted.

1. Chosen premise. The premise on interventions is called a chosen sufficient premise
   throughout. The clause tying sector preservation to the superselection rule is removed:
   the manuscript's Section 3 restricts what measurements distinguish and allows a sector
   swap among controls. The sentence on even operators in Section 3 of the first note now
   speaks of effects.
2. Carrier implementations. The zero extension of a channel from inputs to appearances is one
   implementation among several. The reflection R + P_N gives the same output on every bare
   input and fixes every exact event. Whether a channel satisfies the premise depends on its
   implementation.
3. A test that can exclude implementations. For a trace-preserving operation that preserves
   exact events, a Stinespring argument bounds the defect weight of a unit bare input by one
   half. The bound is zero when the whole output lies in the appearance sector. The
   sector-mixing example attains one half. Eq. (95) bounds the line weight from above. An
   independent observation of a positive line weight, with the whole output in the
   appearance sector, would exclude every implementation that preserves exact events.
4. Presentations. Only a finite Kraus family, such as the one attaining Eq. (95), is said to
   be treatable as presentations. The reset channel's joint row is not closable.
5. Checks. A new check covers the two implementations, the bound on random operations, its
   zero case and its attainment, and the reset channel.

## Points for review, most important first

1. The defect bound of Section 4 (p at most one half, zero for an output in the appearance
   sector) and its use with Eq. (95).
2. The two carrier implementations, and the revised statement on presentations.
