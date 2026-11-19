# `cyclic_ground_data.tex`: pointed representations, cyclicity and properties of a ground

Sources: LG lines 1200 to 1880 ("Pointed cyclic representations and ground-state data",
"The physical cyclic vacuum", "The luminous phenomenal ground"), with
`choice:determination-algebra` (line 4563) and `prop:covariance-content-inner-product`
(line 4749). Base: `e0b1926`. The fragment has three parts, each headed by a marker that
names its place:

| Part | Place | Length |
|---|---|---|
| A | Section 5.1, after the paragraph that ends "whether K still fills H_Φ." | one paragraph |
| B | Section 5.1, after the paragraph that ends "requires additional information." | three paragraphs |
| C | Section 5, a new subsection after "Generating states and reading them on the five-site path" | about one and a half pages, one proposition |

    python manuscripts/consolidated/fragments/build_fragment.py cyclic_ground_data
    python checks/run_all.py -k cyclic_ground
    python consolidation/inventory/proposed/build_cyclic_ground_data_rows.py

Nothing here changes `Universal_Model.tex` or the shared ledger. The preview prints the
three parts one after another with their markers and numbers the subsection 13.1. In
place it is 5.3. `cyclic_ground_data_bibliography.tex` adds one entry, Kadison's 1957
paper on the transitivity theorem.

## Source results a reader can now find

| Source (LG) | Fragment |
|---|---|
| `def:pointed-representation`, the abstract quadruple | A, "Call (𝔐, π, ℋ, y) a pointed representation." |
| `def:vector-state`, `prop:vector-state-positive` | A, `eq:pointed-vector-state` and the sentence after it |
| `choice:determination-algebra`, cyclicity of the selected algebra | A, last three sentences |
| `def:physical-cyclic-subspace` | A, "For any vector ψ ..." |
| `prop:cyclic-is-gns` | B, `eq:gns-intertwiner` and the three sentences after it |
| `prop:covariance-content-inner-product`, the variance | B, `eq:determination-variance` |
| `def:dark-complement`, `prop:no-dark-sector` | B, third paragraph |
| `prop:physical-centered-cyclic` | B, last two sentences |
| `def:ground-datum`, the rank-one criterion | C, `eq:ground-spectral-projection` |
| `rem:ground-properties-independent`, failure of the reverse deductions | C, `eq:ground-gap-condition` and examples (a) to (c) |
| `def:pure-faithful-state` | C, the paragraph before the proposition |
| `q:luminous-purity` | C, `prop:pure-faithful-scalars`, the examples and the last paragraph |

## Change in the clause count

The block is the mathematical clauses of LG lines 1200 to 1880 that are not explicit,
with the two later clauses named above: 11 recoverable and 3 owed. I searched the whole
ledger for other clauses that are not explicit and name cyclicity, the GNS construction,
purity, faithfulness, vector states or spectral projections. The rows found belong to the
fundamental decomposition, the Schur reduction, the ground plane, or interpretive
headings, and none is carried here.

`consolidation/inventory/proposed/cyclic_ground_data_rows.csv` proposes the 14 clauses as
explicit. No row is added.

On the base the ledger goes from 496 rows (324 explicit, 78 recoverable, 94 owed) to 496
(338, 67, 91). LG goes from 181, 61, 63 to 195, 50, 60. Together with the fidelity and
ground-plane rows, which are accepted and waiting for integration, the ledger would have
497 rows (377, 33, 87), with LG at 229, 18, 58.

Clauses near these that the fragment does not carry:

- the second clause of `eq:vacuum-symmetry-phase` (line 1610), that the phase of a unique
  ground is not a model of differences among sentient beings;
- the second clause of `eq:no-nonluminosity-state` (line 1787), the doctrinal reading of
  complete absence of luminosity. Its mathematical half, that the zero vector and the
  zero functional are not states, is one sentence and could be carried with it;
- the identifications `id:matter-vacuum` and `id:fundamental-sectors`.

## What the proofs need

From the consolidated manuscript: R1 and R2; Proposition 5.1 with
`eq:determination-covariance`; the five-site path, its ground and gap, and the cell and site
algebras of Sections 4 and 5; the unequal-cell map and its defect; the gap condition after
`eq:ground-return-energy`.

Used without proof, with citations: the image of a C*-algebra under a representation is
norm closed; a state is pure exactly when its GNS representation is irreducible, and a
representation is irreducible exactly when its commutant is the scalars (Bratteli and
Robinson); Kadison's transitivity theorem (new entry). The spectral theorem is used for
`eq:ground-spectral-projection`.

The proposition is the only place where a prerequisite beyond the manuscript enters in an
essential way. In finite dimensions the transitivity theorem reduces to the statement that
an irreducible *-algebra of matrices is the full matrix algebra, and the proof says so.

## Mathematics added to the sources

Offered for review as additions.

1. `prop:pure-faithful-scalars`: a state that is both pure and faithful exists only on the
   scalars. So for a cyclic ground in a space of dimension at least two, at most one of the
   two properties of LG's open question can hold. LG asks the two questions separately and
   gives one example of a pure state that is not faithful.
2. The five-site examples: the site algebra gives a faithful state that is not pure,
   written as the average of five point states; the full matrix algebra gives a pure state
   that is not faithful; the cell algebra gives a faithful state that is not pure and is
   not cyclic.
3. An example where the abstract and concrete faithfulness differ: `M_2 ⊕ ℂ` on
   `ℂ² ⊗ ℂ²` with the maximally entangled vector.
4. The examples for the independence of uniqueness, cyclicity and a gap: two separate
   paths (cyclic and gapped, not unique) and the operator with eigenvalues `1/n` with the
   algebra of the unilateral shift (unique and cyclic, not gapped). LG states these
   failures without examples.
5. The uniqueness of the GNS intertwiner, not only its existence.
6. The invariance of the complement of the generated space and the location of the defect
   of the unequal-cell path: it is the complement for the cell algebra and lies inside the
   generated space for an algebra with a cyclic ground.

## Sentences that interpret

- "The vector state, the ray of y and the projection P_y describe the selected ground at
  three levels, and none of them fixes a phase of y." LG's diagram with arrows.
- "LG calls K⊥ the dark complement." The name is LG's.
- "Neither the algebra nor the spectrum selects a ray."
- "the examples show that the choice of algebra decides which one, if either."
- The last two sentences, on the expression "pure luminosity". LG has a warning to the
  same effect.

## Names

`(𝔐, π, ℋ, y)` for the pointed representation, `ω_y` for its vector state, and `𝔇 = π(𝔐)`
as in Section 5. `(π_ω, ℋ_ω, Ω_ω)` for GNS data and `U` for the intertwiner. `P_0` for the
spectral projection at the bottom of the spectrum. Section 4 uses `Π_0` for the aligned
graph projection, so `P_0` does not collide, but GPT may prefer another letter. The shift
in (c) has no letter, because `S` is the initial support. "Dark complement" is LG's name.

## Left for integration

- The three marker lines, each preceded by the comment `%%% PREVIEW MARKER`, are to be
  deleted, and the entry of `cyclic_ground_data_bibliography.tex` added to the
  bibliography.
- Part A refers forward to the new subsection for the exception concerning faithfulness.
- The fragment uses no label from the fidelity or ground-plane blocks, so it can be
  integrated before or after them.
- For a trial I inserted the three parts at the places named above in a copy of the
  manuscript at the base, with the new bibliography entry. It builds to 66 pages without
  warnings.
- `checks/fragments/test_cyclic_ground_data.py` has six tests. On this branch the suite
  gives 151 passed and 2 skipped, and 153 with `--full`. With the fidelity and ground-plane
  checks it would be 165 and 167.
- The check cites existing `UM` labels only. The fragment's seven labels can be added once
  they are in the manuscript.
