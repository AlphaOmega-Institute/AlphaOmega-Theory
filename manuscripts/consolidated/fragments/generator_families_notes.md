# `generator_families.tex`: families of generators, words, content equivalence, and a family that keeps the ground

Sources: LG lines 4835 to 5045 ("Vāsanā as conditioning data", "Words as conditioned
histories", "Profile-relative differentiated appearance", "The exhausted profile"), lines
5200 to 5235 ("Affective and perspectival variation"), and lines 6300 to 6345 ("The
exhausted profile representation is one-dimensional"). Base: `2a033c5`. Planned home: one
subsection at the end of Section 5, after "Uniqueness, cyclicity, purity and faithfulness
of a ground". The preview numbers it 13.1. In place it is 5.4.

    python manuscripts/consolidated/fragments/build_fragment.py generator_families
    python checks/run_all.py -k generator_families
    python consolidation/inventory/proposed/build_generator_families_rows.py

Nothing here changes `Universal_Model.tex` or the shared ledger. The scope is the one
agreed for this block: word order on the ground, equality of states on a chosen content
algebra, and the one-dimensional ground representation of the exhausted profile, with the
definition of the generated space that the last needs.

## Source results a reader can now find

| Source (LG) | Fragment |
|---|---|
| `def:phenomenally-complete-profile` | `eq:family-generated-space` and the sentence after it |
| the profile-relative appearance space and its inclusion (no separate row) | second identity of `eq:family-generated-space` |
| `prop:word-order-sensitivity`, the mathematical clause | `eq:word-order` and the paragraph after it |
| `def:content-equivalent-states` | `eq:family-content-equivalence` and the paragraph after it |
| `q:content-affect-algebra-choice` | last two sentences of that paragraph |
| `choice:exhausted-profile-ground-triviality` | the paragraph before the proposition |
| `prop:exhausted-words-ground-preserve` | `prop:ground-keeping-family` |
| `thm:exhausted-profile-character` | `eq:ground-character` and the paragraph after the proof |

## Change in the clause count

`consolidation/inventory/proposed/generator_families_rows.csv` proposes seven clauses of
LG as explicit: two recoverable (`prop:word-order-sensitivity`,
`def:content-equivalent-states`) and five owed (`def:phenomenally-complete-profile`,
`q:content-affect-algebra-choice`, `choice:exhausted-profile-ground-triviality`,
`prop:exhausted-words-ground-preserve`, `thm:exhausted-profile-character`). No row is
added.

The ledger goes from 497 rows (388 explicit, 23 recoverable, 86 owed) to 497 (395, 21,
81). LG goes from 240, 8, 57 to 247, 6, 52, and TAUM is unchanged. The six recoverable LG
rows left are the odd-extension group.

Clauses near these that the fragment does not carry:

- `choice:vasana-profiles`, the disposition profile as the index of a family. The text
  uses families of operators and says that LG indexes them by profiles. The profile
  itself, its strength or form, and the choice of families are part of the cognitive
  model.
- The second clause of `prop:word-order-sensitivity` (dependence of the cognitive result
  on the order), `id:generator-doctrinal-source`, the identification in
  `id:akara-differentiated-luminosity`, `choice:apoha-filter-superimposition` and
  `choice:tracks-not-sectors`.
- The duality-frame group (`def:derived-duality-frame`, `thm:duality-frame-vanishes-ground`,
  `choice:duality-generator`, `cor:no-independent-subject-ket`,
  `q:projection-or-effect-duality-frame` and the two related model choices). It was not
  part of the agreed scope. `thm:duality-frame-vanishes-ground` and
  `cor:no-independent-subject-ket` have mathematical content that could follow in the
  same style if wanted.
- `id:third-yoga` and the other readings of the exhausted profile.

## What the proofs need

From the consolidated manuscript: Proposition 5.1; Eq. `eq:ground-line-kept` with the
paragraph after it (Section 5.1); `prop:pure-faithful-scalars` (Section 5.3); the readout
relation and the states `ξ`, `η` of `eq:path-readout-states` (Section 5.2); the
unequal-cell range `D`; the path Hamiltonian; the recognition sequence of Section 6.

Used without proof: every element of a unital C*-algebra is a combination of four effects;
the GNS representation of a state is one-dimensional exactly when the state is a
character, and a one-dimensional representation is irreducible.

## Mathematics added to the sources

Offered for review as additions.

1. The centered vectors span the generated space minus the ground line for every family.
   LG states the inclusion and equality only when the ground is cyclic for the family.
2. Two orders that give different vectors can give one ray, as for anticommuting Pauli
   matrices, and different pure states can agree on every effect of an algebra. So a
   nonzero commutator on the ground establishes a vector difference and no more. This
   distinction is the reviewer's.
3. A unitary preserves the restriction of every state to the content algebra exactly when
   it commutes with it, and it can change the restriction to a second algebra. The two
   path states and a diagonal unitary are the example.
4. The proposition's statement that every element of the generated algebra commutes with
   the ground projection; the purity of the state on that algebra and its failure to be
   faithful unless the algebra is the scalars; the path Hamiltonian as a family that keeps
   the ground while the site algebra makes it cyclic; the example `|y⟩⟨w|` showing that the
   adjoint condition is needed for the *-algebra.
5. The examples of incomplete families whose union is complete, and the path Hamiltonian
   with one site projection generating every matrix.

## Sentences that interpret

- "LG generates appearances with specified families of operators."
- "LG indexes such families by disposition profiles. The choice of the families and their
  cognitive reading belong to the cognitive model and are not made here."
- "A claim that the order of two determinations changes what can be distinguished
  therefore needs the effects that distinguish the two results".
- "LG takes it, together with a second algebra 𝔉 for affective structure, inside the von
  Neumann algebra generated by 𝔇."
- "LG calls a profile with this property exhausted."

## Names

`𝖦 = {G_a}` for a family, because `𝒢_y` is the relaxation generator in Section 6.
`𝒲(𝖦)` for its unital *-algebra of words, `K_𝖦` for the generated space, `𝔇_𝖦` for the
norm closure. `𝔆` and `𝔉` for the content and affect algebras. `σ_x`, `σ_z` for the Pauli
matrices, because `X` denotes a determination. The label `eq:content-equivalence` is taken
in LG, so the equation is `eq:family-content-equivalence`.

## Left for integration

- The marker line, preceded by the comment `%%% PREVIEW MARKER`, is to be deleted.
- For a trial I inserted the subsection at the end of Section 5 in a copy of the
  manuscript at `2a033c5`. It builds to 75 pages without warnings.
- `checks/fragments/test_generator_families.py` has four tests. On this branch the suite
  gives 176 passed and 2 skipped, and 178 with `--full`.
- The check cites existing `UM` labels only. The fragment's six labels can be added once
  they are in the manuscript.
