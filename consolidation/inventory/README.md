# Statement inventory for the TAUM–LG consolidation

Working input for a merged outline. It lists every formal statement in the two manuscripts so that each result can be
given one home, with what it needs and what it contributes. The manuscripts are not modified by anything here.

The first proposed organization and section-level migration map are in [the consolidation outline](../OUTLINE.md).
That map is not a completed proof-by-proof dependency audit.

```sh
python consolidation/inventory/inventory.py     # from the repository root; rewrites statements.csv
```

Standard library only.

| File | What it is |
| --- | --- |
| `inventory.py` | Reads both manuscript sources and writes `statements.csv` |
| `statements.csv` | Generated: one row per theorem-like environment (281 rows, including model choices and doctrinal identifications) |
| `annotations.csv` | Written by hand: the reviewed columns, keyed by manuscript and label |

## Columns

**Extracted from the source.** `manuscript`, `line`, `section`, `kind`, `title`, `label`, `status`, `has_proof`.
`label` is the first label inside the environment; a few definitions carry only an equation label. `status` is LG's
explicit status label where one follows the statement.

**Automated prompts.** These say where to look. They are not dependency results.

- `auto_refers_to`: statements whose labels this statement or its proof cites.
- `auto_cites_axiom` (LG): whether one of LG's five axioms is cited directly or through other statements.
- `auto_ground_objects`: whether ground vectors, ground projections or ground-plane objects appear in the statement, in
  its proof, or in a statement it refers to. Matched by macro, by raw notation and by the words ground, vacuum, luminous.
- `auto_taum_cited_nearby` (LG): TAUM is cited within a few lines.

**Reviewed.** Filled only from `annotations.csv`, by someone who has read the statement and its proof.

- `proof_minimally_needs`: the hypotheses the proof actually uses.
- `realization_supplies`: how the realization provides those objects.
- `interpretation`: the reading the realization gives the result.
- `basis`: what was read, and by whom the point was made if not the annotator.

## What the automated columns can and cannot tell you

Two examples show both directions of error.

- *A proof can need less than the surrounding axioms.* `thm:ground-graphs-neutral` sits among statements that use both
  grounds, and its `auto_ground_objects` is "statement". Its proof uses only that the selected map between the two lines
  is unitary; uniqueness, cyclicity and energy of the ground vectors play no part. The flag records where the objects
  come from in LG, not what the theorem requires.
- *A pattern can miss objects written another way.* An earlier version of the ground-object pattern matched macros only
  and reported the momentary-continua section as free of ground objects. That section writes its luminous ground
  projections in raw notation. The pattern now includes raw notation and the three words above, and can still miss
  others.

Citations are a weaker guide still. Of LG's 170 statements other than its axioms, 16 cite an axiom directly or through
other statements, while 90 show ground objects. A missing citation is not a missing hypothesis: standing assumptions
and earlier definitions legitimately carry it. The point for the consolidated text is to make each result's hypotheses
identifiable without repeating the axiom system in every theorem.

Neither flag measures dependence on injectivity or on strict excess.

## Current counts

| | Statements | With a proof | With an explicit status label |
| --- | --- | --- | --- |
| TAUM | 106 | 52 | 0 (TAUM uses a summary table) |
| LG | 175 | 71 | 14 |

LG's total includes 29 doctrinal identifications and 17 model choices. Both kinds are part of the consolidation;
excluding them would omit formal statements of the interpretive commitments and selected models. The extractor
refuses any used, declared theorem-like environment that is not in its list of supported kinds.

Six statements are reviewed so far. Add a row to `annotations.csv` after reading a proof; the script refuses annotations
whose label does not exist.

Thirteen TAUM remarks have no label. Their temporary annotation keys are `@` followed by the source line
number (for example, `@100`); the CSV retains a blank label and a separate line column. Such keys are not stable across
source edits. Assign or reconcile stable labels when those entries are migrated, rather than treating line numbers as
permanent identifiers.
