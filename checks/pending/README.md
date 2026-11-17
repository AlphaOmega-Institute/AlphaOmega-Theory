# Pending checks

Checks of results that were worked out and agreed during review but are not yet in the
consolidated manuscript. They cite no `UM` label, because the text they belong to does
not exist yet. When that text is written, move each file beside the other checks and add
its labels.

```bash
python checks/run_all.py -k pending/             # about one second
```

| File | Belongs to | Statements |
| --- | --- | --- |
| `test_port_calibration.py` | Signed response (outline III.3; working title §9) | The smallest product of port norms with which an ordinary unitary group reproduces a response with line weights `w_E` is `Σ|w_E|`. The peak amplitude is only a lower bound for it. The two-level generator with the opposite return sign: energies 1 and 2 and weights 2 and −1 at coupling √2, an exceptional point at 3/2, complex energies at 2. |
| `test_recollection_relations.py` | Part V, planned entry "Grasper–grasped frame and recollection" | The relation between rays; transport of effects through an isometry; the mixing equality under `G′W = WG` and its failure without it; the `w(1−w)` bound for an aspect of a paired state and the cancellation of its two outside parts. |
| `test_deconditioning_owed_clauses.py` | Ledger rows for LG `thm:canonical-deconditioning-semigroup` | The clauses the draft does not yet state (trace-norm estimate with constant 2, generator, continuity at zero) and the three it leaves one step away (fixed ground, exact decay, unique fixed state). The constant 2 is attained at time zero by pure states. |
