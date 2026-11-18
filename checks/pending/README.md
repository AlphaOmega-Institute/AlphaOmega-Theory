# Pending checks

These checks cover results agreed during review and awaiting a full treatment in the
consolidated manuscript. A file can cite a `UM` label for a case already present. Move
each file beside the other checks when its remaining results enter the manuscript.

```bash
python checks/run_all.py -k pending/             # about one second
```

| File | Belongs to | Statements |
| --- | --- | --- |
| `test_port_calibration.py` | Signed response (outline III.3); two-level witness at UM `eq:passive-signed-witness` | The smallest product of port norms with which an ordinary unitary group reproduces a response with line weights `w_E` is `Σ|w_E|`. The peak amplitude is only a lower bound for it. The two-level generator with the opposite return sign: energies 1 and 2 and weights 2 and −1 at coupling √2, an exceptional point at 3/2, complex energies at 2. |
| `test_recollection_relations.py` | Part V, planned entry "Grasper–grasped frame and recollection" | The relation between rays; transport of effects through an isometry; the mixing equality under `G′W = WG` and its failure without it; the `w(1−w)` bound for an aspect of a paired state and the cancellation of its two outside parts. |
