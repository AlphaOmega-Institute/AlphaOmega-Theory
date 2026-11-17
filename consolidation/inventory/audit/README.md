# Ledger audit scripts

Three screens used while reviewing `migration.csv`. Each prints a list for a person to
read; none of them decides anything. Run them from the repository root.

```bash
python consolidation/inventory/audit/multipart_single_rows.py
python consolidation/inventory/audit/dependency_order.py
python consolidation/inventory/audit/untracked_equations.py
python consolidation/inventory/audit/untracked_equations.py LG 7181 7836
```

| Script | Question | Result on the ledger at `1419df8` |
| --- | --- | --- |
| `multipart_single_rows.py` | Which theorem-level statements have one row marked explicit or recoverable although the source statement is enumerated, carries several labelled equations, or is long? | 14 entries. Each was read against the source and the draft and found correct: six explicit, eight recoverable. |
| `dependency_order.py` | Which statements with something owed do other such statements cite, and how many cite none? | TAUM: 54 of 66 cite no other owed entry; LG: 61 of 68. The passive-response section holds the most ready entries (17). |
| `untracked_equations.py` | Which labelled equations outside formal environments and proofs does the ledger not track? | LG 215 of 261, TAUM 57 of 69. Before the pair and orbital supplements were added, this screen showed 29 untracked equations in LG lines 7181–7836; it now shows none there. |

Limits. The first screen is a heuristic for where to read and can miss a short multi-part
statement. The second uses the inventory's automatic cross-references between formal
environments, so it cannot see a definition made in running prose. The third counts
labels, and an unlabelled displayed equation or a prose argument is invisible to it;
`heading_queue.py` covers passages by heading.
