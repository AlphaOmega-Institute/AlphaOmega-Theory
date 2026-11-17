# Proposed ledger rows

Rows written alongside a fragment in `manuscripts/consolidated/fragments/`, in the
columns of `../migration.csv`. Their `draft_labels` and `draft_quote` point into the
fragment, so `../check_migration.py` would reject them. They wait here until the text
is integrated and has its final labels.

    python consolidation/inventory/proposed/build_passive_response_rows.py

The script reads the current rows for the fragment's source lines, applies the proposed
updates, appends the clauses the port separated out, writes
`passive_response_rows.csv`, and applies the validator's rules with the fragment as the
draft passage: allowed statuses, no `absent` on a covered clause, every label present in
the fragment or the draft, a quote wherever only section labels are given, every quote
found exactly once, and no fragment label shared with the draft or either source.

`manuscripts/consolidated/fragments/README.md` gives the change in counts and the
decisions these rows leave open.
