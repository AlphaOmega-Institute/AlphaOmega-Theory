# Proposed ledger rows

Rows written alongside a fragment in `manuscripts/consolidated/fragments/`, in the
columns of `../migration.csv`. Their `draft_labels` and `draft_quote` point into the
fragment, so `../check_migration.py` would reject them. They wait here until the text is
integrated and has its final labels.

    python consolidation/inventory/proposed/build_moving_subspace_rows.py

The script reads the current rows for the fragment's source lines, applies the proposed
updates, keeps rows that are already explicit, appends rows for source passages that had
none, writes `moving_subspace_rows.csv`, and applies the validator's rules with the
fragment as the draft passage. The fragment's notes give the change in counts and the
rulings these rows leave open.
