# Proposed ledger rows

Rows written alongside a fragment in `manuscripts/consolidated/fragments/`, in the
columns of `../migration.csv`. Their `draft_labels` and `draft_quote` point into the
fragment, so `../check_migration.py` would reject them. They wait here until the text is
integrated and has its final labels.

    python consolidation/inventory/proposed/build_symmetry_average_rows.py

The script takes the current rows for the clauses it names, wherever they sit in the
source, applies the proposed updates, appends a row for a source passage that had none,
writes `symmetry_average_rows.csv`, and applies the validator's rules with the fragment
as the draft passage. The fragment's notes give the change in counts.
