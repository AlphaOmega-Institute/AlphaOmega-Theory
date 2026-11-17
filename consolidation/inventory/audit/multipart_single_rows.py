#!/usr/bin/env python3
"""List theorem-level source statements that the ledger records with a single row marked
explicit or recoverable although the statement is enumerated, carries several labelled
equations, or is long.  Such a row may hide clauses; each needs reading against the source.

    python consolidation/inventory/audit/multipart_single_rows.py
"""
import re

from _common import THEOREM_LEVEL, entry_key, read_csv, rows_by_entry, source_lines, statement


def main():
    inventory, ledger = read_csv("statements.csv"), read_csv("migration.csv")
    lines, grouped = source_lines(), rows_by_entry(ledger)
    flagged = []
    for entry in inventory:
        if entry["kind"] not in THEOREM_LEVEL:
            continue
        rows = grouped.get(entry_key(entry), [])
        if len(rows) != 1 or rows[0]["coverage"] == "owed":
            continue
        text = statement(lines[entry["manuscript"]], entry)
        items = sum(line.count("\\item") for line in text)
        equations = sum(len(re.findall(r"\\label\{eq:", line)) for line in text)
        if items >= 2 or equations >= 3 or len(text) >= 22:
            flagged.append((entry["manuscript"], int(entry["line"]), entry["kind"], entry_key(entry)[1],
                            rows[0]["coverage"], items, equations, len(text)))
    print(f"{len(flagged)} single-row, non-owed theorem-level entries with a long or multi-part source statement")
    print("source  line  kind         key                                        coverage     items  eqs  lines")
    for source, line, kind, key, coverage, items, equations, length in sorted(flagged):
        print(f"{source:6s} {line:5d}  {kind:12s} {key:42s} {coverage:12s} {items:5d} {equations:4d} {length:6d}")


if __name__ == "__main__":
    main()
