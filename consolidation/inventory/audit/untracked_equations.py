#!/usr/bin/env python3
"""Count labelled equations that sit outside every formal environment and proof and that
the ledger does not track as supplements, and list them for a range of source lines.

    python consolidation/inventory/audit/untracked_equations.py              # summary
    python consolidation/inventory/audit/untracked_equations.py LG 7181 7836 # list a range

Citation counts say how often the source itself refers back to an equation.  They are a
weak guide: the pair and orbital sections were almost entirely untracked while few of their
equations were cited more than twice.
"""
import collections
import re
import sys

from _common import SOURCES, read_csv


def scan(source, path, formal, tracked):
    text = path.read_text(encoding="utf-8")
    depth, outside = 0, {}
    for number, line in enumerate(text.splitlines(), 1):
        for match in re.finditer(r"\\(begin|end)\{([a-zA-Z*]+)\}", line):
            if match.group(2) in formal:
                depth += 1 if match.group(1) == "begin" else -1
        if depth == 0:
            for label in re.findall(r"\\label\{(eq:[^}]+)\}", line):
                outside[label] = number
    cites = collections.Counter()
    for match in re.finditer(r"\\(?:eqref|cref|Cref|ref)\{([^}]+)\}", text):
        for label in match.group(1).split(","):
            cites[label.strip()] += 1
    untracked = {label: number for label, number in outside.items() if (source, label) not in tracked}
    return outside, untracked, cites


def main():
    inventory, ledger = read_csv("statements.csv"), read_csv("migration.csv")
    formal = {entry["kind"] for entry in inventory} | {"proof"}
    tracked = {(row["source"], row["source_key"]) for row in ledger}
    if len(sys.argv) == 4:
        source, low, high = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
        _, untracked, cites = scan(source, SOURCES[source], formal, tracked)
        for label, number in sorted(untracked.items(), key=lambda item: item[1]):
            if low <= number < high:
                print(f"{source} {number:5d}  {label}  (cited {cites[label]} time(s))")
        return
    for source, path in SOURCES.items():
        outside, untracked, cites = scan(source, path, formal, tracked)
        print(f"{source}: {len(outside)} labelled equations outside formal environments and proofs; "
              f"{len(outside) - len(untracked)} tracked as supplements; {len(untracked)} untracked")
        for threshold in (1, 2, 3):
            print(f"   untracked and cited at least {threshold} time(s): {sum(1 for l in untracked if cites[l] >= threshold)}")


if __name__ == "__main__":
    main()
