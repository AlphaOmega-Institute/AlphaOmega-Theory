#!/usr/bin/env python3
"""Rank source statements that still have something owed by how many other such statements
cite them, and count those whose cited prerequisites are all already present.

The citations are the inventory's `auto_refers_to` field: references between formal
environments.  A definition made in running prose is invisible to it, so a statement
listed as ready may still need notation from an untracked passage.

    python consolidation/inventory/audit/dependency_order.py
"""
import collections
import re

from _common import read_csv


def main():
    inventory, ledger = read_csv("statements.csv"), read_csv("migration.csv")
    coverage = collections.defaultdict(set)
    for row in ledger:
        coverage[(row["source"], row["source_key"])].add(row["coverage"])
    for source in ("TAUM", "LG"):
        entries = {e["label"]: e for e in inventory if e["manuscript"] == source and e["label"]}
        owed = {label for label in entries if "owed" in coverage[(source, label)]}
        cites = {label: {x for x in re.split(r"[;,\s]+", entries[label]["auto_refers_to"]) if x in entries and x != label}
                 for label in entries}
        cited_by = collections.defaultdict(set)
        for label, targets in cites.items():
            for target in targets:
                cited_by[target].add(label)

        def dependants(label):
            seen, stack = set(), [label]
            while stack:
                for other in cited_by[stack.pop()]:
                    if other not in seen:
                        seen.add(other)
                        stack.append(other)
            return seen

        ranking = sorted(((len(dependants(label) & owed), len(cited_by[label] & owed), label) for label in owed), reverse=True)
        ready = [label for label in owed if not (cites[label] & owed)]
        print(f"{source}: {len(entries)} labelled entries, {len(owed)} with something owed, "
              f"{len(ready)} of those cite no other owed entry")
        print("   owed entries most depended upon (all dependants, direct dependants, label, line, section):")
        for total, direct, label in ranking[:10]:
            if total == 0:
                break
            print(f"   {total:3d} {direct:3d}  {label:40s} {entries[label]['line']:>5}  {entries[label]['section'][:50]}")
        by_section = collections.Counter(entries[label]["section"] for label in ready)
        print("   sections holding the most ready entries:")
        for section, count in by_section.most_common(6):
            print(f"   {count:3d}  {section[:70]}")


if __name__ == "__main__":
    main()
