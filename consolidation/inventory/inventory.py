#!/usr/bin/env python3
"""Inventory of the formal statements in TAUM and LG.

Reads the two manuscript sources, never writes to them, and produces
statements.csv: one row per theorem-like environment.

    python consolidation/inventory/inventory.py        # from the repository root

Three kinds of column:

  extracted   manuscript, line, section, kind, title, label, status, has_proof
  automated   refers_to, cites_axiom, ground_objects, taum_cited_nearby
              These are prompts for reading the proofs.  They are not
              dependency results: a statement can rely on an assumption
              through a standing convention or an earlier definition
              without citing it, and a text pattern can miss an object
              written in other notation.
  reviewed    proof_minimally_needs, realization_supplies, interpretation,
              basis.  Filled only from annotations.csv, by someone who has
              read the statement and its proof.
"""
import collections
import csv
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = {"TAUM": ROOT / "manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex",
       "LG": ROOT / "manuscripts/lg/LG_Master_Unified.tex"}
KINDS = ("theorem", "proposition", "lemma", "corollary", "definition", "axiom", "assumption",
         "example", "remark", "openquestion")
ENV = re.compile(r"\\begin\{(" + "|".join(KINDS) + r")\}(\[[^\]]*\])?(.*?)\\end\{\1\}", re.S)
REF = re.compile(r"\\(?:c|C)?(?:ref|eqref)\{([^}]*)\}")
LAB = re.compile(r"\\label\{([^}]*)\}")
# ground vectors, ground projections and ground-plane objects, by macro, by raw notation and by word
GROUND = re.compile(r"\\ket0|\\ket\{0\}|\\bra0|\\Omega_0|\\Kzero|\\Pzero[PF]|\\Nzero|\\Gzero[PF]|\\uzero|P_\{0|"
                    r"\bground\b|\bvacuum\b|\bluminous\b", re.I)


def parse(ms, path):
    text = path.read_text()
    start = text.index("\\begin{document}")
    body = text[start:text.index("\\begin{thebibliography}")]
    line0 = text.count("\n", 0, start)
    secs = [(m.start(), m.group(1)) for m in re.finditer(r"\\section\{([^}]*)\}", body)]
    out = []
    for m in ENV.finditer(body):
        inner, tail = m.group(3), body[m.end():m.end() + 8000]
        pm = re.match(r"\s*\\begin\{proof\}(.*?)\\end\{proof\}", tail, re.S)
        proof = pm.group(1) if pm else ""
        st = re.match(r"\s*\\statuslabel\{([^}]*)\}", tail[pm.end():] if pm else tail)
        labels = LAB.findall(inner) + LAB.findall(proof)
        near = body[max(0, m.start() - 600):m.end() + 600]
        out.append(dict(
            ms=ms, line=line0 + body.count("\n", 0, m.start()) + 1,
            section=next((n for p, n in reversed(secs) if p <= m.start()), ""),
            kind=m.group(1), title=re.sub(r"\s+", " ", (m.group(2) or "")[1:-1]),
            label=LAB.findall(inner)[0] if LAB.findall(inner) else "", labels=labels,
            status=st.group(1) if st else "", has_proof=bool(pm), statement=inner, proof=proof,
            taum_nearby=(ms == "LG" and "TAUM2026" in near)))
    return out


rows = [r for ms, p in SRC.items() for r in parse(ms, p)]
for i, r in enumerate(rows):
    r["key"] = (r["ms"], r["label"] or f"@{r['line']}")
owner = {(r["ms"], lab): r["key"] for r in rows for lab in r["labels"]}
graph = collections.defaultdict(set)
for r in rows:
    for group in REF.findall(r["statement"] + r["proof"]):
        for lab in group.split(","):
            tgt = owner.get((r["ms"], lab.strip()))
            if tgt and tgt != r["key"]:
                graph[r["key"]].add(tgt)


def closure(seed):
    """Statements that reach a member of `seed` by following references."""
    seen, changed = set(seed), True
    while changed:
        changed = False
        for r in rows:
            if r["key"] not in seen and graph[r["key"]] & seen:
                seen.add(r["key"]); changed = True
    return seen


axioms = {r["key"] for r in rows if r["kind"] == "axiom"}
via_axiom = closure(axioms)
g_stmt = {r["key"] for r in rows if GROUND.search(r["statement"])}
g_proof = {r["key"] for r in rows if GROUND.search(r["proof"])}
g_via = closure(g_stmt | g_proof)

notes = {}
ann = HERE / "annotations.csv"
if ann.exists():
    with open(ann, newline="") as fh:
        for a in csv.DictReader(fh):
            notes[(a["manuscript"], a["label"])] = a
known = {r["key"] for r in rows}
unknown = [k for k in notes if k not in known]
assert not unknown, f"annotations refer to statements that do not exist: {unknown}"

with open(HERE / "statements.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["manuscript", "line", "section", "kind", "title", "label", "status", "has_proof",
                "auto_refers_to", "auto_cites_axiom", "auto_ground_objects", "auto_taum_cited_nearby",
                "proof_minimally_needs", "realization_supplies", "interpretation", "basis"])
    for r in rows:
        k = r["key"]
        cites = ("is an axiom" if k in axioms else "directly" if graph[k] & axioms else
                 "through other statements" if k in via_axiom else "no")
        ground = ("statement" if k in g_stmt else "proof" if k in g_proof else
                  "through other statements" if k in g_via else "none found")
        a = notes.get(k, {})
        w.writerow([r["ms"], r["line"], r["section"], r["kind"], r["title"], r["label"], r["status"], r["has_proof"],
                    "; ".join(sorted(t[1] for t in graph[k])), cites if r["ms"] == "LG" else "", ground, r["taum_nearby"] or "",
                    a.get("proof_minimally_needs", ""), a.get("realization_supplies", ""), a.get("interpretation", ""),
                    a.get("basis", "")])

print(f"{len(rows)} formal statements; {len(notes)} reviewed.")
for ms in SRC:
    sub = [r for r in rows if r["ms"] == ms]
    kinds = collections.Counter(r["kind"] for r in sub)
    print(f"  {ms}: {len(sub)} statements, {sum(r['has_proof'] for r in sub)} with a proof, "
          f"{sum(bool(r['status']) for r in sub)} with an explicit status label; " +
          ", ".join(f"{n} {k}" for k, n in kinds.most_common()))
lg = [r for r in rows if r["ms"] == "LG" and r["key"] not in axioms]
print(f"  LG, excluding its {len(axioms)} axioms ({len(lg)} statements): "
      f"{sum(r['key'] in via_axiom for r in lg)} cite an axiom directly or through other statements; "
      f"{sum(r['key'] in g_via for r in lg)} show ground objects in the statement, the proof or a statement referred to.")
