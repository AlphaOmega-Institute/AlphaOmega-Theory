"""List source heading bodies that have no direct entry in migration.csv."""
import argparse
import bisect
import csv
import io
import re

from check_migration import HERE, LABEL, SOURCES


HEADING = re.compile(r"^[ \t]*\\((?:sub){0,2}section)\*?\s*\{", re.MULTILINE)
BRACE = re.compile(r"\\[{}\\]|[{}]")
DISPLAY = re.compile(r"\\\[|\\begin\{(?:equation|align|gather|multline|displaymath)\*?\}")
# Skip a literal percent or backslash while stripping TeX comments; preserve all line breaks.
COMMENT = re.compile(r"\\[%\\]|%[^\n]*")
NAMED = {"LG": {"Ground actualization", "Memory, recognition, and concatenation"},
         "TAUM": {"Fundamental decomposition", "Domain bookkeeping"}}
FIELDS = ["source", "start_line", "end_line", "heading_path", "display_blocks", "source_labels", "priority"]


def title_end(text, start):
    """Find the closing brace, including nested TeX in a title."""
    depth = 1
    for token in BRACE.finditer(text, start):
        if token[0] == "{":
            depth += 1
        elif token[0] == "}":
            depth -= 1
            if depth == 0:
                return token.start()
    raise ValueError("Unclosed source heading")


def heading_bodies(source):
    text = COMMENT.sub(lambda m: "" if m[0].startswith("%") else m[0], source.read_text(encoding="utf-8"))
    begin = text.index(r"\begin{document}") + len(r"\begin{document}")
    end = text.index(r"\begin{thebibliography}", begin)
    headings = list(HEADING.finditer(text, begin, end))
    ancestors = []
    for i, heading in enumerate(headings):
        start = heading.end()
        title_stop = title_end(text, start)
        title = " ".join(text[start:title_stop].split())
        level = heading[1].count("sub")
        ancestors[level:] = [title]
        next_heading = headings[i + 1].start() if i + 1 < len(headings) else end
        body = text[title_stop + 1:next_heading]
        # A label-only parent has no reading of its own. Its children are tested separately.
        if LABEL.sub("", body).replace(r"\appendix", "").strip():
            start_line = text.count("\n", 0, heading.start()) + 1
            end_line = text.count("\n", 0, next_heading)
            yield start_line, end_line, tuple(ancestors), body


def make_queue():
    tracked = {source: set() for source in SOURCES}
    with (HERE / "migration.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["source"] in tracked:
                tracked[row["source"]].add(int(row["source_line"]))
    output = io.StringIO()
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(FIELDS)
    for source, path in SOURCES.items():
        positions = sorted(tracked[source])
        for first, last, ancestry, body in heading_bodies(path):
            pos = bisect.bisect_left(positions, first)
            if pos < len(positions) and positions[pos] <= last:
                continue
            priority = "named" if NAMED[source].intersection(ancestry) else "standard"
            writer.writerow([source, first, last, " > ".join(ancestry), len(DISPLAY.findall(body)),
                             " ".join(LABEL.findall(body)), priority])
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="detect a stale queue without writing it")
    args = parser.parse_args()
    path = HERE / "heading_queue.csv"
    expected = make_queue().encode("utf-8")
    if args.check:
        if not path.exists() or path.read_bytes() != expected:
            parser.exit(1, "Heading queue is stale. Regenerate it with heading_queue.py.\n")
    else:
        path.write_bytes(expected)
    print(f"{len(expected.splitlines()) - 1} heading bodies in the reading queue")


if __name__ == "__main__":
    main()
