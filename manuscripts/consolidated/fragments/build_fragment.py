#!/usr/bin/env python3
"""Build a standalone PDF of a fragment, using the preamble of the consolidated manuscript.

    python manuscripts/consolidated/fragments/build_fragment.py NAME

The wrapper takes every line of Universal_Model.tex before \\begin{document}, so the
fragment is typeset with the same packages, environments and macros. Labels that live
in the main file are read from an auxiliary run of the main file, so cross-references
print their real numbers. They are printed without a hyperlink, because their targets
are not in the standalone PDF. The section counter continues from the last numbered
section before the appendices. Equation numbering continues from the end of the main
file. Nothing is written into Universal_Model.tex.

The script reports LaTeX warnings, bad boxes and pdfTeX warnings such as links to
missing destinations.
"""
import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAIN = HERE.parent / "Universal_Model.tex"


def run(command, cwd):
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True)
    if result.returncode:
        sys.stderr.write(result.stdout[-3000:])
        raise SystemExit(f"{' '.join(command)} failed")


def main():
    parser = argparse.ArgumentParser(description="Build a PDF for a consolidated-manuscript fragment.")
    parser.add_argument("name", help="stem of the fragment's .tex file in this directory")
    name = parser.parse_args().name
    if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
        parser.error("supply a fragment name without a path or extension")
    source = HERE / f"{name}.tex"
    if not source.is_file():
        parser.error(f"fragment does not exist: {source.name}")
    body = source.read_text(encoding="utf-8")
    main_text = MAIN.read_text(encoding="utf-8")
    preamble = main_text.split("\\begin{document}")[0]
    extra_bib = (HERE / f"{name}_bibliography.tex")
    bibliography = re.search(r"\\begin\{thebibliography\}.*?\\end\{thebibliography\}", main_text, flags=re.S).group(0)
    if extra_bib.exists():
        bibliography = bibliography.replace("\\end{thebibliography}", extra_bib.read_text(encoding="utf-8") + "\n\\end{thebibliography}")
    main_sections = re.split(r"^[ \t]*\\appendix\b", main_text, maxsplit=1, flags=re.M)[0]
    sections = len(re.findall(r"^[ \t]*\\section(?=\{|\[)", main_sections, flags=re.M))
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        # a copy of the main file that also reports its last equation number
        (tmp / "main.tex").write_text(
            main_text.replace("\\end{document}", "\\typeout{LASTEQUATION=\\arabic{equation}}\n\\end{document}"),
            encoding="utf-8")
        for _ in range(2):
            run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "main.tex"], tmp)
        # labels of the main file, made available to the wrapper as external references
        aux = (tmp / "main.aux").read_text(encoding="utf-8")
        labels = "\n".join(line for line in aux.splitlines() if line.startswith(r"\newlabel{"))
        external_keys = re.findall(r"^\\newlabel\{([^{}]+)\}", labels, flags=re.M)
        markers = "\n".join(r"\expandafter\def\csname fragmentexternal@" + key + r"\endcsname{}"
                            for key in external_keys)
        # amsmath's eqref also calls ref. Retain ordinary local links and print
        # external references using hyperref's explicit, unlinked form.
        reference_setup = r"""
\NewCommandCopy{\FragmentOriginalRef}{\ref}
\RenewDocumentCommand{\ref}{s m}{%
  \ifcsname fragmentexternal@#2\endcsname
    \FragmentOriginalRef*{#2}%
  \else
    \IfBooleanTF{#1}{\FragmentOriginalRef*{#2}}{\FragmentOriginalRef{#2}}%
  \fi}
"""
        main_log = (tmp / "main.log").read_text(encoding="utf-8", errors="replace")
        last_equation = int(re.search(r"LASTEQUATION=(\d+)", main_log).group(1))
        wrapper = (preamble
                   + "\\makeatletter\n" + labels + "\n" + markers + "\n\\makeatother\n"
                   + "\\begin{document}\n" + reference_setup
                   + f"\\noindent\\textbf{{Fragment for review.}} Section and equation numbers are provisional. "
                     f"References to Sections 1--{sections} are to the consolidated manuscript.\\par\\medskip\n"
                   + f"\\setcounter{{section}}{{{sections}}}\\setcounter{{equation}}{{{last_equation}}}\n"
                   + body
                   + "\n\\begingroup\\small\n" + bibliography + "\n\\endgroup\n\\end{document}\n")
        (tmp / f"{name}.tex").write_text(wrapper, encoding="utf-8")
        for _ in range(2):
            run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"{name}.tex"], tmp)
        log = (tmp / f"{name}.log").read_text(encoding="utf-8", errors="replace")
        problems = [line for line in log.splitlines()
                    if re.match(r"(?:(?:LaTeX(?: Font)?|Package \S+|Class \S+|pdfTeX) warning|Overfull|Underfull)", line, re.I)]
        shutil.copy(tmp / f"{name}.pdf", HERE / f"{name}.pdf")
        pages = re.search(r"Output written on .*?\((\d+) pages?\b", log)
        print(f"{name}.pdf: {pages.group(1) if pages else '?'} pages; warnings: {len(problems)}")
        for line in problems[:20]:
            print("   ", line)


if __name__ == "__main__":
    main()
