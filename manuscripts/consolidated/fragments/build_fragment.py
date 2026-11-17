#!/usr/bin/env python3
"""Build a standalone PDF of a fragment, using the preamble of the consolidated manuscript.

    python manuscripts/consolidated/fragments/build_fragment.py passive_response

The wrapper takes every line of Universal_Model.tex before \\begin{document}, so the
fragment is typeset with the same packages, environments and macros. Labels that live
in the main file are read from an auxiliary run of the main file, so cross-references
print their real numbers. They are printed without a hyperlink, because their targets
are not in the standalone PDF. The fragment's section and equation counters continue
from the end of the main file. Nothing is written into Universal_Model.tex.

The script reports LaTeX warnings, bad boxes and pdfTeX warnings such as links to
missing destinations.
"""
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
    name = sys.argv[1] if len(sys.argv) > 1 else "passive_response"
    body = (HERE / f"{name}.tex").read_text(encoding="utf-8")
    main_text = MAIN.read_text(encoding="utf-8")
    preamble = main_text.split("\\begin{document}")[0]
    extra_bib = (HERE / f"{name}_bibliography.tex")
    bibliography = re.search(r"\\begin\{thebibliography\}.*?\\end\{thebibliography\}", main_text, flags=re.S).group(0)
    if extra_bib.exists():
        bibliography = bibliography.replace("\\end{thebibliography}", extra_bib.read_text(encoding="utf-8") + "\n\\end{thebibliography}")
    sections = len(re.findall(r"^\\section\{", main_text, flags=re.M))
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
        # drop each label's link anchor: the number is printed, and no link points outside this PDF
        labels = "\n".join(re.sub(r"\{[^{}]*\}\{\}\}$", "{}{}}", line)
                           for line in aux.splitlines() if line.startswith("\\newlabel{"))
        main_log = (tmp / "main.log").read_text(encoding="utf-8", errors="replace")
        last_equation = int(re.search(r"LASTEQUATION=(\d+)", main_log).group(1))
        wrapper = (preamble
                   + "\\makeatletter\n" + labels + "\n\\makeatother\n"
                   + "\\begin{document}\n"
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
                    if re.search(r"LaTeX Warning|Overfull|Underfull|undefined|multiply|pdfTeX warning", line)]
        shutil.copy(tmp / f"{name}.pdf", HERE / f"{name}.pdf")
        pages = re.search(r"Output written on .*?\((\d+) pages", log)
        print(f"{name}.pdf: {pages.group(1) if pages else '?'} pages; warnings: {len(problems)}")
        for line in problems[:20]:
            print("   ", line)


if __name__ == "__main__":
    main()
