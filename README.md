# AlphaOmega-Theory

Public work in progress for the AlphaOmega Institute's theoretical research.

The files on `main` are the current shared drafts. They can be edited as the research develops; selected versions may later be released and archived with DOIs.

## Current working manuscript

**Toward a Universal Model: Perception as a correspondence between state spaces**

- [Read the PDF](manuscripts/consolidated/Universal_Model.pdf)
- [Self-contained TeX source](manuscripts/consolidated/Universal_Model.tex)

This is the working consolidation of LG and TAUM into one paper, and it is not yet complete.

## Source manuscripts

LG and TAUM are the two manuscripts being consolidated. They remain the reference for every argument the working manuscript has not yet carried over, and the consolidation does not edit them. The [migration ledger](consolidation/inventory/migration.csv) records what has been carried over and what is still owed. Its [guide](consolidation/inventory/README.md#migration-into-the-consolidated-manuscript) explains how to read it, and the [outline](consolidation/OUTLINE.md) gives the planned structure.

| Manuscript | Self-contained TeX source |
| --- | --- |
| **LG — The Luminous Ground: Emptiness, Awareness, and the Algebra of Appearance** | [LG_Master_Unified.tex](manuscripts/lg/LG_Master_Unified.tex) |
| **TAUM — Formalizing Buddhist (Meta)physics: Toward a Universal Model** | [Toward_a_Universal_Model_Self_Contained.tex](manuscripts/taum/Toward_a_Universal_Model_Self_Contained.tex) |

The initial import on November 11, 2026 contains the November 10 manuscript revisions, unchanged. Later review suggestions have not been silently incorporated into that baseline. The manuscripts distinguish mathematical results, realization assumptions, interpretations, and open problems; their working-draft status applies throughout.

The November 12 housekeeping update clarifies finite-isometric correspondence versus evolution and the pair input Hamiltonian, adds the Gilbert and Kalman attribution to the finite-dimensional response discussion, corrects the TAUM abstract attribution, and updates LG's citation of the public TAUM working manuscript. The original axioms, results and numerical examples are unchanged. Both self-contained sources compile without warnings (LG: 112 pages; TAUM: 59 pages), the edited pages have been visually inspected, and all 41 numerical checks pass, including the full orbital convergence check.

Each manuscript includes its bibliography and is maintained as a single self-contained TeX file suitable for Overleaf. The initial import contains the sources; verification scripts and compiled release PDFs can be added separately.

### Building

Use a TeX Live or equivalent installation with pdfLaTeX, latexmk, KOMA-Script, and the packages named in the preamble.

From `manuscripts/lg/`:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error LG_Master_Unified.tex
```

From `manuscripts/taum/`:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error Toward_a_Universal_Model_Self_Contained.tex
```

From `manuscripts/consolidated/`:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error Universal_Model.tex
```

Three runs of `pdflatex` on the same file give the same result. The committed PDF is built from the committed TeX.

## Verification checks

`checks/` contains reproducible numerical checks of calculations in the working manuscript and in the two sources. Run `python checks/run_all.py` (NumPy and SciPy required). See [checks/README.md](checks/README.md) for what each check establishes and what it leaves open.

## Shared editing

- Fetch the latest `main` before editing and inspect any intervening changes.
- Use a separate branch for simultaneous or overlapping edits, then reconcile the changes before merging.
- Merge a branch into `main` once its changes have been reviewed, then delete it. `main` holds the current state of every manuscript.
- Make focused commits with messages describing the changes. Preserve other contributors' work and avoid force-pushing shared history.
- Keep each manuscript self-contained, including its bibliography. Do not replace sections with dependencies on external companion notes or source files.
- Synchronize Overleaf and Library copies deliberately with a known repository commit. Do not overwrite newer repository content with an older export.

Git records the working history. A dated release is a separate publication step, and can be archived on Zenodo when appropriate.

## Existing formalization

The `PrimitiveSpectra/` directory, `PrimitiveSpectra.lean`, and the Lake/Lean configuration contain the existing PrimitiveSpectra formalization project. Its files and build configuration are preserved alongside the manuscripts.

## Related repositories

- [AlphaOmega-LQCD](https://github.com/AlphaOmega-Institute/AlphaOmega-LQCD): lattice engine and associated computations.
- [magenta-gap](https://github.com/AlphaOmega-Institute/magenta-gap): color-model calculations and supporting material.

