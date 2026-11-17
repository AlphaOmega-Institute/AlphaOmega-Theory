# Reviewing a manuscript fragment

Put the section body in `NAME.tex` beside the builder. Use the main manuscript’s preamble, give new statements and equations their own labels, and cite existing results with `\ref` or `\eqref`. Place any additional bibliography entries in `NAME_bibliography.tex`.

From the repository root, use pdfLaTeX to build the preview.

```sh
python manuscripts/consolidated/fragments/build_fragment.py NAME
```

The preview continues the section and equation numbers from the current manuscript. References to the main manuscript print their numbers. Links within the fragment remain active. Check the final PDF and the warning summary before sending it for review.
