# Validation of the code-and-diagram guides

Validation date: 19 October 2026.

## Content checks

- All five notebooks pass `nbformat.validate`.
- Every displayed Mathematica block matches its recorded contiguous range in the original decoded input cell.
- The two complete source transcriptions still match the original `.nb` input cells and raw `BoxData`.
- Original `.nb` and `.wdx` file hashes are unchanged.
- All 30 diagram regions are included in the notebooks. Their pixel arrays are identical to their recorded rectangles in the archived source PNGs.
- Diagram regions retain the plotted axes, tick labels, curves, colours, original plot titles and mesh details. Surrounding document text and page numbers are excluded.
- Notebook text contains no document comparisons, chapter references, original-document links or full-page-image references.
- All notebook cells are reading material; no synthetic execution output is supplied.

The extraction contact sheet was visually reviewed to check that captions and surrounding text were removed and plot labels were retained. Explanations distinguish sorted-state index, spatial position and energy axes; note missing colour scales; and describe the reversed height scales of the relevant surfaces.

## Reproduce the checks

```bash
python tools/verify_sources.py
```

This is a provenance and presentation check. Mathematica was not rerun, no matrix entries were changed, and no numerical curves were regenerated. The archived diagram inputs are not all present as active plotting routines in the saved programs; the guides say so beside the relevant code.
