# Provenance validation — source-only edition

Checked on 19 October 2026. This edition validates preservation of the supplied source material. It does not claim that the Mathematica calculations were rerun.

## Checks completed

- All five `.ipynb` files pass `nbformat.validate`.
- The notebooks contain **16 original code excerpts**. Each block exactly matches its recorded contiguous range in a decoded original input cell.
- All **14 nonempty calculation input cells** across the two original `.nb` files are transcribed into readable `.wl` files, with original `BoxData` and line information in paired JSON files.
- Plot-label strings, including escaped line continuations and embedded Wolfram formatting, remain intact in the transcription. Unknown executable box types cause the extractor to fail explicitly.
- **27 original thesis pages** are stored as unaltered full-page PNG renders, with **34 placements** across the five notebooks. No axes, curve values, equations, captions or colours have been edited. The PDF-to-printed-page mapping and page-image SHA-256 hashes are recorded in `reference_pages/manifest.json`.
- All page images were inspected in a contact sheet, with equation and figure pages additionally inspected at full size.
- Hashes verify that the original **two `.nb` and five `.wdx` files** are unchanged.
- Every notebook is a Markdown-only reading notebook. It contains no synthetic execution output or freshly generated numerical plot.

Run the verification with:

```bash
python tools/verify_sources.py
```

## What was not validated

Mathematica was not available for a native rerun. The current source inputs contain disabled vortex/impurity sections and do not supply an active workflow for all thesis figures. Consequently, the displayed thesis figures are identified as **archived source results**, not outputs produced by this revision.

The thesis's equations, comments, plot labels and interpretations are preserved as source material. They have not been silently corrected. The earlier Python tutorial's seven numerical tests and its convergence results do not apply to this source-only edition; the Python solver and its numerical plots have been removed from the current repository version.

## Rendering

The complete supplied PDF was downloaded through the user's connected Drive. Selected pages were rendered using Poppler `pdftoppm -r 120 -png`. The source has 74 PDF pages; printed page numbers for the selected pages equal PDF page numbers minus 14. The images retain the complete page layout, including figures, captions, equations, explanatory text and printed page numbers.
