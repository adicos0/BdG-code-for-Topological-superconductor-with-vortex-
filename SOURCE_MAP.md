# Code and diagram provenance

This file records implementation provenance for maintainers. The notebooks themselves explain the code and diagrams directly.

## Code

The authoritative calculation files are `BdG topological.nb` and `added vortex.nb`, unchanged from commit `d88d12cded1a954d53c2d36d515525c8100ded1c`. Original Wolfram data files are also unchanged.

| Input cell | BdG topological.nb line | added vortex.nb line | Content |
| --- | ---: | ---: | --- |
| 1 | 23 | 23 | Parameters, hopping dictionaries, 20 × 20 lattice, periodic matrices and `Hop1x`. |
| 2 | 922 | 948 | Onsite pairing construction `SCs`. |
| 3 | 1049 | 1075 | Occupation function, coupling, target density and initialization. |
| 4 | 1398 | 1437 | Commented alternative matrix construction. |
| 5 | 1441 | 1480 | Matrix assembly, eigensystem, occupations, pairing update, commented vortex block and exports. |
| 6 | 4032 | 45843 | Original sorting expression. |
| 7 | 4073 | 45884 | Original spectrum plotting command. |

The `.wl` transcriptions preserve the original expressions. The paired JSON files retain raw `BoxData`, line locations and decoded text. Notebook metadata records the exact source cell and character range of each displayed code block. New explanatory prose describes the operations in those blocks; it does not change their expressions or parameters.

## Diagram extraction

`diagrams/manifest.json` records the source PNG, rectangle and SHA-256 hash for each of the 30 diagrams. Rectangles use pixel coordinates `(left, top, right, bottom)` in the archived page render. `tools/extract_diagrams.py` copies those regions without resampling, filtering, recolouring, adding labels or regenerating numerical data.

The existing `reference_pages/` directory is retained only as an internal source-asset archive for reproducible extraction. Its manifest records the original document identity, file hash and page mapping. These pages are not displayed or linked in the reading notebooks. The notebook image URLs point only to `diagrams/`.

The extracted regions exclude captions, equations, surrounding paragraphs and page numbers. Original plot-internal labels, tick marks and titles remain. All previously displayed numerical figures are represented among the extracted regions, together with the lattice diagram.

## Scope of the explanations

- Matrix sizes, offsets, iteration order and signs are explained from the active source expressions.
- A commented term is described as inactive. Missing LDOS and impurity plotting routines are not replaced with invented code.
- Diagram explanations identify the actual displayed axes. The site-index LDOS panels are not relabelled as energy plots.
- Heatmaps without a colour bar are not assigned numerical colour values.
- The reversed vertical scales in several surfaces are explained explicitly; visual upward direction is not treated as numerical increase.
- Energy levels near zero are described as energy observations, without assigning an unsupported physical identity.
