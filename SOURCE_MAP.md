# Original-source map

## Authoritative sources

- **Thesis:** Aditya Banerjee, *Understanding the Magnetic Vortex State in Topological Superconductors*, MSc Applied Physics, 2023. [Original PDF](https://drive.google.com/file/d/1d9WheYwEY23cPAlLxxVOUfQDr9jLdqb9/view).
- **Original code:** [`BdG topological.nb`](https://github.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-/blob/d88d12cded1a954d53c2d36d515525c8100ded1c/BdG%20topological.nb) and [`added vortex.nb`](https://github.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-/blob/d88d12cded1a954d53c2d36d515525c8100ded1c/added%20vortex.nb), pinned to original repository commit `d88d12cded1a954d53c2d36d515525c8100ded1c`.

No external derivation, newly chosen model, Python solver or newly computed plot is used in this edition. The notebook text provides navigation and source-status notes. Equations, derivation steps, figure labels and explanatory passages are preserved through full-page thesis images.

## Thesis pages

The PDF has 14 pages before printed page 1. Thus printed page 46 is PDF viewer page 60. The following complete pages are rendered at 120 dpi without cropping, recolouring or redrawing:

| Reading notebook | Printed pages | Material |
| --- | --- | --- |
| 01 | 17–20, 45–46 | Tight-binding derivation (§2.7), Figure 2.1, and topological Hamiltonian Eqs. (4.8)–(4.13). |
| 02 | 30–31, 46–47 | Mean-field/BdG construction, spin matrix and original spectrum Figure 4.20. |
| 03 | 31–33 | Original self-consistency equations, LDOS expression and calculation discussion. |
| 04 | 34–39, 47–48 | Original vortex expression, Eqs. (4.2)–(4.5), conventional vortex figures and Figures 4.20–4.22. |
| 05 | 32, 39–46, 49–52 | Original LDOS/impurity definitions, conventional impurity figures and topological impurity Figures 4.23–4.29. |

Across the notebooks this includes all Chapter 4 figure pages, with the original captions. Some pages appear in more than one notebook so the associated equation is available where needed. [reference_pages/manifest.json](reference_pages/manifest.json) records the PDF hash, printed/PDF page mapping and SHA-256 hash of each of the 27 PNG files.

## Original code cells

The input numbering below counts nonempty `Code`/`Input` calculation cells. It excludes empty cells and the final empty drawing canvas. Source line numbers refer to the text of the pinned `.nb` files.

| Input | BdG topological.nb line | added vortex.nb line | Original content |
| --- | ---: | ---: | --- |
| 1 | 23 | 23 | Parameters, hopping dictionaries, 20 × 20 lattice, site indices, periodic matrices and `Hop1x`. |
| 2 | 922 | 948 | Onsite pairing construction `SCs`. |
| 3 | 1049 | 1075 | Fermi function, coupling, density target and pairing initialization. |
| 4 | 1398 | 1437 | Commented alternative BdG construction. |
| 5 | 1441 | 1480 | Complete active 100-step iteration, observables, commented vortex prototype and exports. |
| 6 | 4032 | 45843 | Original sorting input. |
| 7 | 4073 | 45884 | Original `ListPlot`, including its original plot label, range and grid lines. |

All seven cells from each source are stored in `original_code/*.wl`. The paired JSON files preserve their original `BoxData`, decoded text, original line number and hash. Each notebook code block is an exact contiguous excerpt of that decoded text; metadata records its source cell and character offsets. The extractor flattens `RowBox` string lists and decodes the outer notebook string representation, including escaped line continuations. It does not evaluate or translate Wolfram expressions. Named Wolfram characters and embedded plot-label formatting are retained.

## Preserved source differences and missing workflows

- Original `Nx=20`, `Ny=20`, `iter=100`, coupling and density targets are retained. No small-lattice substitute is used.
- `BdG topological.nb` and `added vortex.nb` keep their separate parameter values and plot labels. Labels in saved plotting cells are not assumed to match all current initialization values.
- The hopping dictionary differs from the printed thesis Hamiltonian. The unused `Hop1x` definition also differs from the active loop. Both originals are preserved.
- The vortex prototype stays inside its original comment delimiters. Its original arctangent and length expression are not replaced or activated.
- `Vimp=0` and the commented `+imp` term are shown as found. No active impurity sweep or active routine producing all thesis LDOS plots was found in these input cells.
- The original `.nb` files contain saved Mathematica graphics. They remain available in the original files; the notebook figure pages come from the thesis PDF, rather than a recreation of those graphics.
- Original equations and interpretations are shown as printed, including unresolved convention issues. Their reproduction is not a claim of independent mathematical or experimental verification.

The source-only revision supersedes the previous tutorial with new Python models, diagnostic examples and generated curves. That earlier version is retained only in Git history.
