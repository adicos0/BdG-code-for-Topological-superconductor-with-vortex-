# BdG code and diagrams

Five reading notebooks explain the original Mathematica code and the existing plots. Each section describes a code block's inputs, operations and outputs. Each diagram has a direct explanation of its axes, colours, curves and visible features.

## Notebook guide

| Notebook | What you will learn |
| --- | --- |
| [01 — Lattice and hopping matrices](notebooks/01_lattice_and_normal_hamiltonian.ipynb) | How site indices, hopping dictionaries, periodic boundaries and spin-mixing matrices are constructed. |
| [02 — BdG matrix and spectra](notebooks/02_bdg_matrix_and_diagonalization.ipynb) | What each matrix block contains, what `Eigensystem` returns, and how `ListPlot` displays energy levels. |
| [03 — Occupations and pairing iteration](notebooks/03_self_consistent_pairing.ipynb) | How occupation weights, density, pairing proposals, mixing and chemical-potential updates work. |
| [04 — Vortex code and spatial diagrams](notebooks/04_vortex_profile_and_spectrum.ipynb) | What the commented vortex block computes, and how to read line profiles, colour maps, surfaces and an energy spectrum. |
| [05 — Impurity code and LDOS diagrams](notebooks/05_impurities_and_ldos.ipynb) | Exact excerpts from the supplied LDOS script, original impurity/vortex blocks, and direct explanations of the saved diagrams. |

The notebooks focus on code and diagrams. There are no full-page scans, document comparisons or chapter/equation-number references in the explanations. All 30 existing diagrams are retained as figure-only regions. Original plotted values, axes, labels and colours are unchanged.

## Code and execution

The `.ipynb` files are reading guides with Mathematica code blocks. To execute the calculation, open [BdG topological.nb](BdG%20topological.nb) or [added vortex.nb](added%20vortex.nb) in Mathematica or Wolfram Desktop and evaluate the input cells in order.

The original `BdG topological.nb` code uses a 20 × 20 lattice and a fixed 100-step iteration. Its machine-specific `SetDirectory[...]` path needs adjustment on another machine. Its export commands can overwrite the `.wdx` files. The topological file comments out its vortex and impurity contributions; the `added vortex.nb` transcription includes an active vortex update. Notebook 05 explains those original blocks without replacing their expressions.

The original `.nb` and `.wdx` files are unchanged. Complete readable input-cell transcriptions are in [original_code](original_code/). The figures are existing archived outputs; this editing pass does not claim a fresh numerical reproduction.

## Original LDOS code walkthrough

The uploaded [Vortex Ldos copy.m](original_code/vortex_ldos_original.m) is preserved byte for byte. Notebook 05 explains its actual parameters, boundary arrays, imported matrices, matrix assembly, eigensystem, spectral sums, export and plotting commands. Every displayed LDOS block is an exact excerpt with source line numbers.

The notebook also shows the original commented impurity terms, the vortex update from `added vortex.nb`, and an original spectrum plotting command. It does not add impurity sweeps, alternative LDOS normalization, a replacement matrix convention, or reconstructed vortex plotting functions. Where the exact plotting command for an archived figure is not supplied, the notebook states that limitation.

See [LDOS source notes](LDOS_WORKFLOW.md). The previously added reconstructed implementation and its examples/tests have been removed from the current tree; the original scientific code and data are retained.

## Maintain the guides

```bash
python -m pip install -r requirements.txt
python tools/extract_original_code.py
python tools/extract_diagrams.py
python tools/build_notebooks.py
python tools/verify_sources.py
```

The verification checks that the displayed code excerpts match the original input cells and that every diagram retains exactly the pixels of its recorded source region. Internal asset provenance is recorded in [SOURCE_MAP.md](SOURCE_MAP.md); the check results are described in [VALIDATION.md](VALIDATION.md).
