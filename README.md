# BdG code and diagrams

Five reading notebooks explain the original Mathematica code and the existing plots. Each section describes a code block's inputs, operations and outputs. Each diagram has a direct explanation of its axes, colours, curves and visible features.

## Notebook guide

| Notebook | What you will learn |
| --- | --- |
| [01 — Lattice and hopping matrices](notebooks/01_lattice_and_normal_hamiltonian.ipynb) | How site indices, hopping dictionaries, periodic boundaries and spin-mixing matrices are constructed. |
| [02 — BdG matrix and spectra](notebooks/02_bdg_matrix_and_diagonalization.ipynb) | What each matrix block contains, what `Eigensystem` returns, and how `ListPlot` displays energy levels. |
| [03 — Occupations and pairing iteration](notebooks/03_self_consistent_pairing.ipynb) | How occupation weights, density, pairing proposals, mixing and chemical-potential updates work. |
| [04 — Vortex code and spatial diagrams](notebooks/04_vortex_profile_and_spectrum.ipynb) | What the commented vortex block computes, and how to read line profiles, colour maps, surfaces and an energy spectrum. |
| [05 — Impurity code and LDOS diagrams](notebooks/05_impurities_and_ldos.ipynb) | The supplied LDOS script, a corrected runnable adaptation, impurity cases, optional pairing updates, spectra and vortex-core plots. |

The notebooks focus on code and diagrams. There are no full-page scans, document comparisons or chapter/equation-number references in the explanations. All 30 existing diagrams are retained as figure-only regions. Original plotted values, axes, labels and colours are unchanged.

## Code and execution

The `.ipynb` files are reading guides with Mathematica code blocks. To execute the calculation, open [BdG topological.nb](BdG%20topological.nb) or [added vortex.nb](added%20vortex.nb) in Mathematica or Wolfram Desktop and evaluate the input cells in order.

The original `BdG topological.nb` code uses a 20 × 20 lattice and a fixed 100-step iteration. Its machine-specific `SetDirectory[...]` path needs adjustment on another machine. Its export commands can overwrite the `.wdx` files. The topological file comments out its vortex and impurity contributions; the `added vortex.nb` transcription includes an active vortex update. Notebook 05 now also provides explicitly adapted functions for these tasks.

The original `.nb` and `.wdx` files are unchanged. Complete readable input-cell transcriptions are in [original_code](original_code/). The figures are existing archived outputs; this editing pass does not claim a fresh numerical reproduction.

## LDOS workflow added from the supplied file

- [Unchanged uploaded Mathematica script](original_code/vortex_ldos_original.m).
- [Adapted executable functions](scripts/vortex_ldos.wl): import checks, impurity potentials, boundary phases, consistent BdG blocks, LDOS, optional pairing relaxation, plots and export.
- [Examples](scripts/ldos_examples.wl), evaluated individually in Mathematica.
- [Workflow and source mapping](LDOS_WORKFLOW.md), including known differences from the archived plots.

Notebook 05 contains the same adapted code in labelled sections. The `.ipynb` remains a Markdown reading guide; the `.wl` is the runnable artifact. No heavy calculation starts when the files are loaded. The program does not infer that the bundled WDX files are a clean impurity-free baseline. Inspect their provenance and provide compatible input data before running a new impurity case.

New LDOS is explicitly spin averaged, uses energy coordinates, and exports site identifiers. Fixed imported pairing cannot produce impurity-dependent changes to a pairing map. Optional relaxation follows the saved constrained-vortex update, records residuals and requires convergence. This is not a verified numerical reconstruction of every historical plot.

## Maintain the guides

```bash
python -m pip install -r requirements.txt
python tools/extract_original_code.py
python tools/extract_diagrams.py
python tools/build_notebooks.py
python tools/verify_sources.py
python tests/check_ldos_algebra.py
```

The verification checks that the displayed code excerpts match the original input cells and that every diagram retains exactly the pixels of its recorded source region. Internal asset provenance is recorded in [SOURCE_MAP.md](SOURCE_MAP.md); the check results are described in [VALIDATION.md](VALIDATION.md).
