# BdG code and diagrams

Five reading notebooks explain the original Mathematica code and the existing plots. Each section describes a code block's inputs, operations and outputs. Each diagram has a direct explanation of its axes, colours, curves and visible features.

## Notebook guide

| Notebook | What you will learn |
| --- | --- |
| [01 — Lattice and hopping matrices](notebooks/01_lattice_and_normal_hamiltonian.ipynb) | How site indices, hopping dictionaries, periodic boundaries and spin-mixing matrices are constructed. |
| [02 — BdG matrix and spectra](notebooks/02_bdg_matrix_and_diagonalization.ipynb) | What each matrix block contains, what `Eigensystem` returns, and how `ListPlot` displays energy levels. |
| [03 — Occupations and pairing iteration](notebooks/03_self_consistent_pairing.ipynb) | How occupation weights, density, pairing proposals, mixing and chemical-potential updates work. |
| [04 — Vortex code and spatial diagrams](notebooks/04_vortex_profile_and_spectrum.ipynb) | What the commented vortex block computes, and how to read line profiles, colour maps, surfaces and an energy spectrum. |
| [05 — Impurity code and LDOS diagrams](notebooks/05_impurities_and_ldos.ipynb) | Which impurity terms are active, how to distinguish energy and index axes, and how to read the saved LDOS and spatial plots. |

The notebooks focus on code and diagrams. There are no full-page scans, document comparisons or chapter/equation-number references in the explanations. All 30 existing diagrams are retained as figure-only regions. Original plotted values, axes, labels and colours are unchanged.

## Code and execution

The `.ipynb` files are reading guides with Mathematica code blocks. To execute the calculation, open [BdG topological.nb](BdG%20topological.nb) or [added vortex.nb](added%20vortex.nb) in Mathematica or Wolfram Desktop and evaluate the input cells in order.

The original code uses a 20 × 20 lattice and a fixed 100-step iteration. Its machine-specific `SetDirectory[...]` path needs adjustment on another machine. Its export commands can overwrite the `.wdx` files. The saved vortex block and impurity contribution are commented out; the guides identify that execution status directly.

The original `.nb` and `.wdx` files are unchanged. Complete readable input-cell transcriptions are in [original_code](original_code/). The figures are existing archived outputs; this editing pass does not claim a fresh numerical reproduction.

## Maintain the guides

```bash
python -m pip install -r requirements.txt
python tools/extract_original_code.py
python tools/extract_diagrams.py
python tools/build_notebooks.py
python tools/verify_sources.py
```

The verification checks that the displayed code excerpts match the original input cells and that every diagram retains exactly the pixels of its recorded source region. Internal asset provenance is recorded in [SOURCE_MAP.md](SOURCE_MAP.md); the check results are described in [VALIDATION.md](VALIDATION.md).
