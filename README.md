# BdG calculations for a topological superconductor with a vortex

This repository contains Mathematica notebooks from my MSc Applied Physics project, [*Understanding the Magnetic Vortex State in Topological Superconductors*](https://drive.google.com/file/d/1d9WheYwEY23cPAlLxxVOUfQDr9jLdqb9/view) (Aditya Banerjee, 2023). The project asks what happens to the low-energy electronic states when superconductivity is induced on the surface of a topological insulator and a vortex is introduced into the superconducting order parameter.

## Readable Python notebooks

Start with [01 — Lattice and normal-state Hamiltonian](notebooks/01_lattice_and_normal_hamiltonian.ipynb). Each notebook is self-contained, with short explanations, equations, defined symbols, visible code sections and saved numerical output. The Python is translated from the original Mathematica calculations where labelled; vortex and impurity examples implement explicitly identified thesis-based extensions.

| Notebook | What it explains |
| --- | --- |
| [01 — Lattice and normal-state Hamiltonian](notebooks/01_lattice_and_normal_hamiltonian.ipynb) | Site indexing, original hopping, thesis Eqs. (4.8)–(4.9), and their differences. |
| [02 — BdG matrix and diagonalization](notebooks/02_bdg_matrix_and_diagonalization.ipynb) | Nambu basis, pairing signs, eigenvectors, particle–hole symmetry and an atomic-limit check. |
| [03 — Self-consistent pairing](notebooks/03_self_consistent_pairing.ipynb) | Original density/gap sums, mixing, chemical-potential update and convergence diagnostics. |
| [04 — Vortex profile and spectrum](notebooks/04_vortex_profile_and_spectrum.ipynb) | Thesis vortex ansatz, corrected phase winding, open boundaries and spatial probability. |
| [05 — Impurities and LDOS](notebooks/05_impurities_and_ldos.ipynb) | Scalar/magnetic perturbations, frozen-profile comparison and LDOS normalization. |

**Read the source distinction:** the active Mathematica hopping is different from the Hamiltonian printed in the thesis. Notebook 03 translates the active code; notebooks 04–05 use the separately labelled thesis model. Tutorial sizes are smaller than the archived 20 × 20 calculation. They do not reproduce the thesis figures. See [SOURCE_MAP.md](SOURCE_MAP.md) for exact provenance, corrections and assumptions, and [VALIDATION.md](VALIDATION.md) for execution and numerical checks.

To run locally from this repository:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
jupyter lab notebooks
```

Run cells from top to bottom. No Wolfram installation or `.wdx` import is needed for the Python notebooks. To run the numerical checks and execute all notebook cells:

```bash
python -m pytest -q tests
python tools/execute_notebooks.py
```

Maintainers can edit the tested functions in `python/bdg_models.py` and prose/examples in `tools/build_notebooks.py`, then run the builder followed by the executor. The builder embeds the function definitions in each notebook so readers can inspect all calculation code directly. Rebuilding clears saved outputs; execution regenerates them.

## The idea

The calculation starts with electrons on a two-dimensional square lattice. The topological-insulator surface model has spin-dependent hopping, a lattice mass term, and a chemical potential. Adding onsite *s*-wave pairing gives a Bogoliubov–de Gennes (BdG) Hamiltonian, which treats electron and hole amplitudes together. Diagonalizing this Hamiltonian gives the excitation energies and wavefunctions. The notebooks also contain an iterative calculation of the pairing field and an adjustment of the chemical potential toward a chosen particle density.

In the thesis, a vortex is modeled by a pairing amplitude that goes to zero at its center and recovers away from it, together with a position-dependent phase. The spectra are then examined for states inside the superconducting gap. The thesis compares conventional and topological superconducting cases and discusses how magnetic and non-magnetic impurities affect vortex-core states.

## What is in this repository

| File | Contents |
| --- | --- |
| [`BdG topological.nb`](BdG%20topological.nb) | Mathematica notebook constructing the lattice and BdG matrices, iterating the pairing calculation, and plotting an energy spectrum. |
| [`added vortex.nb`](added%20vortex.nb) | A related research notebook with a prototype vortex-profile expression and saved numerical output and plots. |
| `Huu.wdx`, `Hdd.wdx`, `Hud.wdx`, `Hdu.wdx` | Wolfram data files for spin-resolved blocks exported by the notebooks. |
| `SCs.wdx` | Wolfram data file for an exported superconducting pairing matrix. |

The `.nb` files were created with Mathematica 12.0. To explore them, open a notebook in Mathematica or Wolfram Desktop and evaluate its input cells in order. The first input cell contains a machine-specific `SetDirectory[...]` path; change it to the directory where you downloaded the repository. Evaluating the calculation can overwrite the `.wdx` files. A 20 × 20 lattice with spin and particle–hole degrees of freedom produces a 1600 × 1600 BdG matrix, so full diagonalization and iteration may take time.

## Scope and interpretation

**Verified from the committed notebooks:** the saved code builds a spinful lattice BdG matrix, calls `Eigensystem`, updates the pairing field iteratively, and exports matrix blocks. A vortex profile is written in the notebooks, but that section is commented out in the committed input cells; the impurity term is also disabled in the main calculation. These files are research snapshots, not a one-command reproduction of every thesis figure.

**Physical interpretation discussed in the thesis:** in-gap states near a vortex can be candidates for vortex-core bound states. A state near zero energy in a finite-system spectrum can be consistent with a Majorana mode, but energy alone does not establish one. Its spatial localization, particle–hole structure, stability with system size and parameters, and the topology and boundary conditions of the model need separate checks. The archived notebook output has not been independently rerun or validated here.

For the full model, calculation choices, and figures, see [the MSc project PDF](https://drive.google.com/file/d/1d9WheYwEY23cPAlLxxVOUfQDr9jLdqb9/view), especially Chapter 4.
