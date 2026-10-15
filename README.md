# BdG calculations for a topological superconductor with a vortex

This repository contains Mathematica notebooks from my MSc Applied Physics project, [*Understanding the Magnetic Vortex State in Topological Superconductors*](https://drive.google.com/file/d/1d9WheYwEY23cPAlLxxVOUfQDr9jLdqb9/view) (Aditya Banerjee, 2023). The project asks what happens to the low-energy electronic states when superconductivity is induced on the surface of a topological insulator and a vortex is introduced into the superconducting order parameter.

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
