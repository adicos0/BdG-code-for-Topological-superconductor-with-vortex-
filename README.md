# BdG calculations for a topological superconductor with a vortex

Original Mathematica code and a guided reading of Aditya Banerjee's 2023 MSc Applied Physics thesis, [*Understanding the Magnetic Vortex State in Topological Superconductors*](https://drive.google.com/file/d/1d9WheYwEY23cPAlLxxVOUfQDr9jLdqb9/view).

## Read the five notebooks

These notebooks use **the original thesis pages, figures and derivations, together with the original Mathematica code**. They preserve the printed notation, captions, original parameters and commented sections. The displayed plots are archived thesis figures. They are not newly simulated examples.

| Notebook | Original material |
| --- | --- |
| [01 — Lattice and normal-state Hamiltonian](notebooks/01_lattice_and_normal_hamiltonian.ipynb) | Thesis §2.7 and §4.4; original hopping dictionaries, site indexing and periodic matrix construction. |
| [02 — BdG matrix and diagonalization](notebooks/02_bdg_matrix_and_diagonalization.ipynb) | Thesis BdG equations; original `SCs`, matrix assembly, `Eigensystem` and both `ListPlot` commands; Figure 4.20. |
| [03 — Self-consistency](notebooks/03_self_consistent_pairing.ipynb) | Original thesis mean-field equations; complete parameter and iteration cells from both Mathematica files. |
| [04 — Vortex profile and spectrum](notebooks/04_vortex_profile_and_spectrum.ipynb) | Thesis §4.2 and vortex figures; original commented vortex code; topological vortex Figures 4.20–4.22. |
| [05 — Impurities and LDOS](notebooks/05_impurities_and_ldos.ipynb) | Original LDOS and impurity equations; the source's disabled impurity terms; original conventional and topological impurity plots. |

The `.ipynb` files are **reading notebooks with Mathematica code blocks**, rather than executable Python translations. The original equations and derivations are shown as full-page images to preserve every sign, index and intermediate step. Open a page's full-size link to zoom. The images are also stored locally in [reference_pages](reference_pages/).

## Original files and execution

| File | Contents |
| --- | --- |
| [BdG topological.nb](BdG%20topological.nb) | Original Mathematica lattice, BdG, iteration and spectrum notebook, including saved outputs. |
| [added vortex.nb](added%20vortex.nb) | Related original notebook, including its commented vortex prototype and saved outputs. |
| `Huu.wdx`, `Hdd.wdx`, `Hud.wdx`, `Hdu.wdx`, `SCs.wdx` | Original exported Wolfram matrix data. |
| [original_code/bdg_topological.wl](original_code/bdg_topological.wl) | Readable transcription of all seven nonempty calculation input cells in `BdG topological.nb`. |
| [original_code/added_vortex.wl](original_code/added_vortex.wl) | Readable transcription of the corresponding seven input cells in `added vortex.nb`. |
| `original_code/*.json` | The original `BoxData` expressions, source line numbers and readable transcriptions for comparison. |

To execute the original calculation, open a `.nb` file in Mathematica or Wolfram Desktop and evaluate its input cells in order. The source was created with Mathematica 12.0 and contains a machine-specific `SetDirectory[...]`; adjust that path in your working copy. Original `Export[...]` commands can overwrite the `.wdx` files. The 20 × 20 lattice gives a 1600 × 1600 BdG matrix.

The original seven `.nb` and `.wdx` files are unchanged. The `.wl` files are readable transcriptions of notebook box expressions; the `.nb` files remain authoritative for native evaluation and saved graphics.

## What is and is not contained in the saved source

The original files contain the active lattice/BdG calculation and a fixed 100-step iteration. The vortex prototype is commented out, and the impurity contribution is disabled in the active matrix assembly. The saved inputs do not provide an active impurity sweep or an active LDOS plotting routine that reproduces every thesis figure. Those figures are therefore shown directly from the thesis, with their original captions.

The printed Hamiltonian and the active hopping dictionaries differ. The notebooks preserve both sources rather than replacing one with a new implementation. Original claims and notation in the page images are attributed to the thesis; reproducing a page is not an independent validation of its physical interpretation.

See [SOURCE_MAP.md](SOURCE_MAP.md) for exact source locations and [VALIDATION.md](VALIDATION.md) for the provenance checks.

## Rebuild or verify the reading notebooks

```bash
python -m pip install -r requirements.txt
python tools/extract_original_code.py
python tools/build_notebooks.py
python tools/verify_sources.py
```

These commands organize and verify the supplied material. They do not run a BdG solver or generate new numerical curves. The previous Python tutorial implementation and its added plots have been removed from the current edition; the earlier commit remains in Git history.
