# Sources, conventions and translation boundaries

These notebooks accompany Aditya Banerjee's 2023 MSc Applied Physics project, [*Understanding the Magnetic Vortex State in Topological Superconductors*](https://drive.google.com/file/d/1d9WheYwEY23cPAlLxxVOUfQDr9jLdqb9/view). Page references below use the printed thesis page numbers, not PDF viewer indices. The 2023 BdG thesis is distinct from the later anyon/Rydberg thesis.

## What comes from the committed code

Source snapshot: [`d88d12cded1a954d53c2d36d515525c8100ded1c`](https://github.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-/tree/d88d12cded1a954d53c2d36d515525c8100ded1c). Original `.nb` and `.wdx` files are preserved. No Python implementation existed in that snapshot. The new Python functions are translations or explicitly labelled extensions, not original Python supplied by the author. Mathematica was not available for an independent rerun of saved output.

The references to cells below count the initial `Code` cell and the subsequent `Input` cells. Empty or saved-output cells are excluded. Line numbers refer to the original text representation of each `.nb` file, at the pinned commit.

| New function or section | Mathematica source | Thesis location | Status |
| --- | --- | --- | --- |
| `coordinates`, `legacy_blocks` | `BdG topological.nb`, cell 1, line 23; `added vortex.nb`, cell 1, line 23 | Related lattice discussion, §4.4 | Translation of original periodic indexing and hopping dictionaries; tutorial sizes reduced. |
| `legacy_bdg` | `SCs` in cell 2, lines 922 / 948; active `hh1`, `hh2`, `H` in cell 5, lines 1441 / 1480 | Eqs. (4.10)–(4.12), pp. 45–46 | Translation of active matrix; canonical sign transformation is added algebra. |
| `source_fermi` | Cell 3, lines 1049 / 1075 | Mean fields, §4.1, pp. 30–33 | Retains source cutoff and thermal energy. Stable exponential evaluation is an implementation change. |
| `source_observables`, `source_iteration` | Cell 5, lines 1441 / 1480 | Mean-field framework, §4.1 | Retains full-state sums, mixing, chemical-potential update and 100 steps. Added residual diagnostics and a final evaluation do not change updates. |
| `thesis_normal` | Contrasted with cell 1; **not** a translation of the active hopping | Eqs. (4.8)–(4.9), pp. 45–46 | New implementation of the stated thesis lattice model. |
| `canonical_bdg`, `canonical_checks` | Compared algebraically with active cell 5 | Pairing/BdG discussion, §§4.1, 4.4 | Explicit canonical basis and verification helpers. |
| `vortex_field` | Commented vortex prototype within cell 5 | Eq. (4.2), p. 34; phase discussion Eq. (4.4) | Extension with corrected `atan2` phase and explicit open boundaries. |
| `add_impurity` | `Vimp=0` in cell 1; `+imp` disabled in cell 5 | Eqs. (4.6)–(4.7), pp. 39–40; §§4.4.1–4.4.2 | New implementation of thesis impurity terms; frozen gap. |
| `electron_ldos` | No active LDOS routine translated | LDOS discussion, §4.1, p. 32 | New full-spectrum, spin-summed convention; derivation and sum rule shown. |

## Differences that must remain visible

1. **Hopping range and coefficients.** The active code uses diagonal neighbours with amplitude `m/4` and nearest axial hopping `-m`. The thesis mass uses axial distance-two hopping `m/8` and nearest axial hopping `-m/2`, with a spin `sigma_z` factor. They are different models.
2. **Spin mixing.** The active code assigns `lambda*(1-i)` to all four axial directions. The thesis model requires directional imaginary hoppings to generate sine terms in momentum space.
3. **Spin-diagonal blocks.** The unused `Hop1x` definition has opposite `H1` signs. The active loop uses the same `H1` in `Huu` and `Hdd`. The translation follows the loop.
4. **Chemical-potential sign.** The source has `H1 - Chempot1[mu]` with `Chempot1=-mu I`, so its active block contains `+mu I`. Thesis Eq. (4.8) contains `-mu I`. Each implementation retains its stated convention; the same numeric `mu` must not be assumed physically equivalent across them.
5. **Nambu signs.** The symmetric source pairing block is mapped to antisymmetric canonical spin-singlet pairing by a down-hole sign change and `Delta=-st`. This exact matrix identity is demonstrated and tested, including complex gaps. It is an algebraic interpretation of the saved matrix, not a documented basis statement from its author.
6. **Vortex status and phase.** The vortex update is commented out in both originals. The absolute-value arctangent in that prototype and thesis Eq. (4.4) cannot wind through all quadrants. The extension uses `atan2`, with its winding explicitly checked.
7. **Length units.** The relation in thesis Eq. (4.3) cannot be used as an SI coherence length without further unit conventions. The tutorial specifies `xi/a` independently.
8. **Time reversal.** The sentence after thesis Eq. (4.8) describes the mass as time-reversal preserving. For physical spin time reversal, an even-in-momentum `M_k sigma_z` changes sign. The notebooks show this correction and cite the primary model paper below.
9. **Boundary and field assumptions.** Source translation uses periodic boundaries. Isolated-vortex extensions use finite open boundaries and imposed pairing, without vector potential or screening. They are not a periodic magnetic-vortex-lattice solution.
10. **Convergence and identification.** One hundred updates are not a convergence criterion. A low-energy eigenvalue or broadened LDOS peak alone does not identify a Majorana mode. The notebooks report residuals and limit their claims accordingly.

## Tutorial parameters versus archived calculations

| Calculation | Executed tutorial | Original / thesis context |
| --- | --- | --- |
| Normal-state comparison | 6 × 6 periodic; m=0.5, lambda=1, mu=0 | Code uses 20 × 20. |
| Matrix/basis demonstration | 4 × 4 periodic; mu=0.13, st=-0.05 | Nonzero example chosen to test signs. |
| Source iterations | 4 × 4 periodic; 100 steps; initial st=-0.01, mu=0 | Original 20 × 20; couplings/targets retained: (-1.45, 1.2), (-1.25, 1.0). |
| Vortex and impurity extensions | 11 × 11 open; m=0.5, lambda=1, mu=-0.7, Delta0=0.2, xi/a=1.5 | The thesis discusses 20 × 20 samples; these finite-size examples are not figure reproductions. |
| Impurity/LDOS sweep | V/E0=0, 0.5, 2; eta/E0=0.02; fixed pairing | Demonstration settings; no self-consistent impurity relaxation. |

All parameters expressed as energies are in a common chosen energy unit `E0`. No material-specific SI calibration is inferred. Computational functions use dense matrices, so memory grows as N² and full diagonalization roughly as N³ for N sites. Increasing from 4 × 4 to 20 × 20 greatly increases the iterative cost.

## Verified external references

- Liang Fu and C. L. Kane, *Superconducting Proximity Effect and Majorana Fermions at the Surface of a Topological Insulator*, Physical Review Letters **100**, 096407 (2008). [arXiv:0707.1692](https://arxiv.org/abs/0707.1692); [publisher DOI: 10.1103/PhysRevLett.100.096407](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.096407). Background for the induced-superconductivity/vortex setting; not evidence that the saved code implements every condition of that theory.
- Vedangi Pathak, Stephan Plugge and Marcel Franz, *Majorana bound states in vortex lattices on iron-based superconductors*, Annals of Physics **435**, 168431 (2021). [arXiv:2012.10588](https://arxiv.org/abs/2012.10588); [full text, §§2.1–2.2](https://arxiv.org/html/2012.10588v2); [DOI: 10.1016/j.aop.2021.168431](https://doi.org/10.1016/j.aop.2021.168431). Used to verify the lattice mass, basis conventions and the need for additional gauge treatment under periodic vortex boundary conditions.

No new result in these tutorials is claimed as a reproduction of a published figure. Derivations and numerical diagnostics added for explanation are labelled explicitly in each notebook.
