# LDOS, impurities and vortex-core views

The unchanged user-supplied input is `original_code/vortex_ldos_original.m`. The revised program is `scripts/vortex_ldos.wl`, with individually evaluated examples in `scripts/ldos_examples.wl`. Notebook 05 explains each block directly and retains the archived figures.

## What changed from the supplied script

| Supplied expression or behaviour | Adaptation and reason |
|---|---|
| Absolute local `SetDirectory` | Explicit input folder; no global working-directory change. |
| Lattice dimensions used before assignment | Explicit dimensions checked against every imported matrix. |
| Diagonal portions of `SCs` used for pairing; `SCs2c` reuses `SCs1` | Extract the off-diagonal onsite pairing block and verify the complete saved layout. |
| `Conjugate[Huuc]` at the same Bloch momentum | Hole block uses the opposite electron momentum before conjugation in the stored hole-sign basis. |
| `Conjugate[huv]` | Explicit conjugate transpose for the lower-left partner. |
| Impurity configuration only present through imported matrices | Explicit scalar or spin-diagonal magnetic potential, starting from a user-verified clean baseline. |
| Full eigensystem sum of one electron and one hole component | Full eigensystem electron spectral weight averaged over the two spins. This has an explicitly different normalization; historical curve heights are not promised. |
| Interleaved site values exported as two-column rows | Energy/site/LDOS CSV; complete metadata and arrays saved as WXF. |
| Final `ListPlot` uses an implicit sample index and omits `Ns^2` | Explicit energy coordinates and the same momentum averaging in all plots and export. |
| Filename uses unassigned Greek `Eta` rather than `eta` | Metadata stores actual numerical broadening; output folder is explicit. |
| Unconditional `DeleteFile`, kernel shutdown and relaunch | No deletion or parallel-kernel side effects; serial momentum loop. |
| No executable impurity/vortex-map plotting routines | Clearly labelled new wrappers for spectra, pairing cuts, maps, surfaces and optional spatial LDOS. |

## Pairing and basis conventions

The onsite pairing block is `B={{0,D},{D,0}}`. With `S=diag(I,-I)`, the full matrix is `{{h(k),B},{B†,-S Conjugate[h(-k)] S}}`. This retains the saved zero-momentum signs while giving a consistent opposite-momentum hole sector. The first two groups of components are electrons. The spin-averaged electron LDOS integrates to one over the whole energy axis for an orthonormal complete eigensystem. The original electron-plus-hole expression integrates to two; this is a normalization/convention change, not evidence by itself that the historical observable was wrong.

Energy values, broadening, coupling and temperature `kT` share model energy units. Old figures mix eV and meV labels; this workflow does not invent a conversion. Coordinates and coherence length use lattice spacings. It imposes a vortex envelope on a periodically repeated cell and does not add a magnetic gauge field or establish magnetic flux boundary conditions.

## Fixed field versus updated field

`PairingMode -> "Imported"` preserves the input pairing. Adding a potential may change spectra and LDOS but cannot change a map of that fixed pairing. `SeedVortex -> True` explicitly replaces the field with a reconstructed envelope; otherwise the imported field is untouched.

`PairingMode -> "Relax"` requires an explicit vortex seed and adapts the onsite anomalous sum, clipped Fermi function, 50/50 mixing and common onsite feedback from `added vortex.nb`. It imposes the vortex factor at every iteration, as the saved code does. This is a constrained update, not a validated new fully unconstrained vortex solver. The `Winding` phase option uses quadrant-aware `ArcTan[dx,dy]`; `OriginalAbsoluteAngle` keeps the original absolute-angle expression for comparisons. The choice is recorded. No convergence or historical agreement is assumed.

To diagnose nonconvergence, call `RelaxPairing` directly and inspect its returned `History` and `Converged` fields. The high-level runner stops when convergence fails. The final matrix is rebuilt after the last updates, avoiding mixed iteration ages for the field and normal blocks.

## Saved-figure settings and missing evidence

The preserved document's figures 4.9 and 4.12–4.19 identify 20 × 20 with 20 supercells for the vortex/impurity LDOS, while the supplied script defaults to 10. Figures 4.12–4.14 use scalar strengths 0.5, 50, 5; figures 4.17–4.19 use magnetic strengths 0.5, 5, 50. The pairing-cut panels (4.10) use 0.5, 1, 5. Identified observation sites are 190 (core) and 180 (historically labelled nearby); the far site's exact index is unspecified. The example uses site 1 and labels it explicitly as a new choice.

With the source loop ordering, site 190=(10,10), 180=(9,20), 1=(1,1). Site 180 is not an immediate neighbour of the core. Plot axes saying "Number of Sites" conflict with the supplied final plotting command's energy-row iteration at a fixed site. The adaptation does not relabel archived images or assert their exact provenance; it supplies energy coordinates for new curves.

Earlier no-vortex LDOS figures use 28 × 28 and 40 supercells and different hopping settings. They require separately prepared compatible matrices. The bundled WDX files are a single unlabelled dataset, not a complete set of inputs for every figure. No full historical reproduction is claimed.

The archived vortex-core heatmaps and surfaces are explained as pairing-field views. Their exact stored field arrays and phase/magnitude/sign plotting choices are not supplied for every case. The new views deliberately show `Abs[Delta]` with explicit coordinates, a colour legend and conventional increasing height; they are reconstructed plotting functions, not recovered historical commands.

## Running and checking

1. Load `scripts/vortex_ldos.wl` in Mathematica/Wolfram Desktop using `Get`.
2. Set the input folder and lattice dimensions. For impurity additions, verify a clean baseline first.
3. Choose imported pairing or explicit seed/relaxation, sites, momentum samples and energy grid.
4. Evaluate a case and check `FailureQ[result]` before plotting.
5. Export to a new output folder. The original matrix files are never overwritten by this workflow.

`tests/ldos_checks.wlt` is a Wolfram `TestReport` suite for a small synthetic model. It checks the matrix symmetry, zero-impurity limit, spin signs, vortex envelope, sum rule, LDOS dimensions and positivity. It has not been executed in this editing environment because no Wolfram runtime is installed.

`python tests/check_ldos_algebra.py` independently checks the documented matrix and spectral conventions with NumPy, including an isolated-site analytic LDOS limit. It also checks balanced Wolfram delimiters/comments/strings, which is not a Wolfram parser or runtime test. Full Wolfram execution and the historical numerical curves remain unverified.
