# Validation of the Python notebooks

Validation date: 19 October 2026. These checks validate the Python translation and tutorial calculations. They do not establish reproduction of the archived Wolfram output or every thesis figure.

## Execution

All five notebooks were executed from top to bottom in separate Python kernels. The committed notebooks include numerical outputs and eight figures. Notebook JSON was checked with `nbformat.validate`; executed cells contained no error outputs. Embedded function definitions were checked against `python/bdg_models.py` using Python's abstract syntax trees.

Commands:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q tests
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python tools/execute_notebooks.py
```

Environment used: Python 3.12; NumPy 2.5.4; SciPy 1.18.1; Matplotlib 3.11.2; nbformat 5.11.1; nbclient 0.11.0; ipykernel 7.4.0; pytest 9.1.1. The version ranges in `requirements.txt` allow compatible environments; these recorded versions identify the actual run. Small floating-point differences are expected across platforms.

The restricted macOS execution environment emitted a process-enumeration permission diagnostic during kernel shutdown. It occurred outside the notebook cells; the executor completed successfully and saved the outputs. No Wolfram kernel was run, and no `.wdx` file was loaded or overwritten.

## Numerical checks — seven passed

1. **Source hopping:** compare the translated matrices against an independent implementation of the source's ordered nine-image periodic lookup, including the transpose.
2. **Thesis Fourier transform:** project the real-space Hamiltonian onto every allowed plane wave of a 5 × 6 periodic lattice. Check the complete 2 × 2 spin blocks against Eq. (4.8), including imaginary spin-mixing signs, and compare the full energy bands.
3. **Basis transformation:** use a spatially varying complex gap to verify the exact signed-hole map between source and canonical BdG matrices. Check Hermiticity, particle–hole symmetry and paired eigenvalues.
4. **Limits:** compare an isolated paired site against ±sqrt(mu²+|Delta|²). At zero pairing, compare the full BdG spectrum with electron energies and their negatives.
5. **Occupation convention:** compare the source down-hole density formula with the electron-sector density from the same full eigensystem. Check physical occupation bounds and f(-E)=1-f(E).
6. **Vortex:** verify positive and negative unit winding on closed sampled loops, zero core amplitude, and reversal by complex conjugation.
7. **Impurity and LDOS:** check the zero-potential limit, scalar/magnetic spin signs, BdG symmetry, nonnegative LDOS and total spin-summed electron spectral weight two per site.

These checks use analytic limits, independently constructed expressions and conservation/normalization constraints. They do not constitute experimental validation or a topological invariant calculation.

## Selected recorded results

| Quantity | Result |
| --- | --- |
| Source matrix versus basis-transformed canonical matrix | Maximum entry difference 0 |
| Relative eigen-equation residual in notebook 02 | 4.114 × 10⁻¹⁵ |
| Eigenvector orthonormality residual | 1.366 × 10⁻¹⁴ |
| Maximum ±energy mismatch | 2.398 × 10⁻¹⁴ E0 |
| Measured vortex winding | +1 |
| Imposed amplitude at the core | 0 |
| Corner amplitude divided by Delta0 | 0.999839 |
| Four-state spatial probability sum | 0.9999999999999998 |
| LDOS integrated spectral weight per site, checked by completeness | 2 |

**Convergence limitation:** neither original 100-step update reaches both 10⁻⁶ tolerances on the executed 4 × 4 lattice. For `BdG topological.nb` parameters, the final gap residual is approximately 9.333 × 10⁻⁵ E0 and the density residual is 7.505 × 10⁻⁴ electrons/site. For `added vortex.nb` parameters, they are 1.141 × 10⁻⁵ E0 and 4.791 × 10⁻⁴ electrons/site. The notebooks report these as not converged at those tolerances.

**Spectral limitation:** the four states nearest zero in the open-sample vortex example have substantial boundary weight. Their small energies are not presented as proof of vortex-core Majorana modes. The scalar and magnetic sweeps share the displayed minimum-absolute-energy values for this example while their core LDOS differs; the notebooks do not infer equivalence or general impurity protection.

## Preserved source files

SHA-256 of the original Mathematica notebooks:

```text
98a24f8f660e64ecab336b95d4f63d82347007a071eb441a7bb52a46f9361389  BdG topological.nb
91ec9c8fe23531115a81ad58c16104109e98fef8f930a7fe3a77a2c36939e9cb  added vortex.nb
```

The original seven `.nb` / `.wdx` files were not edited. Provenance, page references, model differences and explicit extensions are documented in [SOURCE_MAP.md](SOURCE_MAP.md).
