# Validation of code, explanations and diagrams

Validation date: 19 October 2026.

## Completed checks

- All five reading notebooks pass `nbformat.validate`.
- The 16 retained original excerpts match their recorded input-cell character ranges.
- All 11 adapted Mathematica sections in notebook 05 match `scripts/vortex_ldos.wl` exactly.
- The uploaded LDOS script is preserved byte for byte (SHA-256 in `SOURCE_MAP.md`).
- Original `.nb`, `.wdx`, decoded source and archived image hashes are unchanged.
- All 30 original diagrams remain represented; displayed crops retain the exact source pixels.
- The notebook text contains no full-page scans or document-comparison narrative.
- Independent NumPy checks pass for an analytically known isolated-site LDOS, Hermiticity, opposite-momentum particle-hole symmetry, eigensystem residuals, complete spectral weights, the integrated-LDOS sum rule and full phase winding.
- Adapted Wolfram sources have balanced brackets, nested comments and strings. This is a lexical check, not a Mathematica compilation.

## Checks available but not run

No Mathematica or Wolfram runtime is installed in the editing environment. `tests/ldos_checks.wlt` is provided for `TestReport` in Wolfram. Import of the actual WDX dataset, full adapted-program execution, convergence of relaxed 20 × 20 cases and numerical reproduction of archived figures are **unverified**. No synthetic or unexecuted result is presented as an original plot.

The program performs matrix dimension/layout and Hermiticity checks at runtime. It refuses unconfirmed impurity baselines and unconverged relaxation through its high-level entry point. These guards do not replace numerical validation of a scientific calculation.

## Reproduce completed checks

```bash
python -m pip install -r requirements.txt
python tools/build_notebooks.py
python tools/verify_sources.py
python tests/check_ldos_algebra.py
```

The Wolfram suite can be run separately with `TestReport["/absolute/path/to/repo/tests/ldos_checks.wlt"]`.
