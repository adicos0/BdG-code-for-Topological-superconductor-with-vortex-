# Source-explanation validation

Validation date: 19 October 2026.

- All five notebooks pass `nbformat.validate`.
- All displayed Mathematica blocks have original-source provenance; blocks without it are rejected by the verifier.
- The 20 excerpts from the original `.nb` input cells match their exact recorded character ranges.
- The 11 LDOS excerpts match the unchanged user-supplied `.m` file, including original comments and expressions.
- The supplied LDOS file retains SHA-256 `a5c27b03f1d6ea45677a2d47d808d641a5c17e633a25a0965be2ac6adf6ae635`.
- Original `.nb`, `.wdx`, decoded source and archived image hashes are unchanged.
- All 30 diagram crops retain exactly the pixels of their recorded source rectangles.
- The reconstructed solver, normalization changes, phase replacement, plotting wrappers, examples and associated numerical tests have been removed from the current tree.

No numerical execution or historical-figure reproduction is claimed. The notebook explains the supplied code and the visible diagrams. Missing figure-producing routines are identified as missing.

## Run the source checks

```bash
python -m pip install -r requirements.txt
python tools/build_notebooks.py
python tools/verify_sources.py
```
