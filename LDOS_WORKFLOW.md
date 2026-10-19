# Reading the original LDOS code

Notebook 05 explains `original_code/vortex_ldos_original.m`, an unchanged copy of the supplied `Vortex Ldos copy.m`. Every displayed LDOS block has source line numbers and exact character-range metadata.

## What is explained

1. Directory, kernel count, selected sites and sampling parameters.
2. Hopping entries and coordinate loops.
3. The original x and y boundary arrays.
4. Imports of the five saved matrix files.
5. The phase factors and matrix blocks exactly as written.
6. `Eigensystem`, the component selections and the original two-term spectral sum.
7. Export order, `-1/Pi`, and the `Ns^2` divisor.
8. The table plot, its explicit reference lines, and the final fixed-column plot.
9. The inactive density-plot alternative.
10. Original impurity assignments, the vortex update in `added vortex.nb`, and the original spectrum plotting command.

## Limits established by the supplied files

- The LDOS script imports matrices; it does not implement a scalar/magnetic impurity sweep.
- The displayed impurity addition in `BdG topological.nb` is commented out.
- The vortex update computes `stcalc`; it is not a heatmap or surface-plot command.
- The supplied final LDOS plot varies energy-row index `n` at fixed column 184. It does not identify the three coloured curves in the archived impurity figures.
- The exact commands and run inputs for all archived figures have not been established. The figures are described visually without asserting that an unverified source block produced them.

## Preservation

The original spectral expression, component indices, matrix signs, pairing extraction, phase expression, parameters and plotting commands are unchanged. No alternative normalization or scientific calculation is introduced. The added reconstructed implementation and its examples/tests have been removed from the current repository tree.

This is a source explanation and presentation update. The original numerical program was not executed. Verification checks notebook format, exact excerpt matches, unchanged original-file hashes and unchanged diagram pixels.
