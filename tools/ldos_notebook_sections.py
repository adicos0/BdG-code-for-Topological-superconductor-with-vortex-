"""Explain contiguous excerpts from the user's unchanged LDOS file."""
from pathlib import Path
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
SOURCE = 'original_code/vortex_ldos_original.m'


def md(text):
    return nbf.v4.new_markdown_cell(text.strip())


def original(title, start, end=None):
    text = (ROOT/SOURCE).read_text()
    begin = text.index(start)
    stop = len(text) if end is None else text.index(end, begin)
    code = text[begin:stop].rstrip()
    first = text[:begin].count('\n') + 1
    last = first + code.count('\n')
    cell = md(f'**Original code: {title}** — `Vortex Ldos copy.m`, lines {first}–{last}.\n\n```mathematica\n{code}\n```')
    cell.metadata['supplied_code'] = {'path': SOURCE, 'start': begin, 'end': begin+len(code),
                                      'first_line': first, 'last_line': last}
    return cell


def sections():
    return [md('''## 1. The file being explained

The code below is copied directly from your supplied [`Vortex Ldos copy.m`](https://github.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-/blob/main/original_code/vortex_ldos_original.m). The stored file is unchanged. Each code block gives its original line numbers.

The explanations describe the assignments, indices, loops and plotting commands that are present. They do not add another solver, change the matrix convention, replace the LDOS expression, or supply missing plotting routines. The saved diagrams appear after the source walkthrough. A matching figure-producing run is not established simply by placing code beside an image.

The `.ipynb` is a reading notebook. The displayed code is Mathematica, and has not been executed in this editing pass.'''),
    md('''## 2. Set the directory and calculation parameters

`SetDirectory` selects the directory used by the later relative imports and exports. Its path is the one saved in your file.

| Name | Value in the file | How it is used |
|---|---|---|
| `kernum` | 16 | Passed to `LaunchKernels`. |
| `Ns` | 10 | Used in both momentum loops and in the export divisor `Ns^2`. |
| `kkmx`, `kkmy` | `Ns*Nx`, `Ns*Ny` | Set the momentum increments and upper limits. |
| `wmax` | 0.10 | Sets the endpoints of the sampled energy interval. |
| `eta` | 0.005 | Appears as `I*eta` in the two spectral denominators. |
| `nw` | 400 | Divides the energy interval into 400 steps. |
| `ival` | Integers 1 through 200 | Selects the components used to calculate the site-resolved arrays. |

`Nx` and `Ny` are assigned later in the file. These earlier assignments therefore use whatever values or symbolic expressions those names have when evaluated.'''),
    original('directory and sampling parameters', 'SetDirectory[', 'm = 0.50;'),
    md('''## 3. Define the hopping entries and site coordinates

`t[dx][dy]` assigns a number to each listed displacement. `l[dx][dy]` stores the four complex entries shown below. The script then sets `Nx=20`, `Ny=20` and increments `i` through a nested loop, storing `xpos[i]` and `ypos[i]`.

The inner loop runs through y first, so the index is `i=(x-1)*Ny+y`. After the loops, `NN=400`. In this ordering site 184 is `(10,4)`, and site 190 is `(10,10)`.

In this LDOS file, `t` is subsequently used to test the boundary connections. The normal Hamiltonian blocks themselves are imported from `.wdx` files; this section does not rebuild them from these amplitudes.'''),
    original('hopping entries and coordinates', 'm = 0.50;', '\n\n2\n'),
    md('''## 4. Construct the x-boundary array

The nested `If` expressions check whether `t` has a numerical entry for each listed shifted displacement. A connection through an x-shifted copy receives `-1.` or `+1.` according to the branch shown. Cases that do not match a numerical entry do not become stored numerical entries in the final sparse array.

`pxij3[i][j] = pxijmn3[[j,i]]` transposes the temporary table. `Phasex` keeps its numerical entries. `Clear` removes the intermediate objects after they have been used.'''),
    original('Phasex', 'Print["reading phasex3"]', 'Print["reading phasey3"]'),
    md('''## 5. Construct the y-boundary array

This block follows the same pattern, now distinguishing connections through the shifted y copies. `pyijmn3` is the temporary table, `pyij3` holds its transposed entries, and `Phasey` is the sparse array retained for the phase factors.'''),
    original('Phasey', 'Print["reading phasey3"]', 'Print["Reading H matrix now"]'),
    md('''## 6. Import the saved matrices

The five `Import` calls load `Huu.wdx`, `Hud.wdx`, `Hdu.wdx`, `Hdd.wdx` and `SCs.wdx` from the selected directory. The filenames alone do not specify which impurity or vortex calculation created their contents.

The following `empty` assignment is reproduced exactly. Its condition applies `NumericQ` to `2==2`, which evaluates to the Boolean `True`; that Boolean is not a number. Thus the rule does not insert the value 2 at matrix positions, and the unspecified sparse entries are zero.'''),
    original('matrix imports and empty block', 'Print["Reading H matrix now"]', 'Clear[kx, ky]'),
    md('''## 7. Assemble H as written in the file

The first four assignments multiply the imported normal blocks element by element by the exponential boundary factor. `Re[kx]` and `Re[ky]` take real parts. `hh1` puts those four blocks into a larger matrix with `ArrayFlatten`.

`SCs1` takes rows and columns `1` through `NN`; `SCs2` takes rows and columns `NN+1` through `2 NN`. These are the two diagonal portions of imported `SCs`. Both `SCs1c` and `SCs2c` are then calculated from **`SCs1`**, exactly as the file states. Although `SCs2` is assigned, it is not used in those two expressions.

`hh2` uses the displayed minus signs and complex conjugates. `huv` uses `SCs1c` and `SCs2c` in its off-diagonal blocks. Finally `H` combines `hh1`, `huv`, `Conjugate[huv]` and `hh2`.

`H == ConjugateTranspose[H];` is an equality expression with its output suppressed by the semicolon. It does not stop the calculation if equality fails. This walkthrough leaves all these expressions unchanged; it does not substitute another block convention.'''),
    original('phase-dependent matrix assembly', 'Clear[kx, ky]', '(*Clear[Huu,SCs,Hopping12];*)'),
    md('''## 8. Diagonalize at each momentum and form the LDOS array

`CloseKernels[]` closes the parallel kernels, and `LaunchKernels[kernum]` requests the specified number. `ParallelTable` runs the body for the displayed `kx` and `ky` values. The two `For` lines inside comments are inactive.

Inside the table:

1. `ES = Eigensystem[H]` obtains the eigenvalues and eigenvectors.
2. `elist` copies all `4*NN` eigenvalues from `ES[[1]]`.
3. `ulist` takes squared magnitudes of the eigenvector components indexed by `ival`.
4. `vlist` takes squared magnitudes of the components indexed by `ival+2*NN`.
5. For each `w`, the inner sum adds the two displayed denominator terms for every eigenstate, then `Im` takes the imaginary part.

The sampled energy is `-wmax+2*wmax*w/nw`. With `w=0,...,400`, the table contains 401 energies from -0.10 to +0.10. Each energy row contains a value for each of the 200 selected indices.

`Flatten[...,1]` combines the two momentum-table levels into a list of LDOS arrays. The outer `Total` sums those arrays. The result is stored in `ldos`; the normalization factor is applied later when writing the data.'''),
    original('original parallel spectral sum', 'Print["starting ldos calc"]', '      (*Print["kx="'),
    md('''## 9. Write the result to a file

The commented `Export` alternative is inactive. The active code constructs `fl`, deletes the named existing output, opens an append stream, and writes one row for every energy index `n` and selected-site index `i`.

Each row contains the energy and `(-1/Pi)*ldos[[n,i]]/Ns^2`. The site index is not written as a separate column. Because the site loop is inside the energy loop, rows for all selected sites at one energy are written before moving to the next energy.

The filename expression uses the Greek symbol `\\[Eta]`, while the broadening used above is the separately named `eta`. The file is explained as written; neither symbol is replaced here.'''),
    original('original export', '  (*Export["ldos_sup.dat"', '(*Exit[];*)'),
    md('''## 10. Read and plot the exported table

`Import[...,"Table"]` loads the specific filename shown. `ListPlot[Lds, Joined->True]` joins the imported points. This command does not separate the interleaved rows into individual site curves.

`FrameLabel` supplies the displayed labels. `GridLines` places the two red reference lines at -0.04 and +0.04, and `Epilog` places the two text annotations. These values are explicitly supplied plotting settings; this command does not calculate them from the eigenvalues.'''),
    original('import and energy-labelled plot', 'Lds= Import[', '(*ListPlot[Lds[[12500'),
    md('''## 11. Read the final site-selection plot

The active expression selects `ldos[[n,184]]` while `n` runs from 1 to `nw`. Since `ival` is the sequence 1 through 200, column 184 corresponds to selected site 184 here. The alternative column 1 expression is commented out.

This `ListPlot` receives values without explicit x coordinates. Its horizontal coordinate is therefore the list index. The loop varies the energy-row index while holding the selected site fixed; it is not a scan through lattice sites. It uses 400 rows, omitting the last of the 401 energy rows calculated above, and applies `-1/Pi` without the `/Ns^2` factor used by the export.

This explains this particular command. It does not establish that it generated the saved three-colour impurity diagrams.'''),
    original('final active plot', 'ListPlot[ {Table[(-1/Pi)*ldos[[n,184]]'),
    md('''## 12. Commented plotting alternatives

The file also contains the commented `ListDensityPlot` expression below. Because it is inside `(* ... *)`, it does not execute. Its frame labels are `Position` and `Energy`. It is not an active two-dimensional or three-dimensional vortex-core plotting routine.'''),
    original('inactive density-plot alternative', '(*LDos1 = Import[', 'ListPlot[ {Table[(-1/Pi)*ldos[[n,184]]')]
