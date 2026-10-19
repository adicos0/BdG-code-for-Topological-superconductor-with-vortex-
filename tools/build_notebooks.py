"""Build five direct explanations of the original code and existing diagrams."""
from pathlib import Path
import json
import nbformat as nbf
from ldos_notebook_sections import sections as ldos_sections

ROOT = Path(__file__).resolve().parents[1]
REPO = 'https://github.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-'
RAW = 'https://raw.githubusercontent.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-/main'
PIN = 'd88d12cded1a954d53c2d36d515525c8100ded1c'
SOURCES = {stem: json.loads((ROOT/'original_code'/f'{stem}.json').read_text())
           for stem in ['bdg_topological', 'added_vortex']}


def md(text):
    return nbf.v4.new_markdown_cell(text.strip())


def intro(number, title, description):
    return [md(f'''# {number}. {title}

{description}

Read each explanation together with the code directly below it. The code is original Mathematica; the `.nb` files are the executable versions. Diagrams are existing saved results. These reading notebooks do not rerun the calculation.''')]


def excerpt(stem, index, title, start=None, end=None):
    original = SOURCES[stem]['cells'][index-1]
    whole = original['code']
    begin = 0 if start is None else whole.index(start)
    stop = len(whole) if end is None else whole.index(end, begin)
    snippet = whole[begin:stop]
    cell = md(f'''**{title}** — `{SOURCES[stem]['source_file']}`, input {index}.

```mathematica
{snippet}
```''')
    cell.metadata['original_code'] = {'file': SOURCES[stem]['source_file'],
        'source_key': stem, 'input': index, 'start': begin, 'end': stop,
        'source_line': original['line']}
    return cell


def figure(name, title, explanation):
    cell = md(f'''### {title}

![{title}]({RAW}/diagrams/{name}.png)

{explanation}

[Enlarge the diagram]({RAW}/diagrams/{name}.png)''')
    cell.metadata['diagram'] = {'asset': f'diagrams/{name}.png'}
    return cell


def save(filename, cells):
    for index, cell in enumerate(cells):
        cell.id = f'explain-{index:03d}'
    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata = {'edition': 'Code and diagram explanations', 'original_source_commit': PIN}
    nbf.validate(nb)
    nbf.write(nb, ROOT/'notebooks'/filename)


cells = intro('01', 'Lattice and hopping matrices',
              'This notebook explains how the program numbers sites, assigns hopping amplitudes and builds the normal-state matrix blocks.')
cells += [md('''## 1. Understand the lattice diagram'''),
          figure('lattice_connections', 'Sites and connections on a square grid',
                 '''The circles are lattice sites. Horizontal and vertical arrows indicate axial connections; diagonal arrows indicate connections that change both coordinates. A hopping matrix records these connections as entries between pairs of sites.

The numbers in this illustration are schematic labels. The program below uses its own indexing rule, with the `y` loop inside the `x` loop. Use that rule when translating a coordinate into an array index.'''),
          md('''## 2. Set the amplitudes

`m = 0.50` controls the scalar hopping dictionary `t`. The entry `t[dx][dy]` is the amplitude assigned to displacement `(dx,dy)`:

| Displacement | Entry | Value when m=0.5 |
| --- | --- | --- |
| Same site `(0,0)` | `3 m/2` | 0.75 |
| One step horizontally or vertically | `-m` | -0.5 |
| One step along each coordinate, diagonally | `m/4` | 0.125 |

`lam = 1` controls the separate spin-mixing dictionary `l`. All four axial entries receive `lam (1-I)`, where `I` is Mathematica's imaginary unit. These values later form `soc`. `SetDirectory` selects the folder used for subsequent file exports; its saved path is machine-specific.'''),
          excerpt('bdg_topological',1,'Hopping parameters',end='Nx = 20;'),
          md('''## 3. Number the sites

`Nx` and `Ny` are the numbers of sites in the two directions. The loops increment `i` once for each coordinate pair and store that pair in `xpos[i]` and `ypos[i]`. Thus `NN = Nx Ny = 400`.

With Mathematica's one-based indexing, the code gives `i = (x-1) Ny + y`. For `Ny=20`, `(1,1)` is site 1, `(1,20)` is site 20, and `(10,10)` is site 190. `empty` is a 400 × 400 sparse zero matrix. `Vimp` is assigned zero here, but assigning it alone does not insert an impurity into a matrix.'''),
          excerpt('bdg_topological',1,'Coordinates and size',start='Nx = 20;',end='(*Hopping Matrix*)'),
          md('''## 4. Build the scalar hopping matrix H1

For each pair `(i,j)`, the nested `If` statements inspect coordinate differences with shifts of `-Nx`, `0`, `+Nx` and `-Ny`, `0`, `+Ny`. These shifted copies let a connection cross an edge and enter from the opposite edge: periodic boundaries.

`NumericQ` asks whether the dictionary contains a numerical amplitude for a tested displacement. The first matching branch supplies the entry. `tij[i][j] = tijmn[[j,i]]` transposes the temporary table. The final `SparseArray` keeps numerical entries and makes the unspecified entries zero. The result `H1` has size `NN × NN`.'''),
          excerpt('bdg_topological',1,'Periodic scalar hopping',start='(*Hopping Matrix*)',end='tijmn = Table[\n   If[NumericQ[l'),
          md('''## 5. Build the spin-mixing matrix soc

This section repeats the same periodic-image lookup, now using `l` instead of `t`. The result `soc` connects the two spin groups. Its entries are complex because `l` contains `1-I`. `ConjugateTranspose[soc]` means transpose the matrix and reverse the sign of each imaginary part.'''),
          excerpt('bdg_topological',1,'Periodic spin mixing',start='tijmn = Table[\n   If[NumericQ[l',end='Clear[mu]'),
          md('''## 6. Read the chemical-potential and block definitions

`Chempot1[mu]` is diagonal with entries `-mu`. Consequently, `H1-Chempot1[mu]` contains `H1 + mu I`.

`Hop1x[mu]` assembles four `NN × NN` blocks into a `2 NN × 2 NN` matrix. Its upper and lower scalar blocks have opposite `H1` signs. **Execution detail:** the iteration does not call `Hop1x`; it separately defines `Huu` and `Hdd`, both with `H1-Chempot1[mu]`. Notebook 02 follows that active assembly.'''),
          excerpt('bdg_topological',1,'Chemical potential and Hop1x',start='Clear[mu]')]
save('01_lattice_and_normal_hamiltonian.ipynb',cells)

cells = intro('02','BdG matrix, eigenvalues and spectrum plots',
              'Follow the matrix assembly from its smaller blocks to Eigensystem, then read the energy-spectrum diagrams.')
cells += [md('''## 1. Build the onsite pairing block

The nested loops construct `js`: for every site `i`, they find the site with the same `xpos` and `ypos`. With unique coordinates, `js[[i]]` contains `i` itself. This is why the pairing update is onsite.

`st` stores one pairing value per site. `Flatten[st]` turns its nested list into a vector. `SCs2[st]` places that vector on the diagonal. `SCs1` is a zero block. `SCs[st]` places `SCs2[st]` in the two off-diagonal spin blocks, giving a `2 NN × 2 NN` pairing matrix. `SCs3` and the helper `s[i_,j_]` are defined here but are not used in the active iteration shown next.'''),
          excerpt('bdg_topological',2,'Onsite pairing definitions'),
          md('''## 2. Assemble the full matrix

For `NN=400`, each spin block is 400 × 400, `hh1` and `hh2` are 800 × 800, and `H` is 1600 × 1600.

| Variable | What the active code puts into it |
| --- | --- |
| `Huu`, `Hdd` | Same-spin hopping plus the diagonal chemical-potential shift. Both use the same expression here. |
| `Hud` | The spin-mixing matrix `soc`. |
| `Hdu` | The complex conjugate of `soc`. |
| `hh1` | The four electron spin blocks joined into one matrix. |
| `hh2` | The corresponding hole blocks with the signs and conjugations displayed in the code. |
| `Huv` | The pairing block `SCs[st]`. |
| `Hvu` | Its complex conjugate. |
| `H` | The complete electron–hole block matrix. |

`ArrayFlatten` joins blocks; it does not diagonalize them. `Conjugate` reverses the sign of the imaginary part. The exact signs are retained below.'''),
          excerpt('bdg_topological',5,'Matrix assembly and diagonalization',start='  Huu =',end='  fermi ='),
          md('''## 3. What Eigensystem returns

`ES = Eigensystem[H]` returns two lists. `ES[[1]]` contains the eigenvalues, or energies. `ES[[2]]` contains the corresponding eigenvectors, stored as a list of vectors. Each vector has `4 NN` components.

The calculation solves `H . vector = energy vector`. An eigenvalue says where a level lies in energy. Its eigenvector says how the state is distributed across sites and component groups. Plotting the eigenvalues alone does not display that spatial distribution.

In the later occupation sums, `ES[[2]][[;;,i]]` selects component `i` from every eigenvector. The offsets `i+NN`, `i+2NN` and `i+3NN` select the same site in the other three component groups.'''),
          md('''## 4. How ListPlot draws the spectrum

`Sort[ES[[1]]]` puts the energies in ascending order. Since `ListPlot` receives values rather than explicit coordinate pairs, it uses their positions in this sorted list as the horizontal coordinate. The horizontal `n` is therefore a **state index**, not a lattice position or momentum.

`PlotRange` selects the visible energy interval. `Frame` draws a frame; `FrameLabel` labels its axes. `GridLines` places horizontal reference lines at the specified energies, and `Directive[Green,Dashed]` styles them. The lines are chosen plotting parameters, not automatically detected gap edges. `PlotLabel` is saved text; it is not updated automatically from the current values of `mu` or `st`.'''),
          excerpt('bdg_topological',7,'Spectrum plotting command'),
          excerpt('added_vortex',7,'Second saved plotting command'),
          md('''## 5. Read the saved spectrum diagrams

In each diagram, blue points are energy levels and the vertical axis is energy in the displayed units. Several levels may overlap at the resolution of the image. These are saved examples, so their labels and visible ranges may differ from the plotting cells above.'''),
          figure('spectrum_zero_pairing','Spectrum with levels near zero',
                 'The title sets the pairing parameter to zero. The blue levels approach and cross the horizontal zero-energy reference without the wide empty central interval seen in the next diagram. The compressed horizontal spread comes from showing only part of the energy range while keeping a broad state-index axis.'),
          figure('spectrum_gap','Spectrum with an empty central interval',
                 'The green dashed lines mark approximately +0.06 and -0.06 on the energy axis. The visible points form negative- and positive-energy groups, leaving a central interval with no resolved dots. This is how a gap appears in a sorted-energy plot; its horizontal width is not a physical length.'),
          figure('spectrum_central_state','Spectrum with a level close to zero',
                 'A blue level appears close to zero energy, between the green reference lines. Other levels lie close to the two reference boundaries. This image gives their energies; it does not show whether their wavefunctions are at the centre or the edge of the lattice, and does not by itself identify the type of state.')]
save('02_bdg_matrix_and_diagonalization.ipynb',cells)

cells = intro('03','Occupation numbers and the pairing iteration',
              'This notebook explains how each diagonalization produces occupations and a proposed pairing field, and how the next iteration uses those results.')
cells += [md('''## 1. Set the Fermi function and starting values

`f[x]` gives the occupation weight for an energy `x`. The second definition replaces the first: energies above `0.02` receive weight 0, energies below `-0.02` receive weight 1, and intermediate energies use `1/(Exp[Re[x]/kT]+1)`. `Re[x]` selects the real part. `kT` sets the energy scale over which the weights change.

`vs` is the pairing coupling, `mu` is the initial chemical-potential parameter, and `nexp` is the target mean occupation per site. The two files use different values:

| File | vs | nexp | initial mu | initial ds |
| --- | ---: | ---: | ---: | ---: |
| BdG topological.nb | -1.45 | 1.2 | 0 | -0.01 |
| added vortex.nb | -1.25 | 1.0 | 0 | -0.01 |

`st=Table[ds,...]` creates one entry per site. In a fresh session, `ds` is still symbolic at this statement and receives its value on the next line. Running the cells in their original order matters. The large section enclosed in `(* ... *)` is inactive.'''),
          excerpt('bdg_topological',3,'Occupation function and initialization'),
          md('''## 2. Diagonalize using the current fields

`iter=100` fixes the number of passes through the loop. Each pass first constructs `H` using the current `mu` and `st`, then evaluates `Eigensystem[H]`. The program does not stop early when the fields change by a small amount. The next excerpt shows the beginning of that loop; its closing bracket appears in the update section below.'''),
          excerpt('bdg_topological',5,'Beginning of the iteration',end='  fermi ='),
          md('''## 3. Calculate occupation and spin-population differences

`Map[f, ES[[1]]]` creates one occupation weight for every eigenvalue. `Abs[component]^2` gives a nonnegative component weight, and the dot product sums that weight over all eigenstates.

`nutcalc` uses the first component group and `fermi`. `ndtcalc` uses the fourth group and `1-fermi`. `ni` adds the two contributions; `mi` subtracts them. `nt` is the mean occupation per site. `mt` averages `Abs[mi]`, so opposite signs at different sites do not cancel. These sums use the entire returned eigensystem.'''),
          excerpt('bdg_topological',5,'Occupations and averages',start='  fermi =',end='  stcalc='),
          md('''## 4. Form the proposed pairing field

For each site, the code multiplies components from different groups and complex-conjugates the second factor. It adds the two products, weights them by `fermi`, sums over states, and multiplies by `vs`.

`Inner[Times,...,List]` creates the element-by-element products. The following dot product with `fermi` performs their weighted sum. Since `js[[i]]` contains the same site, this produces one onsite proposal per site in `stcalc`.'''),
          excerpt('bdg_topological',5,'Pairing proposal',start='  stcalc=',end='  (*\nvsite ='),
          md('''## 5. Print progress and update the fields

`sta` is the arithmetic mean of the complex entries of `stcalc`. The printed `Abs[sta]` is the magnitude of that mean; it is not the mean of the individual magnitudes.

The chemical-potential update is `mu = mu - 0.1 (nexp-nt)`. If `nt` is below the target, this lowers `mu`; if it is above the target, it raises `mu`. In this code the electron diagonal contains `+mu`, which explains this feedback direction.

The pairing update `st = 0.5 stcalc + 0.5 st` keeps half the previous field and half the proposal. This mixing moderates changes between iterations. A fixed count of 100 updates does not demonstrate convergence by itself.'''),
          excerpt('bdg_topological',5,'Progress output, updates and exports',start='  sta='),
          md('''## 6. Understand the exported arrays

`Huu.wdx`, `Hud.wdx`, `Hdu.wdx` and `Hdd.wdx` store the spin blocks left from the last matrix construction. `SCs.wdx` stores `SCs[st]` after the final pairing update. Because the last update follows the last diagonalization, the saved pairing field and saved normal blocks are not all recomputed together in a final extra pass.

The `Print` line reports the iteration index, current `mu`, `mt`, `nt`, and `Abs[sta]`. These values are progress information. There is no active convergence-curve plotting command in these cells, so this notebook does not add a convergence diagram.''')]
save('03_self_consistent_pairing.ipynb',cells)

cells = intro('04','Vortex code and spatial diagrams',
              'Read the commented vortex block term by term, then interpret the line profiles, colour maps, surfaces and energy levels.')
cells += [md('''## 1. What the vortex block does

This block is enclosed in `(* ... *)`, so it does not execute in the saved program.

| Name or expression | Role in the block |
| --- | --- |
| `vsite` | Index of the site chosen as the centre. |
| `xpos[vsite]`, `ypos[vsite]` | Coordinates of that centre. |
| `l` and `zetanot=l/Pi` | Inputs setting the length scale in the profile. |
| `Sqrt[(x-x0)^2+(y-y0)^2]` | Distance of the current site from the centre. |
| `Tanh[distance/zetanot]` | Amplitude multiplier: zero at the centre, approaching one far away. |
| `phi` | The position-dependent angle computed by the `If`/`ArcTan` statements. |
| `Exp[I*phi]` | A complex factor with unit magnitude that supplies a phase. |

The block multiplies the pairing proposal by the amplitude and phase factors. Its `ArcTan` uses absolute coordinate differences, which discard quadrant signs; when the x coordinates coincide it assigns `Pi/2`. That is the code's actual phase rule. It should not be read as a demonstrated full phase winding.

With the indexing from notebook 01, `vsite=378` on a 20 × 20 lattice means `(19,18)`. A diagram centred at `(10,10)` therefore does not use that saved centre value. The diagrams below are existing outputs, not new output from this disabled block.''')]
for stem in ['bdg_topological','added_vortex']:
    full=SOURCES[stem]['cells'][4]['code'];vpos=full.index('vsite =');begin=full.rfind('(*',0,vpos)
    cells.append(excerpt(stem,5,'Commented vortex block — '+SOURCES[stem]['source_file'],
                         start=full[begin:vpos+len('vsite =')],end='  (*spmtcalcnn'))
cells += [md('''## 2. Read the spatial plots

A line profile shows values along a one-dimensional cut. A colour map places a value at each `(x,y)` grid cell. A surface shows spatial variation through a height coordinate. These are different views; the axes and scales must be read separately.'''),
          figure('pairing_line_profile','Pairing amplitude along a cut',
                 'The horizontal coordinate is labelled r and runs across the lattice cut. The vertical axis is the pairing-amplitude label Δ₀. The curve is approximately 0.03 away from the centre and falls close to zero at r≈10. Its colour varies along the same line; it is not a separate magnetic-field measurement.'),
          figure('pairing_map_smooth','Smooth spatial colour map',
                 'Both axes are lattice coordinates, from 1 to 20. The pale central cell near (10,10) is surrounded by lighter horizontal and vertical bands, while the corners are darker orange. The image has no numerical colour bar, so colours identify spatial differences but do not provide exact amplitude values. The y coordinate increases downward in this map.'),
          figure('pairing_surface_smooth','Surface with a central depression',
                 'The two horizontal axes locate lattice sites. The third axis gives the plotted pairing-profile value. Most of the surface lies around 0.03–0.04, with a narrow central depression toward zero. The mesh follows the sampled lattice; the rendered surface between mesh points does not add extra calculated sites.'),
          figure('pairing_map_structured','Spatial map with a narrow central feature',
                 'The central pale cell remains distinct, while surrounding cells show a more irregular pattern than the smoother map above. Each square represents a lattice cell. Without a colour bar, an exact numerical comparison between the two maps is not available from colour alone.'),
          figure('pairing_surface_structured','Nearly flat outer surface and a narrow core',
                 'The outer surface is comparatively flat and the central depression is sharp. The vertical tick labels place the bottom near zero and the outer region near 0.05. This diagram describes spatial amplitude variation; it does not display the complex phase.'),
          md('## 3. Read the energy plot'),
          figure('spectrum_in_gap','Energy levels between the reference lines',
                 'The horizontal n is a sorted eigenvalue index and the vertical axis is energy. The green dashed lines mark approximately ±0.05. Blue dots inside those lines show levels within that chosen interval. The index axis does not locate these states in the lattice; spatial localization would require their eigenvectors.')]
save('04_vortex_profile_and_spectrum.ipynb',cells)

cells = intro('05','Impurity code and LDOS diagrams',
              'This notebook explains the impurity-related code that is present and gives a direct reading of each saved plot: its axes, curves, colours and visible changes.')
cells[1:] = []
cells[0] = md('# 05. Impurities, LDOS and vortex-core views\n\nRead the original LDOS code section by section, followed by the saved diagrams. Every displayed code block is copied from a supplied source file.')
cells += ldos_sections()
cells += [md('''## 13. What the original files contain about impurities

In `BdG topological.nb`, `Vimp=0` is assigned, but the displayed `+imp` contribution is commented out. The displayed active spin-block assignments therefore do not run a non-magnetic or magnetic impurity sweep. The LDOS file loads saved matrices rather than adding those impurity terms itself. The missing sweep code is not reconstructed here.'''),
          excerpt('bdg_topological',1,'Original impurity parameter',start='Vimp = 0;',end='empty='),
          excerpt('bdg_topological',5,'Original spin-block assignments',start='  Huu =',end='    Huv='),
          md('''## 14. The vortex update present in added vortex.nb

This is the actual vortex-related block in `added vortex.nb`. It sets `vsite=378`, `l=0.3` and `zetanot=l/Pi`. It calculates `phi` using the displayed absolute coordinate differences, forms the eigenvector-product sum, multiplies by `vs`, and multiplies by the `Tanh` and `Exp[I*phi]` factors to form `stcalc`.

`xpos` and `ypos` refer to the coordinate arrays built earlier in that file. Its initial lattice assignments are `Nx=28` and `Ny=28`; these settings must not be silently replaced by the settings printed on a different saved image.

The `Tanh` argument uses the distance from the selected site. `Exp[I*phi]` multiplies the result by a complex phase. The expressions below are left exactly as supplied, including the absolute-value angle formula. This is a field-update block; it contains no heatmap or surface-plot command. The provided files do not establish the exact commands that produced the saved two-dimensional and three-dimensional core views.'''),
          excerpt('added_vortex',5,'Original vortex contribution to stcalc',start='vsite = 378;',end='  (*Print['),
          md('''## 15. The original spectrum plotting command

This command plots `Sort[ES[[1]]]`, the sorted eigenvalues. The title, range and reference lines are specified explicitly. It is copied from the original file; those settings do not demonstrate that this command produced every archived spectrum shown below.'''),
          excerpt('bdg_topological',7,'Original spectrum plot')]

cells += [md('## Saved energy-axis LDOS diagrams\n\nThe following are archived outputs. Their original energy labels and scales are retained; the exact generating run for each image has not been verified from the supplied files.'),
          figure('ldos_central_peak','LDOS with a narrow central peak',
                 'The blue curve has a large peak close to zero energy, with many smaller oscillations on either side. The curve remains appreciable through the central region. The small oscillations are features of the saved curve; the image alone does not establish their numerical origin.'),
          figure('ldos_gap','LDOS suppressed between two sharp peaks',
                 'Red dotted reference lines lie at approximately -0.06 and +0.06 meV. Between them, the LDOS is close to zero; near either line it rises sharply. The two sides are not identical in height and background. These plotted values should be read directly rather than assumed symmetric.'),
          figure('ldos_two_inner_peaks','Two peaks inside the marked energy interval',
                 'The red dotted lines mark approximately ±0.05 meV. Two broad maxima appear inside that interval, with a smaller but nonzero LDOS near zero energy. A central dip and two finite-energy peaks are different from a zero-energy peak.'),
          md('''## Saved pairing line profiles

The three panels below plot pairing amplitude against position. Their horizontal coordinate is spatial, unlike the energy axes immediately above.'''),
          figure('pairing_profiles_impurity','Three spatial profiles with different local distortion',
                 'All three orange curves approach zero near the central coordinate. From left to right, the neighbouring region becomes more uneven, while values far from the centre stay around the same general level. The panel scales are slightly different: the rightmost vertical scale reaches about 0.08, whereas the first two reach about 0.07. The colours identify the plotted line, not three distinct observables.'),
          md('''## Saved impurity LDOS curves

The next six saved images label the horizontal axis **Number of Sites**. That is the archived label. The exact generating command for these three-colour curves is not supplied. The active final command above selects column 184 and contains one active curve; it is not evidence for the site identities of the green, blue and orange curves in these images.

The panels are grouped by their stored scalar/nonmagnetic or magnetic case. The descriptions below state what is visible; they do not assume that changes in peak height prove protection or identify a particular bound state.'''),
          figure('ldos_scalar_weak','Scalar case: green curve above the other two',
                 'The green trace lies above the blue and orange traces through much of the index range. All three curves show a deep minimum around index 100, with larger values on both sides. The green minimum remains slightly higher than the others.'),
          figure('ldos_scalar_middle','Scalar case: a reduced green trace',
                 'The green curve is much smaller than the other two over most indices. The blue and orange curves retain a central depression and prominent shoulders. The orange trace is higher around several shoulder peaks, so colour order is not fixed across all cases.'),
          figure('ldos_scalar_strong','Scalar case: green trace close to the baseline',
                 'The green trace is nearly flat at zero. Blue and orange retain a substantial signal away from the central dip, with the blue trace generally larger. The result is strongly position dependent; it is not a uniform rescaling of every curve.'),
          figure('ldos_magnetic_weak','Magnetic case: three resolved LDOS traces',
                 'All three traces have a central trough. The blue curve has the largest shoulder peaks, while green and orange cross in parts of the outer region. Compare the traces at the same horizontal index to assess their relative values.'),
          figure('ldos_magnetic_middle','Magnetic case: suppression of the green trace',
                 'The green curve is low across the plotted range. Blue and orange show pronounced shoulders surrounding the central trough. The blue shoulders are taller, but the two outer traces remain similar in some regions.'),
          figure('ldos_magnetic_strong','Magnetic case: green trace near zero',
                 'The green result stays close to the baseline. The remaining two traces still have a broad central minimum, with blue generally above orange. The vertical scale differs from other panels, so identical apparent heights on screen need not represent identical values.'),
          md('''## Saved impurity spectra

In each spectrum below, horizontal n is the index in a sorted energy list. Blue dots are eigenvalues. Green dashed lines mark the reference values drawn in each image. The original spectrum command above also specifies reference lines, but its settings are not verified as the settings for every image.'''),
          figure('spectrum_scalar_a','Scalar spectrum: levels inside the reference interval',
                 'Several blue dots lie between the green lines around zero. More distant levels occur above and below that interval. This plot shows energies, not the locations of their wavefunctions.'),
          figure('spectrum_magnetic_a','Magnetic spectrum: small positive and negative levels',
                 'Blue levels occur on both sides of zero, including inside the marked interval. Closely spaced points can overlap at this image resolution, so an exact degeneracy count should not be inferred from visible dots alone.'),
          figure('spectrum_scalar_b','Scalar spectrum with a near-zero level',
                 'A dot lies close to zero and other levels lie near the reference boundaries. The title records the parameters attached to this saved plot. A near-zero energy alone does not establish the state’s spatial distribution or physical identity.'),
          figure('spectrum_magnetic_b','Magnetic spectrum with a near-zero level',
                 'The same types of features are visible: a central near-zero level, nearby levels and an outer sequence of states. The broad resemblance to another spectrum does not establish equality of the underlying matrices or wavefunctions.'),
          md('''## Saved spatial maps and surfaces

Colour maps place one cell at each lattice coordinate. The surface diagrams also show a vertical numerical axis. Several of these surfaces have **zero at the top of the vertical axis and larger positive values lower down**. An apparent upward spike can therefore correspond to a smaller value; read the ticks before describing it as an enhancement.'''),
          figure('pairing_map_scalar_weak','Scalar case: a distinct central cell',
                 'A pale cell near (10,10) stands out from the surrounding orange and yellow cells. The coordinates span a 20 × 20 grid, with y increasing downward. There is no colour bar, so exact amplitudes cannot be recovered from the colours.'),
          figure('pairing_surface_scalar_a','Scalar surface with an inverted vertical scale',
                 'The surface is uneven and a sharp central feature rises toward the tick marked zero. Larger positive tick values are lower in the drawing. Thus its visual upward direction must not be interpreted as an increase of the plotted quantity.'),
          figure('pairing_map_scalar_strong','Scalar case: irregular surrounding cells',
                 'The central region is pale, while surrounding cells fluctuate more strongly between yellow and orange. The map displays spatial structure but supplies no numerical colour scale. It does not show the phase of a complex field.'),
          figure('pairing_surface_scalar_b','Scalar surface with a central feature and rough surroundings',
                 'The vertical ticks run from about 0 near the top toward 0.08 below. The central feature points toward smaller values, while the surrounding surface is uneven. Its vertical range differs from the preceding scalar surface.'),
          figure('pairing_surface_magnetic','Magnetic surface with a pronounced central feature',
                 'The vertical scale again increases downward, from 0 toward about 0.06. The central feature reaches toward the smaller-value end of that scale. Away from the centre, the mesh has substantial site-to-site variation.'),
          figure('pairing_map_magnetic','Magnetic case: patchy spatial pattern',
                 'Orange and yellow cells alternate across much of the grid, and a pale central cell remains visible. The plot gives spatial coordinates but no numerical colour key. It supports a description of the pattern, not a precise amplitude comparison with another colour map.')]
save('05_impurities_and_ldos.ipynb',cells)
print('Built five code-and-diagram guides.')
