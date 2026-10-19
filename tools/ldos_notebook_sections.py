"""Source-linked teaching sections for the supplied and adapted LDOS code."""
from pathlib import Path
import hashlib
import re
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
REPO = 'https://github.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-'


def md(text):
    return nbf.v4.new_markdown_cell(text.strip())


def block(section):
    path = ROOT/'scripts/vortex_ldos.wl'
    text = path.read_text()
    marker = f'(* SECTION: {section} *)'
    start = text.index(marker)
    end = text.find('(* SECTION:', start + len(marker))
    if end < 0:
        end = len(text)
    code = text[start:end].rstrip()
    cell = md(f'**Adapted implementation — `{section}`**\n\n```mathematica\n{code}\n```')
    cell.metadata['adapted_code'] = {'path': 'scripts/vortex_ldos.wl', 'start': start, 'end': start+len(code)}
    return cell


def sections():
    cells = [md(r'''## 1. Files, scope and how to run

The supplied [`vortex_ldos_original.m`](https://github.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-/blob/main/original_code/vortex_ldos_original.m) is preserved byte for byte. The runnable adaptation is [`scripts/vortex_ldos.wl`](https://github.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-/blob/main/scripts/vortex_ldos.wl). Each adapted block below comes directly from that file.

**Verified source behaviour:** the supplied program imports five matrices, adds boundary phases, diagonalizes a matrix at each momentum, computes broadened spectral sums and exports LDOS. It does not run an impurity sweep or generate the saved spatial maps.

**Adaptations:** the new file fixes the pairing-block extraction, uses a consistent hole-momentum convention, supplies explicit energy coordinates, separates sites, uses one normalization, adds impurity potentials and supplies plot functions. An optional pairing iteration is adapted from `added vortex.nb`.

**Not yet verified:** the exact historical settings and input matrices for every archived curve. Saved images remain reference results; adapted code is not described as the original generator of those images. The historical curves include conflicting unit and axis labels, so new plots use model energy units unless a conversion is independently established.

Load the `.wl` file in Mathematica using `Get["/absolute/path/to/repository/scripts/vortex_ldos.wl"]`. Loading defines functions without starting a calculation. These `.ipynb` cells are Markdown reading cells with Mathematica code, not Python execution cells. See the runnable examples in `scripts/ldos_examples.wl`.

### Setup and site numbering

`NN = Nx Ny` is the number of sites. With the saved loop order, $i=(x-1)N_y+y$. Thus site 190 is $(10,10)$ on a $20\times20$ grid, site 180 is $(9,20)$, and site 1 is $(1,1)$. Calling site 180 “nearby” is a historical label; its actual coordinate distance from the core is not one lattice spacing. Site 1 is an explicit example for the farther curve, not a verified historical choice.

`bdgRequire` stops on a failed condition and returns a labelled failure through `Catch`. `bdgNorm` is the largest absolute matrix entry, used for numerical consistency checks.'''), block('setup'), md(r'''## 2. Import the normal and pairing matrices

`Huu` and `Hdd` are the two diagonal spin blocks; `Hud` and `Hdu` connect them. Each must be `NN × NN`. `SCs` is `2 NN × 2 NN`.

In the saved source, `SCs` has the form `{{0,D},{D,0}}`, where `D` is the diagonal matrix of onsite complex pairing values. The supplied LDOS file instead reads its two diagonal blocks; that would discard pairing for this layout. `LoadBdGModel` reads the upper-right block, verifies the entire layout, and sets `Delta = Diagonal[D]`.

The checks reject incompatible dimensions, non-numerical entries, non-Hermitian normal blocks and non-onsite pairing. A square matrix alone does not determine the two lattice dimensions. Supply `Nx` and `Ny` explicitly. Imported matrices do not carry a verified history of impurity settings.'''), block('import'), md(r'''## 3. Put an impurity into the electron blocks

Choose a site $i_0$ and let $P$ contain one at matrix position $(i_0,i_0)$ and zero elsewhere. Let $V$ be the impurity energy.

| Case | Change to `Huu` | Change to `Hdd` |
|---|---|---|
| No additional impurity | 0 | 0 |
| Scalar / non-magnetic | $+VP$ | $+VP$ |
| Magnetic, diagonal in the chosen spin basis | $+VP$ | $-VP$ |

The scalar case shifts both spin energies equally. The magnetic case shifts them with opposite signs. The hole block is constructed later from these changed electron blocks, so an independent hole potential is not added by hand.

All terms have energy units; $P$ is dimensionless. Setting $V=0$ recovers the input blocks exactly. A new case must start from an impurity-free baseline, not from the output of the previous impurity case. The runner requires `CleanBaselineConfirmed -> True` when an impurity is added; this is a data-provenance assertion, not a test performed by the code.'''), block('impurities'), md(r'''## 4. Boundary phases and the full matrix

**Conventions and assumptions:** square-lattice coordinates are counted in lattice spacings ($a=1$); boundaries are periodic with Bloch twists. The imported hopping is onsite, nearest-neighbour or diagonal nearest-neighbour as in the original dictionary. Long-range hopping needs a different image-selection rule. No spacetime metric or Planck-unit convention is used.

`BoundaryPhases` stores `-1`, `0` or `+1` for an edge-crossing displacement. It follows the transposed image-lookup orientation of the supplied program. Reversing the pair of sites reverses the sign. `Exp[I(kx Nx Phasex + ky Ny Phasey)]` multiplies entries elementwise; it is not matrix multiplication. The exponent is dimensionless because each wave number has inverse-lattice-length units.'''), block('boundary'), md(r'''### Electron block, hole block and pairing

Let $h(k)$ be the two-spin electron matrix after the boundary phases are applied, $S=\mathrm{diag}(I,-I)$, and $B=\begin{pmatrix}0&D\\D&0\end{pmatrix}$. The adaptation uses

$$H(k)=\begin{pmatrix}h(k)&B\\B^\dagger&-S h(-k)^* S\end{pmatrix}.$$

Here $I$ is the `NN × NN` identity matrix, $*$ denotes complex conjugation, and $\dagger$ denotes conjugate transpose. The first two groups of `NN` components are electron amplitudes; the last two are hole amplitudes in the stored sign convention. Multiplying the second hole group by minus one converts this to the usual antisymmetric onsite pairing block. At $k=0$ the hole-block signs match the saved matrix layout.

The `-k` argument matters: conjugation reverses the boundary phase. Building both sectors consistently gives the hole hopping the same real-space boundary twist after conjugation. This is a correction to the supplied script's use of `Conjugate[Huuc]` at the same momentum. `ConjugateTranspose[B]` enforces the Hermitian partner for complex pairing. Onsite pairing has no edge-crossing phase.

`shift` adds a common onsite energy to the electron blocks, with the sign convention used by the original chemical-potential update. Every block of $H$ has the same energy units. At zero pairing, the electron and hole sectors decouple. At zero twist, the boundary factors become one.'''), block('matrix'), md(r'''## 5. Construct or retain a vortex profile

With `SeedVortex -> False`, the imported pairing field is retained. With `SeedVortex -> True`, its largest magnitude supplies the amplitude scale and `VortexFactor` supplies a new, explicitly reconstructed seed.

For distance $r$ from the chosen core and positive length $\xi$, measured in lattice sites, the factor is

$$f(r)=\tanh(r/\xi)\exp(i\phi).$$

Both $r/\xi$ and $\phi$ are dimensionless. At the core, $r=0$ gives $f=0$. Far from the core, $\tanh(r/\xi)\to1$, so the magnitude approaches one. The complex exponential changes phase but not magnitude.

`Winding` uses `ArcTan[dx,dy]` to retain quadrant information. `OriginalAbsoluteAngle` retains the original use of absolute coordinate differences, which loses quadrant signs and does not describe a full $2\pi$ phase winding. These choices are labelled rather than silently substituted. The default length `.3/Pi` is inherited numerically from the saved code and is treated as a length in lattice sites; no nanometre conversion is assumed.

This is an imposed profile on a periodically repeated cell. It does not include magnetic Peierls phases, flux quantization, or a separately validated magnetic-periodic single-vortex construction.'''), block('vortex'), md(r'''## 6. Allow the pairing field to respond to the impurity

**Imported mode** changes the Hamiltonian while keeping the supplied pairing field fixed. Consequently its pairing maps cannot develop impurity-dependent distortions. This mode is suitable for a fixed-field spectral comparison.

**Relax mode** adapts the saved constrained-vortex iteration. At each step it:

1. Builds and diagonalizes the zero-twist matrix using the current pairing and common onsite shift.
2. Evaluates the saved clipped occupation function `f` for all eigenvalues.
3. Forms the density from the same component weights as the original solver.
4. Forms the original onsite anomalous sum, multiplies by the coupling and imposed vortex factor, and obtains `proposal`.
5. Measures the maximum difference `Abs[proposal-delta]`, mixes old and proposed fields equally, and updates the common onsite shift.
6. Stops only when both pairing and density residuals meet their tolerances, or when the iteration limit is reached.

The occupation and eigenvector weights are dimensionless. Coupling times the anomalous sum has energy units, as do `proposal`, `delta`, the pairing tolerance and onsite shift. The coefficient `.1` in the shift update carries the model energy scale. Density tolerance is dimensionless. The occupation cutoff `.02` and `kT` use the same energy unit as the matrix.

This retains the saved update convention; it is not a new derivation of a fully unconstrained vortex solution. The vortex factor is imposed at every iteration. The momentum average is applied afterward to LDOS, as in the supplied workflow. Reaching 100 iterations is not evidence of convergence. `RelaxPairing` returns the iteration history even if unconverged; the high-level runner refuses to label an unconverged result as a completed case. The final spectrum and LDOS are rebuilt from the same final field and onsite shift.'''), block('relaxation'), md(r'''## 7. Calculate normalized LDOS

`Eigensystem` returns an energy list and a list of eigenvectors. The code uses **all eigenvalues**, and the electron weights of both spin groups. This avoids adding a second hole contribution to an already complete electron spectral sum.

For eigenstate $n$ and site $i$, define its spin-averaged electron weight

$$w_{ni}=\frac{|u_{n i \uparrow}|^2+|u_{n i \downarrow}|^2}{2}.$$

To see the broadening used in the supplied code, set $x=E-E_n$ and multiply numerator and denominator by the complex conjugate:

$$\frac{1}{x+i\eta}=\frac{x-i\eta}{(x+i\eta)(x-i\eta)}=\frac{x-i\eta}{x^2+\eta^2},$$
$$-\frac1\pi\operatorname{Im}\frac{1}{x+i\eta}=\frac{\eta}{\pi(x^2+\eta^2)}.$$

The adaptation therefore evaluates

$$\rho_i(E)=\frac{1}{N_s^2}\sum_{k}\sum_n w_{ni}(k)\frac{\eta}{\pi[(E-E_n(k))^2+\eta^2]}.$$

`Ns` is the number of samples along each momentum direction, so there are `Ns^2` matrices. `eta` is positive and has energy units. The weights and momentum average are dimensionless; LDOS has inverse-energy units. The output array has one row per energy and one column per selected site.

**Normalization check:** integrating one broadened peak gives

$$\int_{-\infty}^{\infty}\frac{\eta\,dE}{\pi[(E-E_n)^2+\eta^2]}=\frac1\pi\left[\arctan\frac{E-E_n}{\eta}\right]_{-\infty}^{\infty}=1.$$

For orthonormal complete eigenvectors, each electron component has total weight one across all states; the two-spin average thus gives $\sum_n w_{ni}=1$. Averaging `Ns^2` such spectra leaves $\int\rho_i(E)dE=1$. A finite plotted energy window need not contain that entire weight. The spin-summed LDOS would be twice this result and integrate to two.

As $\eta\to0^+$, peaks narrow toward delta functions; on a finite grid they require finer sampling. Far from all eigenvalues, the tails decay as $1/E^2$. `Ns=1` gives the zero-twist result. Increasing `Ns` is a sampling choice, not a normalization change.

The original script instead sums its first electron and first hole component over the full eigensystem with no factor of one half; its full-energy integral is two for a complete orthonormal basis. Its curves must not be assumed numerically identical to this explicitly defined spin average. The new implementation is serial for transparent execution; it does not close or launch the user's parallel kernels.'''), block('ldos'), md(r'''## 8. Plot energy curves, spectra and vortex-core views

`LDOSCurves` supplies explicit `{energy, value}` pairs for every site. It draws green, blue and orange curves in the order supplied in `Sites`. Reference energies are optional inputs; red lines are not automatically measured gaps.

`PairingViews` uses the **same** `Abs[delta]` values for all three views:

- **Line:** a cut along x at the core's y coordinate.
- **Map:** the magnitude at every `(x,y)` coordinate, with a colour legend.
- **Surface:** the same magnitude as height, with a conventional increasing vertical axis.

These are pairing-amplitude views, not LDOS, phase, or magnetic-field maps. Their plotted heights have energy units. Their appearance need not match an archived surface with a reversed vertical axis or a different sign convention. No curve is fitted to an image.

`SpectrumView` plots sorted eigenvalues against sorted-state index. `LDOSMap` instead maps the LDOS at one selected energy across all sites; it requires `Sites -> All`. A large LDOS value at the core is a spectral-weight observation, not proof of a particular bound-state identity.'''), block('plots'), md(r'''## 9. Export and run an explicitly labelled case

The CSV has three columns: energy, site and spin-averaged LDOS. Keeping the site column prevents interleaved curves from being treated as one line. A WXF file retains the matrices, resulting pairing, grid, settings and metadata. A new output directory is required so an earlier result is not overwritten.'''), block('export'), md(r'''### Main options

| Option | Default | Meaning |
|---|---|---|
| `Nx`, `Ny` | 20, 20 | Dimensions expected in the imported data |
| `Ns` | 10 | Samples per momentum direction |
| `Wmax`, `Intervals` | 0.10, 400 | 401 energies from -0.10 to +0.10 |
| `Eta` | 0.005 | Positive broadening in model energy units |
| `Sites` | `{190,180,1}` | Explicit example observation sites |
| `ImpurityKind`, `Strength` | `"none"`, 0 | Additional impurity, if requested |
| `CoreSite` | 190 | Core/impurity site on the chosen lattice |
| `PairingMode` | `"Imported"` | Retain the field or use `"Relax"` |
| `SeedVortex` | `False` | Retain imported field or reconstruct a seed |
| `PhaseMode` | `"Winding"` | Quadrant-aware phase; original formula remains selectable |

The default does not assume the bundled WDX data contain a specific historical vortex or impurity. The matrix checks must pass before calculation. Full $20\times20$ cases involve dense $1600\times1600$ eigensystems; `Ns=20` means 400 such decompositions per LDOS case. The example file never starts a sweep just by loading it.'''), block('run'), md(r'''## 10. Match the types of saved figures without inventing their inputs

| Saved figure family | Relevant code | Available setting / limitation |
|---|---|---|
| Energy-axis LDOS, including the two inner peaks | `CalculateLDOS`, `LDOSCurves` | The vortex curve is labelled with a 20 × 20 lattice and 20 supercells; other historical curves use 28 × 28 and 40. Each needs compatible input matrices. |
| Green/blue/orange impurity LDOS | `AddImpurity`, `CalculateLDOS`, `LDOSCurves` | Scalar and magnetic strengths 0.5, 5 and 50; sites 190 and 180 are identified; exact farther site is unverified. |
| Three pairing line profiles | `RelaxPairing`, `PairingViews[...]["Line"]` | Displayed strengths are 0.5, 1 and 5. Newly relaxed fields must converge separately. |
| Impurity energy spectra | `BdGMatrix`, `SpectrumView` | Diagonalize each case's own final matrix; choose reference lines explicitly. |
| Two-dimensional vortex-core views | `PairingViews[...]["Map"]` | Use the final onsite pairing magnitude and explicit coordinates. |
| Three-dimensional vortex-core views | `PairingViews[...]["Surface"]` | Same pairing data as the map; new plots label increasing height normally. |
| Optional spatial LDOS at a chosen energy | `LDOSMap` | This is a distinct observable, not a relabelling of a pairing map. |

**Axis ambiguity in the archived impurity curves:** the images say “Number of Sites”, but the supplied `ListPlot[Table[ldos[[n,184]],{n,1,nw}]]` varies the energy-row index at a fixed site. Its horizontal index therefore represents energy samples. This makes an energy-sample interpretation plausible for the old curves, but does not establish that this exact command produced every image. Preserve the archived labels and explain the ambiguity; use explicit energy coordinates for new curves.

The code below is linked to the diagram families by what it computes. Exact historical reproduction remains unverified until the corresponding converged matrices, observation sites, energy units and run settings are available. No new numerical curves are embedded as if they were old results.''')]
    return cells
