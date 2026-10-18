"""Build readable, self-contained notebooks from the reviewed Python functions.

Rebuilding clears saved output. Run execute_notebooks.py to regenerate output.
"""
from pathlib import Path
import inspect
import sys
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))
import bdg_models as model

THESIS = 'https://drive.google.com/file/d/1d9WheYwEY23cPAlLxxVOUfQDr9jLdqb9/view'
SOURCE = 'https://github.com/adicos0/BdG-code-for-Topological-superconductor-with-vortex-/blob/d88d12cded1a954d53c2d36d515525c8100ded1c/'
PAPER = 'https://arxiv.org/html/2012.10588v2'


def md(text):
    return nbf.v4.new_markdown_cell(text.strip())


def code(text):
    return nbf.v4.new_code_cell(text.strip())


def functions(*names):
    return code('\n\n'.join(inspect.getsource(getattr(model, name)).rstrip() for name in names))


SETUP = '''import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.special import expit
plt.rcParams.update({"figure.figsize": (7, 4), "font.size": 11})'''


def intro(number, title, scope):
    return [md(fr'''# {number}. {title}

{scope}

**Read this notebook from top to bottom.** Each section explains one calculation before showing its code. All functions needed to run it are visible below; no earlier notebook is required.

**Source labels.** **Thesis** identifies the stated model in Aditya Banerjee's [2023 MSc thesis]({THESIS}). **Source translation** identifies the active Mathematica calculation at commit `d88d12c`. **Extension / algebra** identifies new Python implementation or explicit derivations. Numerical output is evidence only for the displayed parameters. No speculative physical identification is assumed.

**Units and conventions.** Choose one energy unit $E_0$ and lattice spacing $a$. Displayed energies mean $E/E_0$, coordinates mean $x/a$, and $k_x,k_y$ mean momentum times $a$. Thus $m,\lambda,\mu,\Delta,V,k_BT$ all have energy units before rescaling. Eigenvectors have norm one and dimensionless components. This is a nonrelativistic lattice calculation; no spacetime metric signature is used. Boundary conditions and the spin/Nambu basis are stated where each matrix is built.

**Code provenance.** The original files are [BdG topological.nb]({SOURCE}BdG%20topological.nb) and [added vortex.nb]({SOURCE}added%20vortex.nb). See [SOURCE_MAP.md](../SOURCE_MAP.md) for exact cells, differences and references. Python functions are also collected in `python/bdg_models.py` for testing; the definitions here are identical.

## 1. Imports'''), code(SETUP)]


def save(filename, cells):
    notebook = nbf.v4.new_notebook(cells=cells)
    notebook.metadata = {'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
                         'language_info': {'name': 'python', 'version': '3.12'}}
    # Deterministic cell IDs make future source diffs easier to review.
    for index, cell in enumerate(notebook.cells):
        cell.id = f'cell-{index:03d}'
    nbf.validate(notebook)
    nbf.write(notebook, ROOT / 'notebooks' / filename)


cells = intro('01', 'Lattice and normal-state Hamiltonian',
              'Start with the site numbering and hopping in your code, then compare it with the Hamiltonian in thesis §4.4. The two constructions are kept separate.')
cells += [md(r'''## 2. Number the sites — source translation

The source loops over $x$ first and $y$ second. With one-based coordinates it uses
$$i=(x-1)N_y+y.$$
Python uses zero-based coordinates, so
$$i=xN_y+y,\qquad N=N_xN_y.$$
Here $N_x,N_y$ count sites along each direction and $N$ is the total number of sites. The source uses $20\times20$ with periodic boundaries. We use $6\times6$ for a quick demonstration. This changes the finite-size spectrum; it is not a reproduction of a thesis figure.'''),
          functions('coordinates'), code('''nx = ny = 6
xy = coordinates(nx, ny)
n = nx*ny
print("Number of sites:", n)
print("First six (x,y) coordinates:", xy[:6].tolist())
fig, ax = plt.subplots()
ax.scatter(xy[:, 0], xy[:, 1], color="navy")
for i, (x, y) in enumerate(xy):
    ax.annotate(str(i), (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8)
ax.set(xlabel="x/a", ylabel="y/a", title="Python site indices", aspect="equal",
       xlim=(-0.5, nx-0.2), ylim=(-0.5, ny-0.2))
plt.show()'''),
          md(r'''## 3. Build the hopping used by the active source

**Source translation:** first `Code` cell, followed by the active iteration cell. Its scalar hopping values are $3m/2$ onsite, $-m$ for four axial neighbours and $m/4$ for four diagonal neighbours. Its spin-mixing amplitude is $\lambda(1-i)$ for every axial direction.

Let $H_1$ denote the scalar hopping matrix and $S$ the spin-mixing matrix. The active electron block is
$$h_{\rm src}=\begin{pmatrix}H_1+\mu I&S\\S^*&H_1+\mu I\end{pmatrix}.$$
$S$ is symmetric for these dictionaries, so $S^*=S^\dagger$. The separately defined `Hop1x` function has a minus sign on the lower $H_1$, but **the active loop does not use that function**. We preserve the active loop.

**Algebra check:** periodic Fourier transformation of these dictionaries gives
$$H_1(k)=\frac{3m}{2}-2m(\cos k_x+\cos k_y)+m\cos k_x\cos k_y,$$
$$S(k)=2\lambda(1-i)(\cos k_x+\cos k_y).$$
Each pair of opposite axial bonds supplies $e^{ik}+e^{-ik}=2\cos k$. The four diagonal bonds supply
$(e^{ik_x}+e^{-ik_x})(e^{ik_y}+e^{-ik_y})=4\cos k_x\cos k_y$.
All terms have energy units. These expressions differ from thesis Eq. (4.8).'''),
          functions('legacy_blocks'), code('''h1, soc = legacy_blocks(nx, ny)
mu = 0.0
h_source = np.block([[h1 + mu*np.eye(n), soc],
                     [soc.conj(), h1 + mu*np.eye(n)]])
print("Electron matrix shape:", h_source.shape)
print("Hermiticity residual:", np.max(abs(h_source-h_source.conj().T)))'''),
          md(r'''## 4. Build the thesis Hamiltonian — explicitly separate extension

**Thesis:** §4.4, Eqs. (4.8)–(4.9), pp. 45–46:
$$h(k)=\lambda(\sigma_x\sin k_x+\sigma_y\sin k_y)+M_k\sigma_z-\mu I,$$
$$M_k=m\left[(2-\cos k_x-\cos k_y)-\frac14(2-\cos2k_x-\cos2k_y)\right].$$
The $\sigma$ matrices act on physical spin, $\lambda$ controls spin-dependent hopping, $m$ controls the lattice mass and $\mu$ is the chemical potential in this equation's sign convention.

**Algebra, with Fourier convention** $c_r=N^{-1/2}\sum_k e^{-ik\cdot r}c_k$:

1. Collect the constants: $2m-m/2=3m/2$.
2. A nearest axial hopping $-m\sigma_z/2$ plus its adjoint produces $-m\sigma_z\cos k_\alpha$.
3. A distance-two axial hopping $m\sigma_z/8$ plus its adjoint produces $(m/4)\sigma_z\cos2k_\alpha$.
4. A spin hopping $i\lambda\sigma_\alpha/2$ plus its adjoint gives
$$\frac{i\lambda}{2}\sigma_\alpha(e^{-ik_\alpha}-e^{ik_\alpha})=\lambda\sigma_\alpha\sin k_\alpha.$$
Thus each positive nearest bond has $-m\sigma_z/2+i\lambda\sigma_\alpha/2$ and each positive distance-two bond has $m\sigma_z/8$. These are **axial distance-two**, not diagonal, bonds.

**Limiting cases:** at $k=0$, $M_k=0$. Using $\cos k=1-k^2/2+k^4/24+O(k^6)$,
$$1-\cos k=k^2/2-k^4/24+O(k^6),$$
$$\tfrac14(1-\cos2k)=k^2/2-k^4/6+O(k^6),$$
so subtraction gives $k^4/8+O(k^6)$ and $M_k=m(k_x^4+k_y^4)/8+O(k^6)$.

**Source correction:** this spin mass breaks physical time-reversal symmetry, contrary to the sentence after thesis Eq. (4.8). Under spin time reversal, $\sigma_z\to-\sigma_z$ while $M_{-k}=M_k$. See the model discussion in [Pathak, Plugge & Franz, §2.1](https://arxiv.org/html/2012.10588v2).'''),
          functions('thesis_normal'), code('''h_thesis = thesis_normal(nx, ny, periodic=True)
e_source = eigh(h_source, eigvals_only=True)
e_thesis = eigh(h_thesis, eigvals_only=True)
fig, ax = plt.subplots()
ax.plot(e_source, ".", label="Active source translation")
ax.plot(e_thesis, "x", label="Thesis Eq. 4.8 extension")
ax.set(xlabel="Sorted eigenvalue index", ylabel="E/E0", title="Same lattice; different Hamiltonians")
ax.legend()
plt.show()
# Uniform k=0 spinors test the thesis limit h(0) = -mu I.
constant_up = np.r_[np.ones(n), np.zeros(n)] / np.sqrt(n)
print("Thesis k=0 residual at mu=0:", np.linalg.norm(h_thesis @ constant_up))'''),
          md('''## 5. What this establishes

**Computed fact:** the two implementations produce different normal-state spectra. The matrix comparison is a source audit, not a phase diagram. **Model assumption:** the thesis-based implementation uses its stated spin Hamiltonian. Later vortex examples explicitly use this extension. The self-consistency notebook instead translates the original active code.''')]
save('01_lattice_and_normal_hamiltonian.ipynb', cells)

cells = intro('02', 'BdG matrix and exact diagonalization',
              'Follow the four blocks of your active BdG matrix, explain its hole-basis signs, and check the eigensystem numerically.')
cells += [md(r'''## 2. State the basis before writing the matrix

**Thesis:** §4.1 introduces BdG diagonalization; §4.4, Eqs. (4.10)–(4.12), gives pairing and spin blocks. **Algebra:** use an explicit canonical basis to remove sign ambiguity:
$$\Psi=(c_{1\uparrow},\ldots,c_{N\uparrow},c_{1\downarrow},\ldots,c_{N\downarrow},
 c_{1\uparrow}^\dagger,\ldots,c_{N\uparrow}^\dagger,c_{1\downarrow}^\dagger,\ldots,c_{N\downarrow}^\dagger)^T.$$
For an electron matrix $h=h^\dagger$, onsite spin-singlet pairing has
$$A=\begin{pmatrix}0&D_\Delta\\-D_\Delta&0\end{pmatrix},\qquad D_\Delta=\operatorname{diag}(\Delta_1,\ldots,\Delta_N).$$
The minus sign follows from fermion antisymmetry: $A^T=-A$. Then
$$H_{\rm BdG}=\begin{pmatrix}h&A\\A^\dagger&-h^*\end{pmatrix},\qquad
\widehat H=\tfrac12\Psi^\dagger H_{\rm BdG}\Psi+\text{constant}.$$
The factor $1/2$ compensates for Nambu doubling. $H$ has size $4N\times4N$ and each entry has energy units.

**Particle–hole check:** exchange electron and hole blocks with $\tau_x$ and complex conjugate with $K$. For $C=\tau_x K$,
$$\tau_x H^*\tau_x=\begin{pmatrix}-h&A^T\\A^*&h^*\end{pmatrix}
=\begin{pmatrix}-h&-A\\-A^\dagger&h^*\end{pmatrix}=-H.$$
Therefore an eigenpair $(E,w)$ has partner $(-E,Cw)$. This is a BdG redundancy and alone does not prove a topological state.'''),
          functions('coordinates', 'legacy_blocks', 'legacy_bdg', 'canonical_bdg', 'canonical_checks'),
          md(r'''## 3. Translate the source's pairing signs

**Source translation:** the active source uses a symmetric off-diagonal pairing block $s\sigma_x$, where $s$ is its variable `st`. It is consistent with a signed hole basis, rather than canonical holes.

Define $D=\operatorname{diag}(I_N,-I_N)$ and $U=\operatorname{diag}(I_{2N},D)$. $U$ reverses only the down-hole sign. Set the canonical pairing $\Delta=-s$. Then
$$A=\begin{pmatrix}0&-D_s\\D_s&0\end{pmatrix},\quad
AD=\begin{pmatrix}0&D_s\\D_s&0\end{pmatrix},$$
$$-Dh^*D=\begin{pmatrix}-(H_1+\mu I)^*&S^*\\S&-(H_1+\mu I)^*\end{pmatrix}.$$
Consequently $H_{\rm src}=U H_{\rm canonical}(h,\Delta=-s)U^\dagger$.
This is an algebraic identification of the saved matrix; the original code does not explicitly document its Nambu basis. The spin-major ordering is inferred from the block assembly and observable indices. Periodic boundaries are retained here.'''),
          code('''nx = ny = 4
n = nx*ny
h1, soc = legacy_blocks(nx, ny)
mu = 0.13
st = np.full(n, -0.05, complex)
h = np.block([[h1 + mu*np.eye(n), soc], [soc.conj(), h1 + mu*np.eye(n)]])
h_source = legacy_bdg(h1, soc, mu, st)
h_canonical = canonical_bdg(h, -st)
signs = np.r_[np.ones(3*n), -np.ones(n)]
transformed = signs[:, None]*h_canonical*signs[None, :]
print("Source/basis-transformed difference:", np.max(abs(h_source-transformed)))
print("Canonical Hermiticity / particle-hole residuals:", canonical_checks(h_canonical))
assert np.allclose(h_source, transformed)'''),
          md(r'''## 4. Diagonalize and verify the eigenvectors

The source calls `Eigensystem[H]`. Here `scipy.linalg.eigh` solves the same Hermitian eigenproblem:
$$H w_n=E_n w_n,\qquad w_n^\dagger w_m=\delta_{nm}.$$
Mathematica returns eigenvectors as a list of rows; SciPy stores them as columns. Therefore column `vectors[:, n]` is state $n$. SciPy sorts energies in increasing order. Full-spectrum sums are unchanged by this ordering.

We check three distinct statements: the eigen-equation $HW=WE$, normalization $W^\dagger W=I$, and paired energies $E_n=-E_{4N-1-n}$. These tests assess the numerical solution, not the physical suitability of the model.'''),
          code('''energies, vectors = eigh(h_source)
residual = np.linalg.norm(h_source @ vectors-vectors*energies) / max(1, np.linalg.norm(h_source))
normalization = np.max(abs(vectors.conj().T @ vectors-np.eye(4*n)))
pairing_error = np.max(abs(energies+energies[::-1]))
print(f"Relative eigen-equation residual: {residual:.3e}")
print(f"Orthonormality residual: {normalization:.3e}")
print(f"Energy-pair residual [E0]: {pairing_error:.3e}")
fig, ax = plt.subplots()
ax.plot(np.arange(len(energies)), energies, ".", color="navy")
ax.axhline(0, color="gray", linewidth=0.7)
ax.set(xlabel="Sorted state index", ylabel="E/E0", title="Source BdG matrix: complete spectrum")
plt.show()'''),
          md(r'''## 5. Limiting case: an isolated paired site

With no hopping and $h=-\mu I_2$, a $2\times2$ sub-block is
$$B=\begin{pmatrix}-\mu&\Delta\\\Delta^*&\mu\end{pmatrix}.$$
Its characteristic equation is
$$\det(B-EI)=(-\mu-E)(\mu-E)-|\Delta|^2=E^2-\mu^2-|\Delta|^2=0,$$
so $E=\pm\sqrt{\mu^2+|\Delta|^2}$, each twice in the full spinful matrix. Both terms under the square root have units of energy squared. As $\Delta\to0$, the energies become $\pm|\mu|$; at $\mu=0$, they become $\pm|\Delta|$.'''),
          code('''mu_atomic, delta_atomic = 0.3, 0.2
atomic = canonical_bdg(-mu_atomic*np.eye(2), np.array([delta_atomic]))
expected = np.sqrt(mu_atomic**2+delta_atomic**2)*np.array([-1, -1, 1, 1])
print("Atomic energies:", eigh(atomic, eigvals_only=True))
assert np.allclose(eigh(atomic, eigvals_only=True), expected)''')]
save('02_bdg_matrix_and_diagonalization.ipynb', cells)

cells = intro('03', 'Self-consistency in the original code',
              'Translate the original pairing and density iteration, one step at a time. Record residuals so that a fixed iteration count is not mistaken for convergence.')
cells += [md(r'''## 2. Occupation numbers and source conventions

**Thesis:** §4.1, pp. 30–33, describes iterating the BdG spectrum and mean fields. **Source translation:** the equations below are the actual sums in input cell 5, expressed in the signed-hole basis explained in notebook 02. They are not a claim that every thesis mean-field formula has the same normalization.

Write each normalized eigenvector as $(u_{i\uparrow,n},u_{i\downarrow,n},v_{i\uparrow,n},v_{i\downarrow,n})$ across all sites $i$ and states $n$. The Fermi function is
$$f(E)=\frac{1}{e^{E/(k_BT)}+1}.$$
$k_BT$ is a thermal **energy**, not temperature alone. The source sets $k_BT=0.001$ and replaces $f$ by exactly $0$ above $0.02$ and $1$ below $-0.02$, in its chosen energy unit. We keep that cutoff. It differs slightly from the smooth Fermi function. As $T\to0$, $f(E<0)\to1$, $f(E>0)\to0$, and $f(0)=1/2$.

Every sum below is over **all $4N$ states**:
$$n_{i\uparrow}=\sum_n|u_{i\uparrow,n}|^2 f(E_n),\quad
n_{i\downarrow}=\sum_n|v_{i\downarrow,n}|^2[1-f(E_n)],$$
$$n_i=n_{i\uparrow}+n_{i\downarrow},\qquad m_i=n_{i\uparrow}-n_{i\downarrow},$$
$$s_i^{\rm proposed}=v_s\sum_n\left(u_{i\uparrow,n}v_{i\downarrow,n}^*
+u_{i\downarrow,n}v_{i\uparrow,n}^*\right)f(E_n).$$
Here $n_i$ is electrons per site ($0\le n_i\le2$), $m_i$ is a dimensionless spin-population difference (not a magnetic moment in SI units), and $v_s$ is the source pairing coupling. Eigenvector products and $f$ are dimensionless, so $s_i$ has the same energy units as $v_s$. Do not add a positive-spectrum factor of two to these full-spectrum formulas.'''),
          functions('coordinates', 'legacy_blocks', 'legacy_bdg', 'source_fermi', 'source_observables'),
          md(r'''## 3. Update the pairing and chemical potential

**Source translation:** initialize $s_i=-0.01$, $\mu=0$. At each iteration, compute the mean density $\bar n=N^{-1}\sum_i n_i$ and apply
$$\mu_{\rm new}=\mu-\alpha(n_{\rm target}-\bar n),\qquad \alpha=0.1E_0,$$
$$s_{i,\rm new}=\tfrac12s_i^{\rm proposed}+\tfrac12s_i.$$
The factor $1/2$ mixes old and proposed gaps. In the source electron block the chemical potential enters as **$+\mu$**. If the density is too small, $n_{\rm target}-\bar n>0$, so the update lowers $\mu$, lowers electron energies and favours filling. The thesis convention $-\mu$ would require a separate sign analysis.

**Added diagnostic:** $r_s=\max_i|s_i^{\rm proposed}-s_i|$ has energy units and $r_n=|\bar n-n_{\rm target}|$ is dimensionless. We report them independently. Small $r_s$ alone is insufficient. The original loop performs exactly 100 iterations, with no convergence test; this notebook preserves that behavior and adds residual reporting, including a fresh evaluation after the last update.'''),
          functions('source_iteration'),
          code('''# Small periodic lattice for a quick, fully executed tutorial.
# Original source uses nx=ny=20 (1600 x 1600 BdG matrix).
nx = ny = 4
cases = {"BdG topological.nb": (-1.45, 1.2),
         "added vortex.nb": (-1.25, 1.0)}
runs = {}
for label, (vs, target) in cases.items():
    mu, st, history, residuals = source_iteration(nx, ny, vs, target, steps=100)
    runs[label] = (mu, st, history, residuals)
    converged = residuals[0] < 1e-6 and residuals[1] < 1e-6
    print(label)
    print(f"  mu/E0={mu:.6f}; mean |st|/E0={abs(st).mean():.6g}")
    print(f"  final gap residual/E0={residuals[0]:.3e}; density residual={residuals[1]:.3e}")
    print("  Meets both displayed 1e-6 tolerances:", converged)'''),
          md('''## 4. Inspect convergence rather than assume it

The curves below record fields **before each update**. The printed residuals above are recomputed for the final returned fields. A plateau or small gap does not by itself establish convergence or superconductivity. Changing lattice size changes the finite spectrum and can change the iteration behavior.'''),
          code('''fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
for label, (_, _, history, _) in runs.items():
    axes[0].plot(history[:, 0], history[:, 2], label=label)
    axes[1].semilogy(history[:, 0], np.maximum(history[:, 4], 1e-16), label=label)
    axes[2].semilogy(history[:, 0], np.maximum(history[:, 5], 1e-16), label=label)
axes[0].set(ylabel="Mean density [electrons/site]", title="Density")
axes[1].set(ylabel="Gap residual / E0", title="Pairing update")
axes[2].set(ylabel="Density residual [electrons/site]", title="Target density error")
for ax in axes:
    ax.set_xlabel("Iteration")
axes[0].legend(fontsize=8)
fig.tight_layout()
plt.show()'''),
          md('''## 5. Scope of this result

**Verified calculation:** this is a Python translation of the active periodic, spatially uniform initialization and update rules. Dense Hermitian diagonalization replaces Mathematica `Eigensystem`; the saved Mathematica outputs have not been independently reproduced with Wolfram. The prototype vortex block is commented out in both source notebooks. Therefore these iterations do not create a vortex. The next notebook introduces a separately labelled imposed vortex field from the thesis; it is not the converged output of this notebook.''')]
save('03_self_consistent_pairing.ipynb', cells)

cells = intro('04', 'Vortex profile and finite-sample spectrum',
              'Implement the thesis vortex ansatz with a correct phase winding, then diagonalize a finite open sample. This is a thesis-based extension of the commented vortex prototype.')
cells += [md(r'''## 2. Define the vortex before inserting numbers

**Thesis:** §4.2, Eq. (4.2), uses
$$\Delta(x,y)=\Delta_0\tanh(r/\xi)e^{iq\phi},\qquad
r=\sqrt{(x-x_0)^2+(y-y_0)^2}.$$
Here $\Delta_0$ is the far-field pairing energy, $(x_0,y_0)$ is the vortex centre, $\xi>0$ is the recovery length, and $q=1$ describes one positive winding. $r/\xi$ and $q\phi$ are dimensionless.

**Explicit correction to thesis Eq. (4.4) and commented code:** use
$$\phi=\operatorname{atan2}(y-y_0,x-x_0).$$
The printed/source expression $\arctan(|y-y_0|/|x-x_0|)$ loses quadrant information and cannot describe a full $2\pi$ winding. The corrected phase is essential, not a cosmetic plotting change. At $r=0$, phase is undefined but $\Delta=0$, so its assigned numerical value has no physical effect.

**Limits:** since $\tanh z=z+O(z^3)$, near the centre $|\Delta|\simeq\Delta_0r/\xi$. As $r/\xi\to\infty$, $|\Delta|\to\Delta_0$. A closed loop around the centre has phase increment $2\pi q$.

**Length convention:** thesis Eq. (4.3) writes a coherence-length relation involving a penetration length and gap, but without a specified unit conversion it cannot be used as an SI length formula. Here $\xi/a$ is an independently chosen dimensionless input, not a material prediction.

**Boundary assumption:** use open boundaries for this isolated imposed vortex. No vector potential, magnetic screening, flux quantization calculation, or self-consistent relaxation is included. A periodic vortex lattice needs additional gauge construction; see [Pathak, Plugge & Franz, §2.2](https://arxiv.org/html/2012.10588v2).'''),
          functions('coordinates', 'thesis_normal', 'canonical_bdg', 'canonical_checks', 'vortex_field'),
          code('''# Illustrative finite-size parameters, NOT a thesis figure reproduction.
# m and lambda follow the original scale; delta0 and xi are tutorial inputs.
nx = ny = 11
n = nx*ny
center = (5, 5)
delta0, xi, mu = 0.2, 1.5, -0.7
delta = vortex_field(nx, ny, delta0, xi, center)
xy = coordinates(nx, ny)
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, values, cmap, limits, title in [
    (axes[0], abs(delta), "viridis", (0, delta0), "Imposed pairing amplitude"),
    (axes[1], np.ma.masked_where(abs(delta)<1e-12, np.angle(delta)), "twilight", (-np.pi, np.pi), "Pairing phase; core masked")]:
    plot = ax.imshow(values.reshape(nx, ny).T, origin="lower", interpolation="nearest",
                     cmap=cmap, vmin=limits[0], vmax=limits[1])
    ax.set(xlabel="x/a", ylabel="y/a", title=title)
    fig.colorbar(plot, ax=ax, label="|Delta|/E0" if ax is axes[0] else "Phase [rad]")
fig.tight_layout()
plt.show()'''),
          md(r'''## 3. Check the winding and amplitude limits

Take a counterclockwise square loop around the core, without passing through it. For consecutive complex gap values $\Delta_j$, the phase increment is $\arg(\Delta_{j+1}\Delta_j^*)$. Sum these increments and divide by $2\pi$:
$$q_{\rm measured}=\frac{1}{2\pi}\sum_j\arg(\Delta_{j+1}\Delta_j^*).$$
The loop must be sampled finely enough that each true phase step is smaller than $\pi$; the loop below satisfies this. We also compare the core and far-corner amplitudes with the stated ansatz.'''),
          code('''loop_xy = ([(x, 2) for x in range(2, 8)] + [(8, y) for y in range(2, 8)]
           + [(x, 8) for x in range(8, 2, -1)] + [(2, y) for y in range(8, 2, -1)])
loop = np.array([delta[x*ny+y] for x, y in loop_xy])
winding = np.angle(np.roll(loop, -1)*loop.conj()).sum()/(2*np.pi)
core = center[0]*ny+center[1]
print("Measured counterclockwise winding:", winding)
print("Core amplitude/E0:", abs(delta[core]))
print("Corner amplitude / delta0:", abs(delta[0])/delta0)
assert np.isclose(winding, 1)
assert abs(delta[core]) < 1e-12'''),
          md(r'''## 4. Diagonalize the imposed-vortex Hamiltonian

**Extension:** use the thesis normal-state Hamiltonian from notebook 01, the canonical spin-singlet BdG matrix from notebook 02, and the imposed profile above. The pairing is held fixed. We compare against the same open sample with uniform $\Delta_0$, so the only changed input in this comparison is the gap profile.

For a state $w_n$, the total Nambu probability at site $i$ is
$$p_{i,n}=|u_{i\uparrow,n}|^2+|u_{i\downarrow,n}|^2+|v_{i\uparrow,n}|^2+|v_{i\downarrow,n}|^2.$$
Its site sum is one. We average over the four states closest to zero energy; this makes the plot insensitive to rotations within a fully included degenerate subspace. It is a probability map, **not LDOS**. The selected group must not cut through a degeneracy; we report the separation to the next state.'''),
          code('''electron = thesis_normal(nx, ny, mu=mu, periodic=False)
h_vortex = canonical_bdg(electron, delta)
h_uniform = canonical_bdg(electron, np.full(n, delta0))
energies, vectors = eigh(h_vortex)
uniform_energies = eigh(h_uniform, eigvals_only=True)
print("Hermiticity / particle-hole residuals:", canonical_checks(h_vortex))
order = np.argsort(abs(energies))
selected = order[:4]
separation = abs(energies[order[4]])-abs(energies[order[3]])
print("Selected energies/E0:", energies[selected])
print("Separation from next |E| / E0:", separation)
probability = (abs(vectors[:, selected])**2).reshape(4, n, 4).sum(axis=0).mean(axis=1)
print("Spatial probability sum:", probability.sum())
boundary = (xy[:, 0] == 0) | (xy[:, 0] == nx-1) | (xy[:, 1] == 0) | (xy[:, 1] == ny-1)
core_region = np.linalg.norm(xy-np.array(center), axis=1) <= 2
print("Mean selected-state weight on outermost boundary:", probability[boundary].sum())
print("Mean selected-state weight within radius 2a of core:", probability[core_region].sum())
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for spectrum, label, marker in [(uniform_energies, "Uniform pairing", "."), (energies, "Imposed vortex", "x")]:
    middle = spectrum[len(spectrum)//2-12:len(spectrum)//2+12]
    axes[0].plot(np.arange(-12, 12), middle, marker, label=label)
axes[0].set(xlabel="State index relative to spectrum centre", ylabel="E/E0", title="24 central eigenvalues")
axes[0].legend()
plot = axes[1].imshow(probability.reshape(nx, ny).T, origin="lower", interpolation="nearest", cmap="magma")
axes[1].set(xlabel="x/a", ylabel="y/a", title="Mean probability of 4 closest states")
fig.colorbar(plot, ax=axes[1], label="Probability per site")
fig.tight_layout()
plt.show()'''),
          md('''## 5. Interpretation and limits

**Computed facts:** the output gives finite-sample energies, winding and spatial weights for an imposed profile. With the displayed parameters, the selected four-state average has most of its weight near the outer boundary, rather than the core. The printed boundary and core-region weights quantify this observation; they are not a topological classification. **Not established here:** self-consistency of the vortex, a thermodynamic gap, or a Majorana zero mode. A small eigenvalue alone does not settle these questions. The thesis discusses such interpretations, while these notebooks keep the calculation and identification separate. [Fu & Kane (2008)](https://arxiv.org/abs/0707.1692) provides the underlying theoretical context for induced superconductivity and vortices on a topological-insulator surface.''')]
save('04_vortex_profile_and_spectrum.ipynb', cells)

cells = intro('05', 'Impurities and local density of states',
              'Add the impurity terms written in the thesis to the same imposed-vortex model. Compute a properly normalized electron LDOS and compare scalar and magnetic potentials.')
cells += [md(r'''## 2. Add one impurity — thesis-based extension

**Thesis:** §4.3, Eqs. (4.6)–(4.7), and §4.4.1–4.4.2 distinguish
$$H_{\rm imp}^{\rm scalar}=V(n_{j\uparrow}+n_{j\downarrow}),\qquad
H_{\rm imp}^{\rm magnetic}=V(n_{j\uparrow}-n_{j\downarrow}).$$
$j$ is the impurity site, $V$ is an energy, and $n_{j\sigma}=c_{j\sigma}^\dagger c_{j\sigma}$ is an occupation operator. A scalar impurity adds $(V,V)$ to the two spin diagonals. A magnetic impurity polarized along $z$ adds $(V,-V)$.

Construct the electron matrix first, then rebuild the complete BdG matrix. This automatically puts the conjugated, sign-reversed impurity into the hole sector. Adding the same sign to both particle and hole blocks would break this convention.

**Scope:** impurity terms are disabled in the committed Mathematica loop. The Python calculations below implement the thesis equations as an extension. The pairing profile and chemical potential are held fixed at every impurity strength. Thus this is a **frozen-profile comparison**, not a self-consistent reproduction of the thesis impurity figures. Open boundaries and the canonical basis from notebook 02 are used.'''),
          functions('coordinates', 'thesis_normal', 'canonical_bdg', 'canonical_checks', 'vortex_field', 'add_impurity', 'electron_ldos'),
          md(r'''## 3. Define the electron LDOS without double counting

**Thesis:** §4.1, p. 32, introduces local density of states. **Explicit convention and algebra here:** use all $4N$ BdG states and the electron components only:
$$\rho_i(E)=\sum_{n=1}^{4N}\left(|u_{i\uparrow,n}|^2+|u_{i\downarrow,n}|^2\right)L_\eta(E-E_n),$$
$$L_\eta(x)=\frac{\eta}{\pi(x^2+\eta^2)},\qquad \eta>0.$$
$\eta$ is a numerical broadening energy. $L_\eta$ and $\rho_i$ have inverse-energy units. This spin-summed convention differs by a factor of two from a spin-averaged LDOS. Do not also add hole weights when using this complete signed spectrum.

The normalization follows in steps:
$$\int_{-\infty}^{\infty}L_\eta(x)\,dx=
\left[\frac1\pi\arctan(x/\eta)\right]_{-\infty}^{\infty}=1,$$
$$\int\rho_i(E)\,dE=\sum_n|u_{i\uparrow,n}|^2+\sum_n|u_{i\downarrow,n}|^2=1+1=2.$$
The last equality follows from completeness of the normalized eigenvectors. The finite plotting window does not contain all spectral weight; the code checks the exact weight sum instead. As $\eta\to0^+$, each Lorentzian approaches a delta distribution. A finite-width zero-energy peak does not prove an exact zero eigenvalue.'''),
          code('''nx = ny = 11
n = nx*ny
center = (5, 5)
site = center[0]*ny+center[1]
delta0, xi, mu = 0.2, 1.5, -0.7
electron = thesis_normal(nx, ny, mu=mu, periodic=False)
delta = vortex_field(nx, ny, delta0, xi, center)
grid = np.linspace(-0.6, 0.6, 601)
eta = 0.02
strengths = [0.0, 0.5, 2.0]  # tutorial sweep in E0; not fitted thesis results
results = {}
for magnetic in [False, True]:
    label = "Magnetic (z)" if magnetic else "Scalar"
    for strength in strengths:
        changed = add_impurity(electron, site, strength, magnetic)
        h = canonical_bdg(changed, delta)
        assert max(canonical_checks(h)) < 1e-12
        energies, vectors = eigh(h)
        rho, weight = electron_ldos(energies, vectors, site, grid, eta)
        assert np.isclose(weight, 2)
        results[label, strength] = (energies, rho)
        print(f"{label:12s} V/E0={strength:3.1f}: min |E|/E0={min(abs(energies)):.6f}; total electron weight={weight:.6f}")'''),
          md('''## 4. Compare spectra and core LDOS

The impurity sits at the imposed vortex centre. All panels use the same energy window and broadening, making the curves comparable. The zero-potential results must agree because the scalar and magnetic perturbations both vanish at $V=0$.'''),
          code('''assert np.allclose(results["Scalar", 0.0][0], results["Magnetic (z)", 0.0][0])
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for ax, label in zip(axes, ["Scalar", "Magnetic (z)"]):
    for strength in strengths:
        ax.plot(grid, results[label, strength][1], label=f"V/E0 = {strength:g}")
    ax.set(xlabel="E/E0", title=f"{label}: core electron LDOS")
    ax.legend()
axes[0].set_ylabel("LDOS [1/E0], spin summed")
fig.tight_layout()
plt.show()
fig, ax = plt.subplots()
for label, marker in [("Scalar", "o"), ("Magnetic (z)", "s")]:
    ax.plot(strengths, [min(abs(results[label, v][0])) for v in strengths], marker+"--", label=label)
ax.set(xlabel="Impurity strength V/E0", ylabel="Smallest |E| / E0", title="Three sampled strengths; lines guide the eye")
ax.legend()
plt.show()'''),
          md('''## 5. What can be concluded

**Computed facts:** the tables and curves show the response of this specified finite, open, fixed-profile model. The smallest-absolute-energy curves overlap to the displayed precision for this particular sweep, while the LDOS curves differ. Equality of one spectral statistic does not make the two perturbations equivalent. **Model assumptions:** onsite scalar or z-polarized magnetic perturbation; no Hartree correction; fixed chemical potential, gap amplitude and vortex phase; phenomenological Lorentzian broadening.

**Interpretation not established:** neither a zero-bias LDOS feature nor weak movement of a finite-size eigenvalue is sufficient to identify a protected Majorana mode. Such an identification requires additional spatial, particle–hole, finite-size and topology checks. This notebook does not claim that the tutorial sweep reproduces the thesis figures or experimentally validates the model.

For sources and the precise distinction between original code, algebraic conventions and extensions, return to [SOURCE_MAP.md](../SOURCE_MAP.md).''')]
save('05_impurities_and_ldos.ipynb', cells)
print('Built five notebooks.')
