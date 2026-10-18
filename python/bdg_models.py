"""Auditable Python translation and explicitly separate thesis-model extensions.

Basis: spin-major electrons, then spin-major holes. See SOURCE_MAP.md.
All energies are in one chosen unit E0; coordinates are in lattice spacings.
"""
import numpy as np
from scipy.linalg import eigh
from scipy.special import expit


def coordinates(nx, ny):
    """Zero-based version of the source's x-outer, y-inner site order."""
    return np.array([(x, y) for x in range(nx) for y in range(ny)])


def legacy_blocks(nx, ny, m=0.5, lam=1.0):
    """Translate first Code cell: periodic dictionary hopping, not thesis Eq. 4.9.

    Require >=3 sites per direction so nearest-neighbour images are distinct.
    The source selects the first numeric periodic image, rather than summing
    duplicate images. These coincide for the supported sizes.
    """
    if min(nx, ny) < 3:
        raise ValueError('Legacy translation requires nx, ny >= 3.')
    xy = coordinates(nx, ny)
    n = nx * ny
    hopping = {(0, 0): 1.5*m}
    for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        hopping[dx, dy] = -m
    for dx in [-1, 1]:
        for dy in [-1, 1]:
            hopping[dx, dy] = m/4
    mixing = {d: lam*(1-1j) for d in [(1, 0), (-1, 0), (0, 1), (0, -1)]}

    def from_dictionary(entries):
        matrix = np.zeros((n, n), complex)
        for i, (x, y) in enumerate(xy):
            for (dx, dy), value in entries.items():
                j = ((x-dx) % nx)*ny + (y-dy) % ny
                matrix[i, j] = value
        return matrix.T  # source: tij[i][j] = tijmn[[j,i]]

    return from_dictionary(hopping), from_dictionary(mixing)


def legacy_bdg(h1, soc, mu, st):
    """Translate active hh1, hh2, SCs and H of Mathematica input cell 5."""
    n = len(st)
    diagonal = h1 + mu*np.eye(n)
    electron = np.block([[diagonal, soc], [soc.conj(), diagonal]])
    hole = np.block([[-diagonal.conj(), soc.conj()],
                     [soc, -diagonal.conj()]])
    zero = np.zeros((n, n), complex)
    pairing = np.block([[zero, np.diag(st)], [np.diag(st), zero]])
    return np.block([[electron, pairing], [pairing.conj(), hole]])


def canonical_bdg(electron, delta):
    """Standard spin-singlet BdG in (c_up,c_down,c_up†,c_down†)."""
    n = len(delta)
    zero = np.zeros((n, n), complex)
    pairing = np.block([[zero, np.diag(delta)], [-np.diag(delta), zero]])
    return np.block([[electron, pairing],
                     [pairing.conj().T, -electron.conj()]])


def source_fermi(energies, kT=0.001):
    """Original cutoff at +/-0.02, with a stable exponential evaluation."""
    return np.where(energies > 0.02, 0.0,
                    np.where(energies < -0.02, 1.0, expit(-energies/kT)))


def source_observables(energies, vectors, vs):
    """Original full-spectrum sums; numpy eigenvectors are COLUMNS."""
    n = vectors.shape[0]//4
    up, down, hole_up, hole_down = np.split(vectors, 4, axis=0)
    f = source_fermi(energies)
    n_up = abs(up)**2 @ f
    n_down = abs(hole_down)**2 @ (1-f)
    st_new = vs*((up*hole_down.conj() + down*hole_up.conj()) @ f)
    return n_up+n_down, n_up-n_down, st_new


def source_iteration(nx=4, ny=4, vs=-1.45, target=1.2, steps=100):
    """Original fixed-count iteration; residual logging is added for diagnosis."""
    h1, soc = legacy_blocks(nx, ny)
    st = np.full(nx*ny, -0.01, complex)
    mu = 0.0
    history = []
    for iteration in range(steps):
        energies, vectors = eigh(legacy_bdg(h1, soc, mu, st))
        density, magnetization, proposed = source_observables(energies, vectors, vs)
        mean_density = density.mean()
        history.append([iteration+1, mu, mean_density, abs(st).mean(),
                        np.max(abs(proposed-st)), abs(mean_density-target),
                        abs(magnetization).mean()])
        mu = mu - 0.1*(target-mean_density)
        st = 0.5*proposed + 0.5*st
    # Extra diagnostic evaluation of the returned, UPDATED fields.
    energies, vectors = eigh(legacy_bdg(h1, soc, mu, st))
    density, magnetization, proposed = source_observables(energies, vectors, vs)
    residuals = np.array([np.max(abs(proposed-st)), abs(density.mean()-target)])
    return mu, st, np.array(history), residuals


def thesis_normal(nx, ny, m=0.5, lam=1.0, mu=0.0, periodic=False):
    """Extension from thesis Eqs. 4.8-4.9, not the active source hopping.

    Fourier convention c_r = sum_k exp(-i k.r)c_k / sqrt(N).
    Build site-major 2x2 blocks, then permute to spin-major order.
    """
    identity = np.eye(2)
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1., -1.])
    n = nx*ny
    h = np.kron(np.eye(n), 1.5*m*sz-mu*identity).astype(complex)
    for i, (x, y) in enumerate(coordinates(nx, ny)):
        for dx, dy, spin in [(1, 0, sx), (0, 1, sy)]:
            for distance, hop in [(1, -m/2*sz + 0.5j*lam*spin),
                                  (2, m/8*sz)]:
                xx, yy = x+distance*dx, y+distance*dy
                if not periodic and not (0 <= xx < nx and 0 <= yy < ny):
                    continue
                j = (xx % nx)*ny + (yy % ny)
                h[2*i:2*i+2, 2*j:2*j+2] += hop
                h[2*j:2*j+2, 2*i:2*i+2] += hop.conj().T
    order = np.r_[np.arange(0, 2*n, 2), np.arange(1, 2*n, 2)]
    return h[np.ix_(order, order)]


def vortex_field(nx, ny, delta0=0.2, xi=1.5, center=None, winding=1):
    """Imposed open-sample extension of Eq.4.2 with quadrant-correct phase."""
    if xi <= 0:
        raise ValueError('xi must be positive.')
    if center is None:
        center = ((nx-1)/2, (ny-1)/2)
    xy = coordinates(nx, ny)-np.asarray(center)
    radius = np.hypot(xy[:, 0], xy[:, 1])
    phase = np.arctan2(xy[:, 1], xy[:, 0])
    return delta0*np.tanh(radius/xi)*np.exp(1j*winding*phase)


def add_impurity(electron, site, strength, magnetic=False):
    """Frozen-profile extension: Eq.4.6 scalar or Eq.4.7 z-magnetic impurity."""
    h = electron.copy()
    n = h.shape[0]//2
    h[site, site] += strength
    h[n+site, n+site] += -strength if magnetic else strength
    return h


def electron_ldos(energies, vectors, site, grid, eta=0.02):
    """Spin-summed electron LDOS using ALL BdG eigenstates exactly once."""
    if eta <= 0:
        raise ValueError('eta must be positive.')
    n = vectors.shape[0]//4
    weights = abs(vectors[site])**2 + abs(vectors[n+site])**2
    lorentzian = eta/np.pi/((grid[:, None]-energies[None, :])**2+eta**2)
    return lorentzian @ weights, weights.sum()


def canonical_checks(matrix):
    """Hermiticity and particle-hole residuals (max absolute entries)."""
    half = len(matrix)//2
    exchange = np.r_[np.arange(half, 2*half), np.arange(half)]
    return (np.max(abs(matrix-matrix.conj().T)),
            np.max(abs(matrix.conj()[np.ix_(exchange, exchange)]+matrix)))
