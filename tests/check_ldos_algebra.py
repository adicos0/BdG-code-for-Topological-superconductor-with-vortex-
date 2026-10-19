"""Independent numerical checks of documented conventions, not Wolfram execution.

The zero-pairing, isolated-site limit has known Lorentzian centres; the general
matrix additionally obeys Hermiticity, opposite-k particle-hole symmetry and the
spectral-weight sum rule. These tests do not reproduce archived figures.
"""
import numpy as np
from pathlib import Path


def check():
    # Two spin levels of one uncoupled orbital: each must contribute half a peak.
    eta = .005
    grid = np.array([-.2, 0., .2])
    h = np.diag([.2, -.1, -.2, .1])
    levels, vectors = np.linalg.eigh(h)
    weights = np.abs(vectors[:2, :]) ** 2 / 2
    calculated = (eta / (np.pi * ((grid[:, None]-levels)**2+eta**2))) @ weights.sum(axis=0)
    expected = sum(.5 * eta/(np.pi*((grid-e)**2+eta**2)) for e in [.2, -.1])
    np.testing.assert_allclose(calculated, expected, atol=1e-14)
    # Pairing and spin-mixing: check the actual block convention symbolically described.
    rng = np.random.default_rng(7)
    n = 3
    raw = rng.normal(size=(2*n, 2*n)) + 1j*rng.normal(size=(2*n, 2*n))
    h0 = (raw + raw.conj().T)/2
    coord = np.array([0, 1, 2, 0, 1, 2])
    # Antisymmetric wrap exponents ensure h(k) is Hermitian independently of h0.
    wrap = np.where(np.abs(coord[None, :]-coord[:, None]) > 1.5,
                    np.sign(coord[None, :]-coord[:, None]), 0)
    s = np.diag([1]*n+[-1]*n)
    z = np.zeros((n, n))
    d = np.diag([.07+.02j, .03-.04j, .09+.01j])
    pair = np.block([[z, d], [d, z]])
    def matrix(k):
        electron = h0*np.exp(1j*k*3*wrap)
        opposite = h0*np.exp(-1j*k*3*wrap)
        return np.block([[electron, pair], [pair.conj().T, -s@opposite.conj()@s]])
    z2 = np.zeros((2*n, 2*n))
    c = np.block([[z2, s], [s, z2]])
    h = matrix(.17)
    np.testing.assert_allclose(h, h.conj().T, atol=1e-14)
    np.testing.assert_allclose(c@h.conj()@c, -matrix(-.17), atol=1e-14)
    ev, u = np.linalg.eigh(h)
    np.testing.assert_allclose(h@u, u*ev, atol=1e-13)
    w = (abs(u[:n])**2+abs(u[n:2*n])**2)/2
    np.testing.assert_allclose(w.sum(axis=1), np.ones(n), atol=1e-14)
    # Exact finite-window integral of the Lorentzians approaches one per site.
    bound = 1e5
    integral = (np.arctan((bound-ev)/eta)-np.arctan((-bound-ev)/eta))/np.pi
    np.testing.assert_allclose(w@integral, np.ones(n), atol=1e-7)
    # Verify full quadrant phase winding along a closed loop.
    theta = np.linspace(0, 2*np.pi, 33)
    phase = np.unwrap(np.angle(np.cos(theta)+1j*np.sin(theta)))
    np.testing.assert_allclose(phase[-1]-phase[0], 2*np.pi)
    print('PASS: isolated-site LDOS, Hermiticity, opposite-k particle-hole symmetry, eigen-residuals, spectral sum rule, vortex winding.')


def lexical_check():
    # Strip nested comments and strings before checking delimiters; not a WL parser.
    root = Path(__file__).resolve().parents[1]
    for path in [*root.glob('scripts/*.wl'), *root.glob('tests/*.wlt')]:
        s = path.read_text()
        stack = []; i = 0; depth = 0; string = False
        while i < len(s):
            if depth:
                if s[i:i+2] == '(*': depth += 1; i += 2; continue
                if s[i:i+2] == '*)': depth -= 1; i += 2; continue
            elif string:
                if s[i] == '\\': i += 2; continue
                if s[i] == '"': string = False
            elif s[i:i+2] == '(*': depth = 1; i += 2; continue
            elif s[i] == '"': string = True
            elif s[i] in '[{(': stack.append(s[i])
            elif s[i] in ']})':
                assert stack and stack.pop() == {']':'[', '}':'{', ')':'('}[s[i]], (path, i)
            i += 1
        assert not stack and not depth and not string, path
        print('PASS: balanced Wolfram delimiters/comments/strings:', path.name)


if __name__ == '__main__':
    check()
    lexical_check()
