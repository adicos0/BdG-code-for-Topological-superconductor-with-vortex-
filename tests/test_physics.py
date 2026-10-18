"""Independent limits and structural checks, not a claim of Wolfram reproduction."""
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import eigh
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'python'))
from bdg_models import (coordinates, legacy_blocks, legacy_bdg, canonical_bdg,
                        canonical_checks, thesis_normal, vortex_field,
                        add_impurity, electron_ldos, source_observables, source_fermi)


def test_legacy_dictionary_against_source_image_search():
    # Independent literal periodic-image lookup in the order of the source Ifs.
    nx, ny, m, lam = 4, 5, 0.5, 1.0
    xy = coordinates(nx, ny)
    h1, soc = legacy_blocks(nx, ny, m, lam)
    for i, (x, y) in enumerate(xy):
        for j, (xx, yy) in enumerate(xy):
            expected_h, expected_s = None, None
            for shift_x in [-nx, 0, nx]:
                for shift_y in [-ny, 0, ny]:
                    dx, dy = x-(xx+shift_x), y-(yy+shift_y)
                    if expected_h is None:
                        if (dx, dy) == (0, 0):
                            expected_h = 1.5*m
                        elif abs(dx)+abs(dy) == 1:
                            expected_h = -m
                        elif abs(dx) == abs(dy) == 1:
                            expected_h = m/4
                    if expected_s is None and abs(dx)+abs(dy) == 1:
                        expected_s = lam*(1-1j)
            assert h1[j, i] == (expected_h or 0)
            assert soc[j, i] == (expected_s or 0)


def test_thesis_fourier_blocks_and_bands():
    # Plane waves with stated exp(-ik.r) convention test signs, spin, range.
    nx, ny, m, lam, mu = 5, 6, 0.47, 1.13, -0.23
    xy = coordinates(nx, ny)
    n = nx*ny
    h = thesis_normal(nx, ny, m, lam, mu, periodic=True)
    predicted_bands = []
    for ix in range(nx):
        for iy in range(ny):
            kx, ky = 2*np.pi*ix/nx, 2*np.pi*iy/ny
            mass = m*((2-np.cos(kx)-np.cos(ky))-(2-np.cos(2*kx)-np.cos(2*ky))/4)
            mixing = lam*(np.sin(kx)-1j*np.sin(ky))
            hk = np.array([[mass-mu, mixing], [mixing.conjugate(), -mass-mu]])
            wave = np.exp(-1j*(kx*xy[:, 0]+ky*xy[:, 1]))/np.sqrt(n)
            plane = np.zeros((2*n, 2), complex)
            plane[:n, 0] = wave
            plane[n:, 1] = wave
            np.testing.assert_allclose(h@plane, plane@hk, atol=2e-14)
            predicted_bands.extend(eigh(hk, eigvals_only=True))
    np.testing.assert_allclose(eigh(h, eigvals_only=True), sorted(predicted_bands), atol=2e-14)


def test_complex_pairing_basis_change_and_particle_hole():
    rng = np.random.default_rng(17)
    n = 12
    h1, soc = legacy_blocks(3, 4)
    st = rng.normal(size=n)*0.1 + 1j*rng.normal(size=n)*0.1
    mu = 0.17
    h = np.block([[h1+mu*np.eye(n), soc], [soc.conj(), h1+mu*np.eye(n)]])
    canonical = canonical_bdg(h, -st)
    signs = np.r_[np.ones(3*n), -np.ones(n)]
    np.testing.assert_allclose(legacy_bdg(h1, soc, mu, st),
                               signs[:, None]*canonical*signs[None, :])
    assert max(canonical_checks(canonical)) < 1e-14
    energies = eigh(canonical, eigvals_only=True)
    np.testing.assert_allclose(energies, -energies[::-1], atol=1e-13)


def test_atomic_and_normal_state_limits():
    mu, delta = 0.31, 0.13+0.07j
    atomic = canonical_bdg(-mu*np.eye(2), np.array([delta]))
    expected = np.sqrt(mu**2+abs(delta)**2)*np.array([-1, -1, 1, 1])
    np.testing.assert_allclose(eigh(atomic, eigvals_only=True), expected, atol=1e-14)
    h = thesis_normal(4, 5, mu=mu, periodic=False)
    electron = eigh(h, eigvals_only=True)
    normal = canonical_bdg(h, np.zeros(20))
    np.testing.assert_allclose(eigh(normal, eigvals_only=True),
                               np.sort(np.r_[electron, -electron]), atol=1e-13)


def test_source_observables_and_density_completeness():
    h1, soc = legacy_blocks(3, 4)
    energies, vectors = eigh(legacy_bdg(h1, soc, 0.13, np.full(12, -0.07)))
    density, magnetization, _ = source_observables(energies, vectors, -1.45)
    f = source_fermi(energies)
    # Physical electron-sector density agrees with source hole-sector formula.
    electron_density = (abs(vectors[:12])**2 + abs(vectors[12:24])**2) @ f
    np.testing.assert_allclose(density, electron_density, atol=1e-13)
    assert np.all((density >= -1e-13) & (density <= 2+1e-13))
    assert np.all(abs(magnetization) <= 1+1e-13)
    np.testing.assert_allclose(source_fermi(-energies), 1-f, atol=1e-15)


def test_vortex_winding_core_and_reversal():
    nx = ny = 9
    loop_xy = ([(x, 1) for x in range(1, 7)] + [(7, y) for y in range(1, 7)]
               + [(x, 7) for x in range(7, 1, -1)] + [(1, y) for y in range(7, 1, -1)])
    for q in [-1, 1]:
        delta = vortex_field(nx, ny, center=(4, 4), winding=q)
        loop = np.array([delta[x*ny+y] for x, y in loop_xy])
        measured = np.angle(np.roll(loop, -1)*loop.conj()).sum()/(2*np.pi)
        np.testing.assert_allclose(measured, q)
        assert delta[4*ny+4] == 0
        assert abs(delta[0]) < 0.2
    np.testing.assert_allclose(vortex_field(nx, ny, winding=-1),
                               vortex_field(nx, ny, winding=1).conj())


def test_impurity_zero_limit_and_ldos_sum_rule():
    n = 16
    h = thesis_normal(4, 4)
    delta = vortex_field(4, 4)
    grid = np.linspace(-1, 1, 201)
    for magnetic in [False, True]:
        np.testing.assert_array_equal(add_impurity(h, 5, 0, magnetic), h)
        changed = add_impurity(h, 5, 0.43, magnetic)
        assert np.isclose(changed[5, 5]-h[5, 5], 0.43)
        assert np.isclose(changed[n+5, n+5]-h[n+5, n+5], -0.43 if magnetic else 0.43)
        matrix = canonical_bdg(changed, delta)
        assert max(canonical_checks(matrix)) < 1e-14
        energies, vectors = eigh(matrix)
        for site in [0, 5, 15]:
            rho, weight = electron_ldos(energies, vectors, site, grid)
            assert np.all(rho >= 0)
            np.testing.assert_allclose(weight, 2, atol=1e-13)
