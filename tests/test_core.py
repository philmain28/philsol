"""Structural tests for the operators returned by philsol.core.eigen_build.

These pin down the grid ordering index = i + j*nx, the forward-difference
stencils and the periodic wrap diagonals without needing an eigen-solve,
plus an exact anisotropic mode check and a convergence test against an
analytical solution.
"""

import numpy as np
import scipy.sparse as sps

from philsol.core import eigen_build
from philsol.solve import solve


def simple_n(nx: int, ny: int) -> np.ndarray:
    """Index tensor with a distinct value at every grid point and component."""
    n: np.ndarray = np.empty((nx, ny, 3))
    for i in range(nx):
        for j in range(ny):
            for c in range(3):
                n[i, j, c] = 1.0 + i + 10 * j + 100 * c
    return n


def test_operator_shapes() -> None:
    """P acts on stacked (Ex, Ey) and every operator on the N grid points."""
    nx, ny = 4, 3
    n: np.ndarray = simple_n(nx, ny)
    P, mats = eigen_build(1.0, n, 0.1, 0.1)
    N = nx * ny
    assert sps.issparse(P)
    assert P.shape == (2 * N, 2 * N)
    for key in ("epsx", "epsy", "epszi", "ux", "uy", "vx", "vy"):
        assert mats[key].shape == (N, N)


def test_eps_diagonals_follow_index_ordering() -> None:
    """The eps diagonals must follow index = i + j*nx ordering; the periodic
    wrap diagonal in Uy is derived from the same convention."""
    nx, ny = 4, 3
    n: np.ndarray = simple_n(nx, ny)
    _, mats = eigen_build(1.0, n, 0.1, 0.1)
    expected_x = np.array([n[i, j, 0] ** 2 for j in range(ny) for i in range(nx)])
    expected_y = np.array([n[i, j, 1] ** 2 for j in range(ny) for i in range(nx)])
    expected_zi = np.array([1.0 / n[i, j, 2] ** 2 for j in range(ny) for i in range(nx)])
    assert np.allclose(mats["epsx"].diagonal(), expected_x)
    assert np.allclose(mats["epsy"].diagonal(), expected_y)
    assert np.allclose(mats["epszi"].diagonal(), expected_zi)


def test_uy_forward_difference_and_wrap() -> None:
    """Uy is the forward difference in y; with a periodic boundary the row
    j = ny-1 wraps onto column j = 0, i.e. diagonal -(ny-1)*nx."""
    nx, ny = 3, 2
    dx, dy = 0.1, 0.2
    n: np.ndarray = simple_n(nx, ny)
    _, mats = eigen_build(1.0, n, dx, dy, x_boundary="periodic", y_boundary="periodic")
    Uy: np.ndarray = mats["uy"].toarray()
    expected = np.zeros((nx * ny, nx * ny))
    for i in range(nx):
        for j in range(ny):
            m = i + j * nx
            expected[m, m] = -1.0 / dy
            if j < ny - 1:
                expected[m, m + nx] = 1.0 / dy
            else:
                expected[m, i] = 1.0 / dy
    assert np.allclose(Uy, expected)


def test_ux_forward_difference_and_wrap() -> None:
    """Ux is block diagonal over the y rows; the periodic wrap sits at
    [nx-1, 0] inside each block."""
    nx, ny = 3, 2
    dx, dy = 0.1, 0.2
    n: np.ndarray = simple_n(nx, ny)
    _, mats = eigen_build(1.0, n, dx, dy, x_boundary="periodic", y_boundary="periodic")
    Ux: np.ndarray = mats["ux"].toarray()
    expected = np.zeros((nx * ny, nx * ny))
    for j in range(ny):
        for i in range(nx):
            m = j * nx + i
            expected[m, m] = -1.0 / dx
            if i < nx - 1:
                expected[m, m + 1] = 1.0 / dx
            else:
                expected[m, j * nx] = 1.0 / dx
    assert np.allclose(Ux, expected)


def test_no_wrap_without_periodic_boundaries() -> None:
    """Without 'periodic' the last row of each direction is one-sided and
    carries no wrap entry."""
    nx, ny = 3, 2
    dx, dy = 0.1, 0.2
    n: np.ndarray = simple_n(nx, ny)
    _, mats = eigen_build(1.0, n, dx, dy)
    Uy: np.ndarray = mats["uy"].toarray()
    expected = np.zeros((nx * ny, nx * ny))
    for i in range(nx):
        for j in range(ny):
            m = i + j * nx
            expected[m, m] = -1.0 / dy
            if j < ny - 1:
                expected[m, m + nx] = 1.0 / dy
    assert np.allclose(Uy, expected)
    Ux: np.ndarray = mats["ux"].toarray()
    for j in range(ny):
        assert Ux[j * nx + nx - 1, j * nx] == 0.0


def test_v_operators_are_negative_transposes() -> None:
    """Vx and Vy are the negated transposes of Ux and Uy."""
    n: np.ndarray = simple_n(4, 3)
    _, mats = eigen_build(1.0, n, 0.1, 0.1, x_boundary="periodic", y_boundary="periodic")
    assert np.allclose(mats["vx"].toarray(), -mats["ux"].toarray().T)
    assert np.allclose(mats["vy"].toarray(), -mats["uy"].toarray().T)


def test_uniform_anisotropic_principal_modes() -> None:
    """In a homogeneous diagonal tensor medium with periodic boundaries the
    constant fields along the principal axes are exact modes with
    n_eff = n_component, which checks the eps tensor wiring end to end."""
    nx, ny = 6, 4
    k0: float = 2 * np.pi / 1.55e-6
    n_x, n_y = 1.5, 2.0
    n: np.ndarray = np.empty((nx, ny, 3))
    n[:, :, 0] = n_x
    n[:, :, 1] = n_y
    n[:, :, 2] = 1.0
    P, _ = eigen_build(k0, n, 50e-9, 50e-9, x_boundary="periodic", y_boundary="periodic")
    beta, _, _ = solve(P, beta_trial=k0 * n_x, neigs=4)
    n_eff: np.ndarray = np.abs(beta / k0)
    assert np.min(np.abs(n_eff - n_x)) < 1e-5
    assert np.min(np.abs(n_eff - n_y)) < 1e-5


def test_parabolic_profile_convergence() -> None:
    """Convergence of beta against an exact analytical solution.

    A parabolic permittivity eps(y) = eps0 - c*y^2, invariant along x, turns
    the Ex-polarized mode equation into the quantum harmonic oscillator
    -Ex'' + (k0**2*c)*y^2*Ex = (k0**2*eps0 - beta**2)*Ex, whose fundamental
    eigenvalue is exactly sqrt(k0**2*c) = k0*sqrt(c). Because the profile is
    x-invariant and periodic in x, the discrete operator reproduces this
    scalar problem exactly, so beta must converge to the analytical value at
    second order in the grid spacing.
    """
    k0: float = 2 * np.pi / 1.55
    eps0: float = 9.0
    alpha: float = 9.0  # oscillator strength (um^-2)
    c: float = alpha / k0**2
    beta_exact: float = np.sqrt(k0**2 * eps0 - np.sqrt(alpha))

    errors: list[float] = []
    for ny in (40, 80, 160):
        dy: float = 8.0 / ny
        y: np.ndarray = -4.0 + dy * np.arange(ny)
        profile: np.ndarray = np.sqrt(eps0 - c * y**2)
        n: np.ndarray = np.empty((4, ny, 3))
        n[:, :, 0] = profile
        n[:, :, 1] = profile
        n[:, :, 2] = profile
        P, _ = eigen_build(k0, n, dy, dy, x_boundary="periodic", y_boundary="periodic")
        beta, _, _ = solve(P, beta_trial=beta_exact, neigs=2)
        # The analytical mode must appear in the spectrum; measure how close
        # the nearest computed mode gets (the Ey-polarized twin satisfies a
        # slightly different discrete equation, so take the minimum).
        err: float = float(np.min(np.abs(np.abs(beta) - beta_exact)) / beta_exact)
        errors.append(err)

    # Second order: refining dy by 4x must shrink the error by ~16x.
    order: float = float(np.log2(errors[0] / errors[2]) / 2.0)
    assert 1.5 < order < 2.5
    assert errors[2] < 1e-4
