"""Structural tests for the operators returned by philsol.core.eigen_build.

These pin down the grid ordering index = i + j*nx, the forward-difference
stencils and the periodic wrap diagonals without needing an eigen-solve,
plus one exact anisotropic mode check.
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
