"""Regression tests for the y-periodic wrap diagonal in philsol.core.eigen_build.

On a vacuum grid with periodic boundaries every mode is a plane wave, so the
fundamental mode sits at exactly n_eff = 1 (constant field, zero transverse
wavevector) and no mode can exceed n_eff = 1. A wrong y-wrap diagonal breaks
the constant-field eigenvector on non-square grids and drags n_eff below 1.
"""

import numpy as np

from philsol.core import eigen_build
from philsol.solve import solve


def vacuum_n_eff(nx: int, ny: int) -> np.ndarray:
    """Solve a vacuum grid, periodic in x and y, and return the mode n_eff values."""
    k0: float = 2 * np.pi / 1.55e-6
    dx: float = 50e-9
    dy: float = 50e-9
    n: np.ndarray = np.ones((nx, ny, 3))
    P, _ = eigen_build(k0, n, dx, dy, x_boundary="periodic", y_boundary="periodic")
    beta, _, _ = solve(P, beta_trial=k0, neigs=2)
    return beta / k0


def test_vacuum_periodic_nonsquare() -> None:
    """Non-square grid: the wrap bug showed up here as n_eff = 0.83 on 40 x 10."""
    n_eff: np.ndarray = vacuum_n_eff(40, 10)
    assert np.min(np.abs(n_eff - 1.0)) < 1e-5
    # Vacuum with periodic boundaries only supports plane waves, n_eff <= 1.
    assert np.max(np.abs(n_eff)) <= 1.0 + 1e-6


def test_vacuum_periodic_square() -> None:
    """Square grid: correct both before and after the wrap fix, guards regression."""
    n_eff: np.ndarray = vacuum_n_eff(40, 40)
    assert np.min(np.abs(n_eff - 1.0)) < 1e-5
    assert np.max(np.abs(n_eff)) <= 1.0 + 1e-6
