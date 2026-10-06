# philsol library reference (core workflow)

This reference covers the main public entry points most users need:
`eigen_build`, `solve.solve`, `construct.extra_feilds`, and the
`phil_class` wrapper. For the full source, see the `philsol/` directory.

## `philsol.core.eigen_build`

```python
def eigen_build(
    k0: float,
    n: np.ndarray,
    dx: float,
    dy: float,
    x_boundary: Optional[Literal['periodic']] = None,
    y_boundary: Optional[Literal['periodic']] = None,
) -> Tuple[sps.csr_matrix, Dict[str, sps.csr_matrix]]
```

Build the sparse operator `P` for the vector wave-equation eigenproblem.

### Parameters
- `k0` — Free-space wavenumber, $k_0 = 2\pi / \lambda_0$.
- `n` — Refractive-index array with shape `(nx, ny, 3)`. The last axis is the diagonal tensor components `[n_x, n_y, n_z]`.
- `dx`, `dy` — Grid spacing in x and y.
- `x_boundary`, `y_boundary` — Pass `'periodic'` to enable periodic boundary conditions in that direction.

### Returns
- `P` — Assembled sparse operator as a `scipy.sparse.csr_matrix`.
- `operators` — Dictionary of component matrices/tensors:
  `epsx`, `epsy`, `epszi`, `ux`, `uy`, `vx`, `vy`.

## `philsol.solve.solve`

```python
def solve(
    P: sparse.spmatrix,
    beta_trial: float,
    E_trial: Optional[np.ndarray] = None,
    neigs: int = 1,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]
```

Solve the sparse eigenproblem near `beta_trial`.

### Parameters
- `P` — Sparse operator from `eigen_build`.
- `beta_trial` — Initial guess for the propagation constant.
- `E_trial` — Optional initial guess for the eigenvector (passed to the solver as `v0`).
- `neigs` — Number of eigenvalues/vectors to compute.

### Returns
- `beta` — Propagation constants, shape `(neigs,)`. Effective index is `neff = beta / k0`.
- `Ex` — x-component electric fields, shape `(neigs, nx*ny)`.
- `Ey` — y-component electric fields, shape `(neigs, nx*ny)`.

Field arrays are flattened in row-major (C) order; the grid order used internally is `j*nx + i` (y outer, x inner).

For the optional SLEPc/PETSc solver, see `philsol.solve.solve_fancy`.

## `philsol.construct.extra_feilds`

```python
def extra_feilds(
    k0,
    beta,
    Ex,
    Ey,
    matrices,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
```

> **Note:** the function name is intentionally left as the historic misspelling `extra_feilds`.

Reconstruct the remaining vector field components from the transverse E-fields.

### Parameters
- `k0` — Free-space wavenumber.
- `beta` — Propagation constant for the mode being reconstructed.
- `Ex`, `Ey` — Flattened transverse E-field components.
- `matrices` — Operator dictionary returned by `eigen_build`.

### Returns
- `Ez`, `Hx`, `Hy`, `Hz` — Flattened z-component E-field and full H-field components.

## `philsol.classy.phil_class`

A high-level wrapper that bundles grid definition, matrix assembly, solving,
and field reconstruction.

```python
class phil_class:
    def __init__(
        self,
        n: np.ndarray,
        k0: float,
        x_max: Optional[float] = None,
        y_max: Optional[float] = None,
        dx: Optional[float] = None,
        dy: Optional[float] = None,
    ) -> None
```

Either `(x_max, y_max)` or `(dx, dy)` must be supplied to define the grid.

### Methods

- `build_stuff(x_bound=None, y_bound=None, kx_bloch=0, ky_bloch=0, matrices=False)`  
  Assemble the eigenproblem. Set `matrices=True` if you later need full vector fields.

- `solve_stuff(neigs, beta_trial, extra_fields=False, poynting_vector=False)`  
  Solve for `neigs` modes near `beta_trial`. Set `extra_fields=True` to compute `Ez, Hx, Hy, Hz`; set `poynting_vector=True` (which also requires `extra_fields=True`) to compute the Poynting vector.

- `destroy_crap(fields=False)`  
  Free memory by clearing large internal arrays.

## Minimal example

```python
import numpy as np
import philsol as ps

nx, ny = 120, 120
dx = dy = 0.05
lam = 1.55
k0 = 2 * np.pi / lam

n = np.ones((nx, ny, 3)) * 1.0
n[40:80, 40:80, :] = 1.45

P, _ = ps.eigen_build(k0, n, dx, dy)
beta, Ex, Ey = ps.solve.solve(P, k0 * 1.45, neigs=5)
```

See `docs/theory.md` for the underlying equations and conventions.
