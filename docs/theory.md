# Waveguide equation cheat sheet

Assumptions:

- Time-harmonic fields with $e^{-i\omega t}$.
- Non-magnetic material, $\mu = \mu_0$.
- Relative permittivity tensor $\boldsymbol{\varepsilon}(x,y) = \operatorname{diag}\bigl(\varepsilon_x(x,y),\,\varepsilon_y(x,y),\,\varepsilon_z(x,y)\bigr)$, with $\varepsilon_i = n_i^2$.

## Vector wave equation for an anisotropic dielectric

From Maxwell's equations in a source-free region,

$$
\nabla \times \mathbf{E} = i\omega\mu_0 \mathbf{H},\qquad
\nabla \times \mathbf{H} = -i\omega\varepsilon_0\,\boldsymbol{\varepsilon}\,\mathbf{E},
$$

the electric field satisfies

$$
\nabla\times\nabla\times\mathbf{E} - k_0^2\,\boldsymbol{\varepsilon}\,\mathbf{E} = 0,
\qquad k_0 = \frac{\omega}{c} = \frac{2\pi}{\lambda_0}.
$$

For a $z$-invariant waveguide we use the modal ansatz

$$
\mathbf{E}(x,y,z) = \mathbf{e}(x,y)\,e^{i\beta z},\qquad
\mathbf{H}(x,y,z) = \mathbf{h}(x,y)\,e^{i\beta z}.
$$

Substituting this into the wave equation gives a 2-D eigenvalue problem for the transverse profiles $\mathbf{e}_m$, $\mathbf{h}_m$ and the propagation constants $\beta_m$.

## Mode superposition

Any guided field can be written as a superposition of modes:

$$
\mathbf{E}(x,y,z) = \sum_m c_m\,\mathbf{e}_m(x,y)\,e^{i\beta_m z},
$$

$$
\mathbf{H}(x,y,z) = \sum_m c_m\,\mathbf{h}_m(x,y)\,e^{i\beta_m z}.
$$

The coefficients $c_m$ are fixed by the field at $z=0$.

## Effective index

For each mode,

$$
n_{\mathrm{eff},m} = \frac{\beta_m}{k_0}.
$$

## Modal inner product and projection (Poynting normalisation)

A convenient orthogonality relation for lossless, reciprocal waveguides uses the longitudinal component of the time-averaged Poynting vector. With transverse fields $\mathbf{e}_m = (e_{m,x},e_{m,y},0)$ and $\mathbf{h}_m = (h_{m,x},h_{m,y},0)$,

$$
\langle m \mid n \rangle
= \frac{1}{2}\int_A \bigl[\mathbf{e}_m \times \mathbf{h}_n^*\bigr]\cdot\hat{\mathbf{z}}\,dA
= \frac{1}{2}\iint \bigl(e_{m,x}\,h_{n,y}^* - e_{m,y}\,h_{n,x}^*\bigr)\,dx\,dy.
$$

If the modes are normalised,

$$
\langle m \mid n \rangle = \delta_{mn}.
$$

The expansion coefficients of an arbitrary input field $|\psi\rangle$ are then

$$
c_m = \frac{\langle m \mid \psi\rangle}{\langle m \mid m\rangle}
= \frac{\displaystyle\frac{1}{2}\int_A \bigl[\boldsymbol{\psi}\times\mathbf{h}_m^*\bigr]\cdot\hat{\mathbf{z}}\,dA}
{\displaystyle\frac{1}{2}\int_A \bigl[\mathbf{e}_m\times\mathbf{h}_m^*\bigr]\cdot\hat{\mathbf{z}}\,dA}.
$$

If only the electric field $\mathbf{E}_0$ is known, philsol also constructs the corresponding magnetic mode $\mathbf{h}_m$, so the same formula can be used.

### Self-adjointness and the grid dot product

The finite-difference operator built by philsol is self-adjoint under the Poynting-flux inner product, so its eigenvectors are orthogonal with respect to that inner product. On the discrete grid the Poynting integral reduces to a weighted sum over the sampled eigenvector components:

$$
\langle m \mid \psi \rangle \approx \Delta x\,\Delta y \sum_{ij} \mathbf{E}_{m,ij}^* \cdot \boldsymbol{\psi}_{ij}.
$$

This is exactly what is being done in **Cell 7** of `examples/Basic Waveguide Physics.ipynb` (`np.sum(emoji_flat * Ex[i,:].conjugate()) * Ex[i,:]`) to project the input image onto mode $i$.

## Useful extra relations

- Time-averaged Poynting vector:
  $$
  \langle \mathbf{S} \rangle = \frac{1}{2}\operatorname{Re}\bigl[\mathbf{E}\times\mathbf{H}^*\bigr].
  $$

- Out-of-plane magnetic field from Faraday's law:
  $$
  H_z = \frac{i}{k_0}\left(\frac{\partial E_y}{\partial x} - \frac{\partial E_x}{\partial y}\right).
  $$

- Beat length between two modes:
  $$
  L_{\pi} = \frac{\pi}{|\beta_m - \beta_n|}.
  $$
