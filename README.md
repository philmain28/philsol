# philsol
## Modes for the Masses (Massless?)
Fed up with relying on expensive proprietary software for your electromagnetic waveguide research? philsol might just be the package for you. In a world where high performance hardware is cheaper than specialist software, philsol throws elegance and sophistication out of the window and replaces them with brute force.

This is a fully vectorial finite-difference waveguide mode solver and a direct Python implementation of the algorithm found in the paper:
['Full-vectorial finite-difference analysis of microstructured optical fibres', by Zhu and Brown.](https://doi.org/10.1364/OE.10.000853)

Warning: I haven't thoroughly tested so be wary and check the results are sensible...

New Warning: The original paper by Zhu and Brown is in Gaussian rather than SI units.
To correct this, use the conversion table [here](https://en.wikipedia.org/wiki/Gaussian_units).

## Installation
- Install using pip with the command `pip install philsol`
- If you can't be bothered, the important part is the function `eigen_build` in `core.py`.

## Examples
- Commented example projects can be found in the *examples* directory.
- To run the examples, first install philsol in your Python environment (see above).

## Features
### Solver
- Solves vector Maxwell (Helmholtz) equations in 2-D for arbitrary refractive-index profiles.
- Returns the x and y components of the electric field.
- philsol can handle anisotropic refractive indices with diagonal tensor.
- Choice of solving routines: the default scipy.sparse solver or SLEPc (slepc4py and petsc4py). These libraries can be fiddly to set up but are very heavily featured, including some limited GPU support.
- Additional field components Ez, Hx, Hy, Hz can be calculated from the construct module.
- Perfect electric conductor, periodic and absorbing boundary conditions.

### Geometry building
- The quickest way of importing geometry is with a bitmap image.
- See *examples/example_image.py* and *examples/Hollow_Core_Fibre.ipynb* for examples loading .bmp images.
- See *examples/example_build.py* for an example in building geometry using PIL/Pillow.
- See *examples/Boundary Interpolation Example.ipynb* for an attempt to handle curved boundaries with Pillow's anti-aliasing capabilities.

## Final Note
I wrote this code a while ago during my [PhD](https://purehost.bath.ac.uk/ws/portalfiles/portal/200802035/philip_main_thesis.pdf) but sometimes I get nostalgic about photonics, so if you are doing anything cool with philsol I would love to know about it.


