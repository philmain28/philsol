# philsol
## Modes for the Masses (Massless?)
Fed up with relying on expensive proprietary software for your electromagnetic waveguide research?  philsol might just be the package for you. 
This is a fully vectorial finite difference waveguide mode solver. In a world where high performance hardware is cheaper than specialist software (including certain commercial packages which may or may not rhyme with my code), philsol throws elegence and sophistication out of the window and replaces it with brute force. 

This is a direct Python implimentaion of the algorithm found in the paper: 
['Full-vectorial finite-difference analysis of microstructured optical fibres', by Zhu and Brown.](https://doi.org/10.1364/OE.10.000853)

Warning: I haven't thoroughly tested so be wary and check the results are sensible...

New Warning: Original paper by Zhu and Brown is in gaussian not S. I. units. 
This means that the constructed H fields must be corrected by a factor of $\mu_0 c$ for calculations in S.I units. 

## Installation
- Intall using pip with command 'pip install philsol'
- If you can be bothered, the important part is the function eigenbuild in core.py. 

## Examples
- Commented example projects can be found in the *examples* directory.
- To run the examples, first install philsol to your Python environment (see above)

## Features
### Solver
- Solves vector Maxwell(Helmholtz) equations in 2D for arbitary refractive index profile.
- Return x and y componants of electric field.
- philsol can handle anisotropic refractive indices with diagonal tensor.
- Currently hard coded with conductive boundary.
- Now includes choice of eigensolver: the default scipy.sparse solver based on ARPACK and Petsc (but you will need a working install of slepc4py and petsc4py)
- Extra field componants Ez, Hx, Hy, Hz can be calculated from construct module
- Periodic boundary conditions 

### Geometry building
- The quickest way of importing geometry is with a bitmap image 
- See *examples/example_image.py* for an example in loading .bpm images
- See *examples/example_build.py* for an example in building geometry using PIL/Pillow

## To do 
- More intelligent geometry aproximation (e.g pixel interpolation on curved boundaries)
- More boundry condition options Bloch, PML...
- GPU eigensolving


