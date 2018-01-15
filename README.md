# philsol
## Modes for the Masses (Massless?)
Fed up with relying on expensive proprietary for your waveguide research?  philsol might just be the package for you. 
This is a fully vectorial finite difference waveguide mode solver. In a world where high performance hardware is cheaper than high end software, philsol throws elegence and sophistication out of the window and replaces it with brute force. 

This is a direct Python implimentaion of the algorithm found in the paper: 
['Full-vectorial finite-difference analysis of microstructured optical fibres', by Zhu and Brown.](https://doi.org/10.1364/OE.10.000853)

Warning: I haven't thoroughly tested so be wary and check the results are sensible...

## Installation
- Clone or download, and run `python setup.py install`

## Examples
- Commented example projects can be found in the *examples* directory.
- To run the examples, first install philsol to your Python environment (see above)

## Features
### Solver
- Solves vector Maxwell equations in 2D for arbitary refractive index profile. 
- philsol can handle anisotropic refractive indices with diagonal tensor.
- Currently hard coded with conductive boundary.

### Geometry building
- The quickest way of importing geometry is with a bitmap image (For a tutorial on geometry building https://www.youtube.com/watch?v=QjRi0Mq3G2g)
- See *examples/example_image.py* for an example in loading .bpm images
- See *examples/example_build.py* for an example in building geometry using PIL/Pillow

## To do 
- Solve for z componant and H fields
- More boundry condition options
- Example of RGB as x,y,z refractive indices in imported images
