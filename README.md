# philsol
Fed up with relying on expensive proprietary for your waveguide research?  philsol might just be the software for you. 
philsol is a fully vectorial finite difference waveguide mode solver which throws elegence and sophistication out of the window and replaces it with brute force. 
This is a direct implimentaion of the algorithm found in the paper: 
'Full-vectorial finite-difference analysis of microstructured optical fibres' Zhu and Brown. 

## (lack of) features
- Solves vector Maxwell equations in 2D for arbitary refractive index profile. 
- philsol can handle anisotropic refractive indixes with diagonal tensor.
- Currently hard coded with conductive boundary.
- see main.py for minimal example

## To do 
- Solve for other fields 
- Cpu optimise (sparse matrix storage etc)
- Add option of pyCuda solvers 
- More boundry condition options
