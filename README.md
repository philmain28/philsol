# philsol
Fed up with relying on expensive proprietary for your waveguide research?  philsol might just be the software for you. 
This is a fully vectorial finite difference waveguide mode solver. In a world where high spec hardware is cheaper than high spec software philsol throws elegence and sophistication out of the window and replaces it with brute force. 
This is a direct python implimentaion of the algorithm found in the paper: 
'Full-vectorial finite-difference analysis of microstructured optical fibres' Zhu and Brown. 
warning: I haven't thoroughly tested so don't trust the results...

## (lack of) features
- Solves vector Maxwell equations in 2D for arbitary refractive index profile. 
- philsol can handle anisotropic refractive indices with diagonal tensor.
- Currently hard coded with conductive boundary.
- see main.py for minimal example
- Best way of importing geometry is with a bitmap image (For a tutorial on geometry building https://www.youtube.com/watch?v=QjRi0Mq3G2g)
- see test1.py for example using SEM image

## To do 
- Solve for z componant and H fields 
- Add option of pyCuda solvers 
- More boundry condition options
