from typing import Optional, Tuple, Union, List
import philsol as ps
import numpy as np


class phil_class:
    """Helper class to wrap up some of the functionality of the philsol library and avoid lots of boilerplate"""
    def __init__(self, n: np.ndarray, k0: float, x_max: Optional[float] = None, y_max: Optional[float] = None, dx: Optional[float] = None, dy: Optional[float] = None) -> None:
        """Initialize the electromagnetic solver with material and dimensional parameters.
        
        Args:
            n: Refractive index array with shape (num_x, num_y, ...)
            k0: Free-space wavenumber
            x_max: Maximum x dimension (used with y_max to calculate grid spacing)
            y_max: Maximum y dimension (used with x_max to calculate grid spacing)
            dx: Grid spacing in x direction (alternative to x_max)
            dy: Grid spacing in y direction (alternative to y_max)
            
        Raises:
            Exception: If neither (x_max, y_max) nor (dx, dy) are provided
        """
        self.k0 = k0
        self.n = n
        self.num_x, self.num_y, _ = np.shape(n)
        self.Eigs = 1
        self.E_trial = None
        self.P = None
        self.mats = None
        
        if x_max == None and y_max == None: 
            if dx != None and dy != None:
               self.dx = dx
               self.dy = dy
               self.x = np.array(range(self.num_x)) * self.dx
               self.y = np.array(range(self.num_y)) * self.dy
            else: 
               raise Exception('Gonna need some dimensions yo!')
               
        if dx == None and dy == None: 
            if x_max != None and y_max != None:
               '''
               self.dx = dx
               self.dy = dy
               self.x = np.array(range(self.num_x)) * self.dx
               self.y = np.array(range(self.num_y)) * self.dy
               '''
               self.dx = x_max / float(self.num_x - 1)
               self.dy = y_max / float(self.num_x - 1)
            else: 
               raise Exception('Gonna need some dimensions yo!')


    def build_stuff(self, x_bound: Optional[str] = None, y_bound: Optional[str] = None, kx_bloch: float = 0, ky_bloch: float = 0, matrices: Optional[bool] = None) -> None:
        """Build the eigenvalue problem matrices.
        
        Constructs the matrices needed to solve the eigenvalue problem using
        the core.eigen_build function.
        
        Args:
            x_bound: Boundary condition type in x direction
            y_bound: Boundary condition type in y direction
            kx_bloch: Bloch wave vector in x direction
            ky_bloch: Bloch wave vector in y direction
            matrices: If True, store additional matrices needed for field calculations
        """
        
        
        if matrices == None:
            self.P, _ = ps.core.eigen_build(    self.k0, 
                                                self.n, 
                                                self.dx, 
                                                self.dy, 
                                                x_boundary = x_bound, 
                                                y_boundary = y_bound     
                                                )

        else:
            self.P, self.mats  =  ps.core.eigen_build(    self.k0, 
                                                          self.n, 
                                                          self.dx, 
                                                          self.dy, 
                                                          x_boundary = x_bound, 
                                                          y_boundary = y_bound
                                                          )
    
    def solve_stuff(self, neigs: int, beta_trial: float, extra_fields: bool = False, poynting_vector: bool = False) -> None:
        """Solve the eigenvalue problem and optionally calculate field components.
        
        Computes eigenvalues (propagation constants) and eigenvectors (field distributions)
        for the electromagnetic modes, with options to calculate the full vector fields
        and Poynting vector.
        
        Args:
            neigs: Number of eigenvalues/eigenvectors to compute
            beta_trial: Initial guess for the propagation constant
            extra_fields: If True, calculate the full E and H field components
            poynting_vector: If True, calculate the Poynting vector (requires extra_fields=True)
            
        Raises:
            Exception: If P is not built or if fields are requested without matrices
        """
        
        # some insults if everything isn't set up
        if self.P == None: 
            raise Exception("Build your eigenproblem you idiot")
            
        if self.mats == None and extra_fields == True:
            raise Exception("If full vector fields are required the build_stuff step needs matrices = True")
        

        self.Eigs = neigs
        self.beta, self.Ex, self.Ey = ps.solve.solve(   self.P, 
                                                        beta_trial, 
                                                        E_trial = self.E_trial, 
                                                        neigs = self.Eigs    
                                                        )

        if extra_fields == True:

            self.E = np.empty((self.Eigs, self.num_x * self.num_y, 3), dtype = complex)
            self.H = np.empty_like(self.E , dtype = complex )

            for i in range(self.Eigs):
                self.E[i, :, 0] = self.Ex[i,:]
                self.E[i, :, 1] = self.Ey[i,:]
                ez, hx, hy, hz = ps.construct.extra_feilds(self.k0, self.beta[i], self.Ex[i,:], self.Ey[i,:], self.mats)
                self.E[i, :, 2] = ez
                self.H[i, :, 0] = hx
                self.H[i, :, 1] = hy        
                self.H[i, :, 2] = hz  
        
        if poynting_vector == True:
            if not extra_fields:
                raise Exception("Poynting vector requires full E and H fields")
        
            self.S = np.empty((self.Eigs, self.num_x * self.num_y, 3), dtype = complex)
            
            # Calculate Poynting vector components: S = E × H*
            self.S[:, :, 0] = self.E[:, :, 1] * np.conj(self.H[:, :, 2]) - self.E[:, :, 2] * np.conj(self.H[:, :, 1])
            self.S[:, :, 1] = self.E[:, :, 2] * np.conj(self.H[:, :, 0]) - self.E[:, :, 0] * np.conj(self.H[:, :, 2])
            self.S[:, :, 2] = self.E[:, :, 0] * np.conj(self.H[:, :, 1]) - self.E[:, :, 1] * np.conj(self.H[:, :, 0])
                

    
    def destroy_crap(self, fields: bool = False) -> None:
        """Free memory by setting large data structures to None.
        
        Args:
            fields: If True, also clear field-related data structures
        """
        
        self.n = None
        self.P = None
        self.mats = None
        self.Ex = None 
        self.Ey = None
        
        if fields == True: 
            self.E = None
            self.H = None
        