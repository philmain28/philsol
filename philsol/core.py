import numpy as np
import scipy.sparse as sps
import time
from typing import Dict, Tuple, Optional, Literal

def eigen_build(k0: float, n: np.ndarray, dx: float, dy: float, 
               x_boundary: Optional[Literal['periodic']] = None, 
               y_boundary: Optional[Literal['periodic']] = None) -> Tuple[sps.csr_matrix, Dict[str, sps.csr_matrix]]:
    """Build eigenvalue problem matrices for electromagnetic field simulation.
    
    Constructs the operator matrices for solving electromagnetic wave equations on a 2D grid
    using finite difference frequency domain method. The function handles different boundary conditions
    and builds the necessary differential operators and material property tensors.
    
    Args:
        k0: Wavenumber in free space
        n: 3D array of refractive indices with shape (nx, ny, 3), where the last dimension
           represents the tensor components
        dx: Grid spacing in x-direction
        dy: Grid spacing in y-direction
        x_boundary: Boundary condition for x-direction, 'periodic' for periodic boundaries
        y_boundary: Boundary condition for y-direction, 'periodic' for periodic boundaries
        
    Returns:
        Tuple containing:
            - P: The assembled operator matrix as a sparse CSR matrix
            - operators: Dictionary of component operators and tensors used in the calculation
              (epsx, epsy, epszi, ux, uy, vx, vy)
    """
    
    # detect dimensions of refractive index grid    
    nx, ny, _ = np.shape(n)
    print('Assembling matrix for {} grid points...\n'.format(nx*ny))
    

    
    # construct finite difference operators single row of FD grid 
    Ux_temp = ( - np.eye(nx, k = 0) + np.eye(nx, k = 1) ) / dx
    
    Uy  = ( - sps.eye(nx*ny, k = 0) + sps.eye(nx*ny, k = nx) ) / dy
  
    if x_boundary == 'periodic':
        Ux_temp[nx-1, 0] = 1. / dx
        
    # Wrap around periodic boundary. Something something bloch functions.
    if y_boundary == 'periodic':
        Uy = Uy + sps.eye(nx*ny, k= -(nx-1)*ny) / dy
    
    
    #This statement is kind of confusing but is the equivalent to doing a tensor 
    #contraction. So each row operation is appended to the diagonals of a larger 
    #matrix so we can operate on the whole grid at once.  
    Ux = sps.block_diag([Ux_temp for _ in range(ny)], format = 'csr')
    
    # Construct other differential operators 
    Vx = - Ux.transpose()
    Vy = - Uy.transpose()
    I =  sps.eye(nx*ny)


    # We then build relative permitivity tensors  
    epsx = np.empty(nx*ny, dtype=n.dtype)
    epsy = np.empty(nx*ny, dtype=n.dtype)
    epszi = np.empty(nx*ny, dtype=n.dtype)
    count = 0 
    for j in range(0,ny):
        for i in range(0, nx):            
            epsx[count] = n[i, j, 0]**2   
            epsy[count] = n[i, j, 1]**2
            epszi[count] = 1. / n[i, j, 2]**2
            count = count + 1 
            
            
    epsx =  sps.spdiags(epsx, 0, nx*ny, nx*ny, format = 'csr')
    epsy =  sps.spdiags(epsy,0, nx*ny, nx*ny, format = 'csr')   
    epszi = sps.spdiags(epszi,0, nx*ny, nx*ny, format = 'csr')
    

    # Now we need to construct the full operator matrices
    t = time.time()     
    Pxx = (- Ux * epszi * Vy * Vx * Uy / k0**2 
             + (k0**2 * I + Ux * epszi * Vx) * (epsx + Vy * Uy / k0**2)    )                   
    
    Pyy = (- Uy * epszi * Vx * Vy * Ux / k0**2 
             + (k0**2 * I + Uy * epszi * Vy) * (epsy + Vx * Ux / k0**2)    )
             
    Pxy = (  Ux * epszi * Vy * (epsy + Vx * Ux / k0**2)
              -  (k0**2 * I + Ux * epszi * Vx) * Vy * Ux / k0**2   )
          
    Pyx = (  Uy * epszi * Vx * (epsx + Vy * Uy / k0**2 ) 
             - (k0**2 * I + Uy * epszi * Vy) * Vx * Uy / k0**2 )  
    
    print('and we are done (after {} secs).'.format(time.time() - t))     
          
    
    # Ok we should be able to do the final assembly now !!!
    P = sps.vstack( [ sps.hstack([Pxx, Pxy]) , sps.hstack([Pyx, Pyy])] ) 

    return P, {'epsx': epsx, 'epsy':epsy, 'epszi': epszi, 
                                       'ux': Ux, 'uy': Uy, 'vx': Vx, 'vy': Vy }

