import numpy as np
import scipy.sparse as sps
import time

def eigen_build(k0, n, dx, dy, operators=False):
    # lets find out size of grid and construct some finite difference operators    
    nx, ny, dummy = np.shape(n)
    print('Assembling matrix for {} grid points...\n'.format(nx*ny))
    
    Ux = (-sps.eye(nx*ny, k= 0) + sps.eye(nx*ny, k=1)) / dx  
    Uy = (-sps.eye(nx*ny, k= 0) + sps.eye(nx*ny, k=nx)) / dy 
    Vx = - Ux.transpose()
    Vy = - Uy.transpose()
    I =  sps.eye(nx*ny)    

    
    # We then build relative permitivity tensors  
    epsx = np.empty(nx*ny)
    epsy = np.empty(nx*ny)    
    epszi = np.empty(nx*ny) 
    count = 0 
    for j in range(0,ny):
        for i in range(0, nx):            
            epsx[count] = n[i, j, 0]**2   
            epsy[count] = n[i, j, 1]**2
            epszi[count] = 1. / n[i, j, 2]**2
            count = count + 1 
            
            
    epsx =  sps.spdiags(epsx, 0, nx*ny,nx*ny)
    epsy =  sps.spdiags(epsy,0, nx*ny,nx*ny)   
    epszi = sps.spdiags(epszi,0, nx*ny,nx*ny)
    

    # Now we need to construct the full operator matrices
    t = time.time()     
    Pxx = (- Ux * epszi * Vy * Vx * Uy / k0**2 
             + (k0**2 * I + Ux * epszi * Vx) * (epsx + Vy * Uy / k0**2)    )                   
    
    Pyy = (- Uy * epszi * Vx * Vy * Ux / k0**2 
             + (k0**2 * I + Uy * epszi * Vy) * (epsy + Vx * Ux / k0**2)    )
             
    Pxy = (  Ux * epszi * Vy * (epsy + Vx * Ux / k0**2)
              -  (k0**2 * I + Ux * epszi * Vx) * Vy * Ux / k0**2   )
          
    Pyx = (  Uy * epszi * Vx * (epsy + Vy * Uy / k0**2 ) 
             - (k0**2 * I + Uy * epszi * Vy) * Vx * Uy / k0**2 )  
    
    print('and we are done (after {} secs).'.format(tempus.time() - t))     
          
    
    # Ok we should be able to do the final assembly now !!!
    P = sps.vstack( 
           [ sps.hstack([sps.coo_matrix(Pxx), sps.coo_matrix(Pxy)]) , 
                     sps.hstack([sps.coo_matrix(Pyx), sps.coo_matrix(Pyy)])]
                  ) 

    return P, {'epsx': epsx, 'epsy':epsy, 'epszi': epszi, 
                   'ux': Ux, 'uy': Uy, 'vx': Vx, 'vy': Vy }

