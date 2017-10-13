"""
Created on Wed Oct 11 09:36:36 2017

@author: pbm24

fd mode solver based on: 
Full-vectorial finite-difference anaylysis of microstructured oprical fibres 
Zhaoming Zhu and Thomous G Brown

"""



import numpy as np
#import numpy.linalg as la
import scipy.constants as cst
import scipy.sparse.linalg as crunch
import scipy.sparse as sps
import time as tempus
        

def eigen_build(k0, n, dx, dy):    
    # lets find out size of grid and construc some finite difference operators    
    nx, ny, dummy = np.shape(n)
    print('Assembling matrix for {} grid points...'.format(nx*ny))
    
    Ux = (-sps.eye(nx*ny, k= 0) + np.eye(nx*ny, k=1)) / dx  
    Uy = (-sps.eye(nx*ny, k= 0) + np.eye(nx*ny, k=nx)) / dy 
    Vx = - Ux.transpose()
    Vy = - Uy.transpose()
    I =  sps.eye(nx*ny)    

    
    # We can build relitice permitivity tensors each point is averaged with its 
    # neighbours.
    epsx = np.zeros((nx*ny, nx*ny))
    epsy = np.zeros((nx*ny, nx*ny))    
    epszi = np.zeros((nx*ny, nx*ny)) 
    count = 0 
    for j in range(0,ny):
        for i in range(0, nx):            
            epsx[count, count] = n[i, j, 0]**2   
            epsy[count, count] = n[i, j, 1]**2
            epszi[count, count] = 1. / n[i, j, 2]**2
            count = count + 1 
            
            
    epsx =  sps.coo_matrix(epsx)
    epsy =  sps.coo_matrix(epsy)   
    epszi = sps.coo_matrix(epszi)
    

    # Now we need to construct the full operator matrices
    t = tempus.time()     
    Pxx = (- Ux * epszi * Vy * Vx * Uy / k0**2 
             + (k0**2 * I + Ux * epszi * Vx) * (epsx + Vy * Uy / k0**2)    )                   
    
    Pyy = (- Uy * epszi * Vx * Vy * Ux / k0**2 
             + (k0**2 * I + Uy * epszi * Vy) * (epsy + Vx * Ux / k0**2)    )
             
    Pxy = (  Ux * epszi * Vy * (epsy + Vx * Ux / k0**2)
              -  (k0**2 * I + Ux * epszi * Vx) * Vy * Ux / k0**2   )
          
    Pyx = (  Uy * epszi * Vx * (epsy + Vy * Uy / k0**2 ) 
             - (k0**2 * I + Uy * epszi * Vy) * Vx * Uy / k0**2 )  
    
    print('and we are done (after {} secs).'.format(tempus.time() - t))     
          
    #print(sps.dia_matrix(Pxx).get_shape())         
    #print(sps.dia_matrix(Pxy).get_shape())
    
    P = sps.vstack( 
           [ sps.hstack([sps.coo_matrix(Pxx), sps.coo_matrix(Pxy)]) , 
                     sps.hstack([sps.coo_matrix(Pyx), sps.coo_matrix(Pyy)])]
                  ) 
   # Ok we should be able to do the final assembly now !!!
    #P = 
    #sps.vstack( [ sps.hstack([Pxx, Pxy], format = "dia"), sps.hstack([Pxx, Pxy], format = "dia") ], format = "dia" )    
    
    return P 

def SolveE(P, beta_trial):
    # Solves eigenproblem
    
    beta, E = crunch.eigsh(P, 1, sigma = beta_trial**2)
    
    Ex, Ey = np.split(E, 2) 
    
    return beta**0.5, Ex, Ey
    
    
      
