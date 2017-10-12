"""
Created on Wed Oct 11 09:36:36 2017

@author: pbm24

fd mode solver based on: 
Full-vectorial finite-difference anaylysis of microstructured oprical fibres 
Zhaoming Zhu and Thomous G Brown

"""



import numpy as np
import numpy.linalg as la
import scipy.constants as cst
import scipy.sparse.linalg as crunch
import matplotlib.pyplot as plt

def big_mul(mats):        
    # multiplies a big stack of matrices
    a,b = np.shape(mats[0])
    product = np.eye(a)
    for mat in mats:
        product = np.matmul(product, mat)
    
    return product
        

def eigen_build(k0, n, dx, dy):
    # lets find out size of grid and construc some finite difference operators    
    nx, ny, dummy = np.shape(n)
    Ux = (-np.eye(nx*ny, k= 0) + np.eye(nx*ny, k=1)) / dx #(- np.eye(nx*ny, k=-1) + np.eye(nx*ny, k=1)) / (2*dx) 
    Uy = (-np.eye(nx*ny, k= 0) + np.eye(nx*ny, k=nx)) / dy #(- np.eye(nx*ny,k=-nx) + np.eye(nx*ny, k=nx)) / (2.*dy)
    Vx = - np.transpose(Ux)
    Vy = - np.transpose(Uy)
    I =  np.eye(nx*ny)    

    
    # We can build relitice permitivity tensors each point is averaged with its 
    # neighbours.
    epsx = np.zeros((nx*ny, nx*ny))
    epsy = np.zeros((nx*ny, nx*ny))    
    epszi = np.zeros((nx*ny, nx*ny)) 
    count = 0 
    for j in range(0,ny):
        for i in range(0, nx):
            if j == 0 and i == 0:               
                epsx[count, count] = n[i, j, 0]**2   
                epsy[count, count] = n[i, j, 1]**2
                epszi[count, count] = 1. / n[i, j, 2]**2

            elif j == 0 and i != 0:
                epsx[count, count] = n[i, j, 0]**2    
                epsy[count, count] = (n[i, j, 1]**2 + n[i-1, j, 1]**2) / 2.
                epszi[count, count] = 2. / ( n[i, j, 2]**2 + n[i-1 , j, 2]**2 )
                
            elif j != 0 and i== 0:
                epsx[count, count] = (n[i, j, 0]**2 + n[i, j - 1, 0]**2) / 2.   
                epsy[count, count] = n[i, j, 1]**2
                epszi[count, count] = 2. / ( n[i, j, 2]**2 + n[i, j-1, 2]**2 )
                
            else:
                epsx[count, count] = (n[i, j, 0]**2 + n[i, j - 1, 0]**2) / 2.   
                epsy[count, count] = (n[i, j, 1]**2 + n[i-1, j, 1]**2) / 2.
                epszi[count, count] = 4. / ( n[i, j, 2]**2 + n[i-1 , j-1, 2]**2 
                                        +  n[i, j-1, 2]**2 + n[i-1 , j, 2]**2 )

            count = count + 1 

    # Now we need to construct the full operator matrices
    Pxx = (- big_mul([Ux, epszi, Vy, Vx, Uy]) / k0**2 
           + big_mul( [k0**2 * I + big_mul([Ux, epszi, Vx]),              
                                     epsx + big_mul([Vy, Uy]) / k0**2 ])  )
    #Pxx = k0**2 * epsx + np.matmul(Vy, Vy) 
    #Pxy = - np.matmul(Vy, Vx)
    #Pyy = k0**2 * epsx + np.matmul(Vy, Vy)
    
    
    
    #plt.pcolor(Pxx) 
    #plt.show()
    Pyy = (- big_mul([Uy, epszi , Vx, Vy, Ux]) / k0**2 
           + big_mul( [k0**2 * I + big_mul([Uy, epszi, Vy]),              
                                     epsy + big_mul([Vx, Ux]) / k0**2 ])  )
    #plt.pcolor(Pyy)
    #plt.show()
    Pxy = ( big_mul(
            [ big_mul([Ux, epszi, Vy]), epsy + big_mul([Vx, Ux]) / k0**2])
            - big_mul( 
               [ k0**2 * I + big_mul([Ux, epszi, Vx]), big_mul([Vy, Ux]) ]
                                                                   ) / k0**2 )
    #plt.pcolor(Pxy)
    #plt.show()                                                               
    Pyx = ( big_mul(
            [ big_mul([Uy, epszi, Vx]), epsy + big_mul([Vy, Uy]) / k0**2])
            - big_mul( 
               [ k0**2 * I + big_mul([Uy, epszi, Vy]), big_mul([Vx, Uy]) ]
                                                                   ) / k0**2 )    
    #plt.pcolor(Pyx)
    #plt.show()                 
        
    # Ok we should be able to do the final assembly now !!!
    P = np.concatenate(
           (  np.concatenate((Pxx, Pxy), axis = 1)
              ,  np.concatenate((Pyx, Pyy), axis = 1)) , axis = 0 )
    #plt.figure()
    #plt.pcolor(P)    
    
    return P 

def SolveE(P, beta_trial):
    # Solves eigenproblem
    
    beta, E = crunch.eigs(P, 1, sigma = beta_trial**2)
    
    Ex, Ey = np.split(E, 2) 
    
    return beta**0.5, Ex, Ey
    
    
      