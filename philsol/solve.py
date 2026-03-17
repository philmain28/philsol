import scipy.sparse.linalg as linalg
import time
import numpy as np
from typing import Tuple, List, Optional, Union
import scipy.sparse as sparse

def solve(P: sparse.spmatrix, beta_trial: float, E_trial: Optional[np.ndarray] = None, neigs: int = 1) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solves eigenproblem and returns beta and the transverse E-fields
    
    Args:
        P: Sparse matrix for the eigenvalue problem
        beta_trial: Initial guess for eigenvalue
        E_trial: Initial guess for eigenvector (optional)
        neigs: Number of eigenvalues/vectors to compute
        
    Returns:
        Tuple containing:
            - Propagation constants (square root of eigenvalues)
            - Ex field components
            - Ey field components
    """
    print('Solving eigenmodes on CPU')
    t = time.time()

    beta_squared, E = linalg.eigs(P, neigs, sigma=beta_trial ** 2, v0 = E_trial)

    Ex, Ey = np.split(E, 2)
    
    Ex, Ey = np.transpose(Ex), np.transpose(Ey)
    print('{} secs later we have the final solution.'.format(time.time() - t))

    return beta_squared ** 0.5, Ex, Ey

def solve_fancy(P: sparse.spmatrix, beta_trial: float, E_trial: Optional[np.ndarray] = None, neigs: int = 1) -> Tuple[List[complex], List[np.ndarray], List[np.ndarray]]:
    """
    Solves eigenproblem with PETSc and SLEPc solvers, in theory these should
    give better performance but require a bit more set up. 
    
    Args:
        P: Sparse matrix for the eigenvalue problem
        beta_trial: Initial guess for eigenvalue
        E_trial: Initial guess for eigenvector (optional)
        neigs: Number of eigenvalues/vectors to compute
        
    Returns:
        Tuple containing:
            - List of propagation constants (square root of eigenvalues)
            - List of Ex field components
            - List of Ey field components
    """
    from petsc4py import PETSc
    from slepc4py import SLEPc
    

    t = time.time()
    
    # convert eigenproblem into a form that petsc can understand
    fancy_P = PETSc.Mat().createAIJ(size=P.shape,
                                      csr=(P.indptr, P.indices,
                                           P.data))
    
    
    # initalise solver object
    E = SLEPc.EPS(); E.create()
    
    # lets set the solver options
    E.setOperators(fancy_P)
    E.setProblemType(SLEPc.EPS.ProblemType.NHEP)
    E.setDimensions(nev = neigs)
    if E_trial != None:
        E.setInitialSpace(E_trial)

    # now we set up the spectral region to look in 	
    E.setTarget(beta_trial**2) 
    E.setWhichEigenpairs(7) #look for closest in absolute value
    print('Solving eigenmodes using fancy solver')
    E.solve()
    
    #now we can start unpacking	
    nconv = E.getConverged()
    beta = []
    Ex = []
    Ey = []
    vr, _ = fancy_P.getVecs()
    vi, _ = fancy_P.getVecs()
    for i in range(nconv):
        beta.append(E.getEigenpair(i, vr, vi)**0.5)
        #beta.append(E.getEigenvalue(i)**0.5)
        Exr, Eyr = np.split(np.array(vr),2)
        Exi, Eyi = np.split(np.array(vi),2)
        Ex.append(Exr + 1j * Exi) 
        Ey.append(Eyr + 1j * Eyi)
        
        
     
    E.destroy()
    fancy_P.destroy()
    
    return beta, Ex, Ey 
