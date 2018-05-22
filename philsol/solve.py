import scipy.constants as cst
import scipy.sparse.linalg as linalg
import time
import numpy as np

def solve(P, beta_trial, E_trial=None, neigs=1):
    """
	Solves eigenproblem and returns beta and the transverse E-feilds
	"""
    print('Solving eigenmodes on CPU')
    t = time.time()

    beta_squared, E = linalg.eigs(P, neigs, sigma=beta_trial ** 2, v0 = E_trial)

    Ex, Ey = np.split(E, 2)
    
    Ex, Ey = np.transpose(Ex), np.transpose(Ey)
    print('{} secs later we have the final solution.'.format(time.time() - t))

    return beta_squared ** 0.5, Ex, Ey

def solve_fancy(P, beta_trial, E_trial=None, neigs=1):
    """
	Solves eigenproblem with fancier tools
	"""
    from petsc4py import PETSc
    from slepc4py import SLEPc
    

    t = time.time()
    
    # convert eigenproblem into a form that petsc can understand
    fancy_P = PETSc.Mat().createAIJ(size=P.shape,
                                      csr=(P.indptr, P.indices,
                                           P.data))
    
    #define a spectral transform
    ''''
    trans = SLEPc.ST(); trans.create(MPI_COMM_WORLD)
    trans.setShift(beta_trial**2)
    '''
    
    
    # initalise solver object
    E = SLEPc.EPS(); E.create()
    
    # lets set the solver options
    #E.setST(trans)
    E.setOperators(fancy_P)
    E.setProblemType(SLEPc.EPS.ProblemType.NHEP)
    E.setDimensions(nev = neigs)
    if E_trial != None:
        E.setInitialSpace(E_trial)

    E.setTarget(beta_trial**2)
    E.setWhichEigenpairs(7)
    #E.setTarget(1.E20)
    #now we set up the spectral transform
    
    
    #beta_trial**2)
    print('Solving eigenmodes using fancy solver')
    E.solve()
    
    nconv = E.getConverged()
    beta = []
    Ex = []
    Ey = []
    vr, wr = fancy_P.getVecs()
    vi, wi = fancy_P.getVecs()
    for i in range(nconv):
        beta.append(E.getEigenpair(i, vr, vi)**0.5)
        #beta.append(E.getEigenvalue(i)**0.5)
        Exr, Eyr = np.split(np.array(vr),2)
        Exi, Eyi = np.split(np.array(vi),2)
        Ex.append(Exr + 1j * Exi) 
        Ey.append(Eyr + 1j * Eyi)
        
        
     
    E.destroy()
    fancy_P.destroy()
    #now we can start unpacking
    return beta, Ex, Ey #beta_squared ** 0.5, Ex, Ey
