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

    beta_squared, E = linalg.eigs(P, neigs, sigma=beta_trial ** 2)

    Ex, Ey = np.split(E, 2)
    
    Ex, Ey = np.transpose(Ex), np.transpose(Ey)
    print('{} secs later we have the final solution.'.format(time.time() - t))

    return beta_squared ** 0.5, Ex, Ey
