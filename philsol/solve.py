#import numpy.linalg as la
import scipy.constants as cst
import scipy.sparse.linalg as crunch
import time as tempus
import numpy as np

from philsol.cusolver_complex import csreigvsi

def solve_Et(P, beta_trial, neigs=1):
    # Solves eigenproblem and returns beta and the transverse E-feilds
    print('Solving eigenmodes on CPU')
    t = tempus.time()

    beta_squared, E = crunch.eigs(P, neigs, sigma=beta_trial ** 2)

    Ex, Ey = np.split(E, 2)
    print('{} secs later we have the final solution.'.format(tempus.time() - t))

    return beta_squared ** 0.5, Ex, Ey
