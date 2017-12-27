#import numpy.linalg as la
import scipy.constants as cst
import scipy.sparse.linalg as crunch
import time as tempus
import numpy as np

from philsol.cusolver import csreigvsi

def solve_Et(P, beta_trial, neigs=1):
    # Solves eigenproblem and returns beta and the transverse E-feilds
    print('Solving eigenmodes on CPU')
    t = tempus.time()

    beta_squared, E = crunch.eigs(P, neigs, sigma=beta_trial ** 2)

    Ex, Ey = np.split(E, 2)
    print('{} secs later we have the final solution.'.format(tempus.time() - t))

    return beta_squared ** 0.5, Ex, Ey


def solve_Et_cuda(P, beta_trial, neigs=1):
    # Solves eigenproblem and returns beta and the transverse E-feilds
    print('Solving eigenmodes on GPU')
    print("WARNING! THIS WILL CURRENTLY ONLY SOLVE A SINGLE MODE CLOSEST TO TRIAL.\
          BATCH ANALYSIS WILL HOPEFULLY BE IMPLEMENTED SOON!")

    Pcsr = P.tocsr()
    t = tempus.time()

    beta_squared, E = csreigvsi(Pcsr, beta_trial ** 2, neigs=neigs, tol=1e-6, maxite=100)

    Ex, Ey = np.split(E, 2)
    print('{} secs later we have the final solution.'.format(tempus.time() - t))

    return beta_squared ** 0.5, Ex, Ey
