# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 11:16:50 2017

@author: joel
"""

 ### Interface cuSOLVER PyCUDA


import pycuda.gpuarray as gpuarray
import pycuda.driver as cuda
import pycuda.autoinit
import numpy as np
import scipy.sparse as sp
import ctypes

"""
EVEN MORE BROKEN
"""

# #### wrap the cuSOLVER cusolverSpDcsreigvsi() using ctypes

# cuSparse
_libcusparse = ctypes.cdll.LoadLibrary('cusparse64_91.dll')

class cusparseMatDescr_t(ctypes.Structure):
    _fields_ = [
        ('MatrixType', ctypes.c_int),
        ('FillMode', ctypes.c_int),
        ('DiagType', ctypes.c_int),
        ('IndexBase', ctypes.c_int)
        ]
_libcusparse.cusparseCreate.restype = int
_libcusparse.cusparseCreate.argtypes = [ctypes.c_void_p]

_libcusparse.cusparseDestroy.restype = int
_libcusparse.cusparseDestroy.argtypes = [ctypes.c_void_p]

_libcusparse.cusparseCreateMatDescr.restype = int
_libcusparse.cusparseCreateMatDescr.argtypes = [ctypes.c_void_p]


# cuSOLVER
_libcusolver = ctypes.cdll.LoadLibrary('cusolver64_91.dll')



_libcusolver.cusolverSpCreate.restype = int
_libcusolver.cusolverSpCreate.argtypes = [ctypes.c_void_p]

_libcusolver.cusolverSpDestroy.restype = int
_libcusolver.cusolverSpDestroy.argtypes = [ctypes.c_void_p]



_libcusolver.cusolverSpZcsreigvsi.restype = int

"""
# Enforce data types
_libcusolver.cusolverSpZcsreigvsi.argtypes= [ctypes.c_void_p, #cuso_handle
                                            ctypes.c_int, #m
                                            ctypes.c_int, #nnz
                                            ctypes.c_void_p, #descrA
                                            ctypes.c_void_p, #dcsrVal
                                            ctypes.c_void_p, #dcsrIndPtr
                                            ctypes.c_void_p, #dcsrColInd
                                            ctypes.c_double, #mu0
                                            ctypes.c_void_p, #x0
                                            ctypes.c_int, #maxite
                                            ctypes.c_double, #tol
                                            ctypes.c_void_p, #mu
                                            ctypes.c_void_p] #x 
"""

def csreigvsi(Acsr, mu0, neigs=1, tol=1e-10, maxite=1000):
    
    mu0 = np.complex64(mu0) # Convert mu0 to c_type
    
    mu = np.complex64(mu0) # Create float64 to store mu, seed with mu0
    dmu = cuda.mem_alloc(mu.nbytes) # Allocate device memory to mu
    
    # Initialised sparse matrix, and store csr properties
    dcsrVal = gpuarray.to_gpu(Acsr.data)
    dcsrColInd = gpuarray.to_gpu(Acsr.indices)
    dcsrIndPtr = gpuarray.to_gpu(Acsr.indptr)
    descrA = ctypes.c_void_p()
    
    m = ctypes.c_int(Acsr.shape[0])  # Get row length of sparse matrix
    nnz = ctypes.c_int(Acsr.nnz) # Get number of non-zero components
    
    x0 = np.ones(Acsr.shape[0]).astype(np.complex64) # Seed x0 with ones
    dx0 = gpuarray.to_gpu(x0) # Copy x0 array to device
    
    x = np.empty_like(x0) # Create empty array to store x
    dx = gpuarray.to_gpu(x) # Copy x array to device

    tol = ctypes.c_double(tol) # Convert convergence tollerance to ctype
    maxite = ctypes.c_int(maxite) # Convert maximum iterations to ctype
    
    # Create cusparse handle
    _cusp_handle = ctypes.c_void_p()
    status = _libcusparse.cusparseCreate(ctypes.byref(_cusp_handle))
    assert(status == 0)
    print("CuSparse handle created successfully.")
    cusp_handle = _cusp_handle.value
    
    # Create MatDescriptor
    status = _libcusparse.cusparseCreateMatDescr(ctypes.byref(descrA))
    assert(status == 0)
    print("CuSparse MatDescriptor created successfully.")
    
    # Create cusolver handle
    _cuso_handle = ctypes.c_void_p()
    status = _libcusolver.cusolverSpCreate(ctypes.byref(_cuso_handle))
    assert(status == 0)
    print("CuSolver handle created successfully.")
    cuso_handle = _cuso_handle.value
    
    print('cusp handle: ' + str(cusp_handle))
    print('cuso handle: ' + str(cuso_handle))
    
    
    ### Call solver
    status = _libcusolver.cusolverSpZcsreigvsi(cuso_handle,
                                     m,
                                     nnz,
                                     descrA,
                                     int(dcsrVal.gpudata),
                                     int(dcsrIndPtr.gpudata),
                                     int(dcsrColInd.gpudata),
                                     mu0, #mu0
                                     int(dx0.gpudata), #x0
                                     maxite, #maxite
                                     tol, #tol
                                     int(dmu), #mu
                                     int(dx.gpudata) ) #x
    print("CuSolver ran with status code {}".format(status))
    if status != 0:
        print("A potential error occurred with the solver. \
              Check before trusting these results.")
    
    x = dx.get() # Get eigenvector result as numpy array
    
    # Get mu from device
    # This is a bit tedious, but it works. Will find a better way soon.
    mu = np.array([0.]) # Create temporary array to copy dmu to
    cuda.memcpy_dtoh(mu, dmu) # Copy device mu to host array
    mu = mu[0] # Extract mu float from array
    
    # Destroy handles
    status = _libcusolver.cusolverSpDestroy(cuso_handle)
    print("CuSolver handle destroyed successfully.")
    status = _libcusparse.cusparseDestroy(cusp_handle)
    print("CuSparse handle destroyed successfully.")
    
    return mu, x


if __name__ == "__main__":
    import time
    from scipy.sparse.linalg import eigs
    print("Run speedtest demo")
    #### Prepare the matrix and parameters, copy to Device via gpuarray
    N = 80000
    
    val = np.linspace(10,10*N,N).astype(float)
    col = np.arange(0,N,dtype=np.int32)
    row = np.arange(0,N,dtype=np.int32)
    
    Acoo = sp.coo_matrix((val,(row,col))) # Create sparse matrix in coo format
    Acsr = Acoo.tocsr() # Convert coo to csr
    
    mu0 = np.complex64(425.1) # Initial eigenvalue guess
    
    t0 = time.time()
    mu_gpu, x_gpu = csreigvsi(Acsr, mu0) # Solve
    dt1 = time.time() - t0
    print("CuSolver time: %s" %dt1)

    # Test on CPU
    t0 = time.time()
    mu_cpu, x_cpu = eigs(Acsr, 1, sigma=mu0)
    dt2 = time.time() - t0
    print("SciPyu eigs time: %s" %dt2)

    ratio = dt2/dt1
    if ratio > 1:
        print("CUDA is %s times faster than CPU." %ratio)
    else:
        print("CUDA is %s times slower than CPU." %(1./ratio))
