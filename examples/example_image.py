"""
Created on Thu Oct 12 12:12:29 2017

@author: pbm24

In this example the Ex, Ey feilds are solved 
for 10 modes of a waveguide. 

"""

import philsol as ps
import numpy as np
import scipy.constants as cst
import matplotlib.pyplot as plt
from PIL import Image

def bmp_to_array(file):
    """
    Convert from a gretscale bpm to a 2D numpy array of normalised floats.
    """
    # TODO: Option to use RGB channels and x, y, z refractive indices
    img = Image.open(file).convert('L')
    img_array = np.asarray(img.getdata(),dtype=np.float64).reshape((img.size[1],img.size[0]))
    img_array /= 255
    
    return img_array

# Lets generate some scale
points = 200

#lets import an sem image
index_array = bmp_to_array('image_example.bmp')

# Check that the array has been correctly scaled
assert np.max(index_array) <= 1.0
assert np.min(index_array) >= 0.0

# Create 3*x*y array of the refractive index profile
# In the format [nx[x,y], ny[x,y], nz[x,y]]
n = np.array([index_array,index_array,index_array]) 

# Reshape to an x*x*3 array.
# Effectively a 2D array, where each element is a 3-long list of [nx, ny, nz]
n = np.transpose(n)

# Scale the index array from 1 to 1+dn, where here dn = 0.45
n = np.ones_like(n) + 0.45 * (n - np.min(index_array)) / (np.max(index_array) - np.min(index_array)) 

# Now we do some scaling so the waveguide guides light at telecom frequecies. 
lam = 1.55E-6 # Optical wavelength in m
k = 2 * cst.pi / lam # Optical wavevector calculated from wavelength

# Set the size of a single pixel/mesh-point
dx = 20.E-6 / 90 # Size of a pixel in x-axis, in m
dy = 20.E-6 / 90 # Size of a pixel in y-axis, in m

# Now plot the index profile 
plt.figure()

x = np.array(range(points)) * dx # Create x-axis scaled by dx
y = np.array(range(points)) * dy # Create y-axis scaled by dy

plt.pcolor(x*1.E6,y*1.E6, n[:,:,1], cmap="gray")
plt.show()

#%%

# Assemble finite difference matrices, and discard the operators to _
P, _ = ps.eigen_build(k, n, dx, dy, operators=False)

#%% Now we solve

# Initial guess of effective waveguide refractive index
neff = 1.45

# Calculate effective beta from effective index guess, and wavelength
beta_in = 2 * cst.pi * neff / lam

# Number of eigenmodes to solve for
neigs = 10

# Call philsol solver, to solve the constructed finite-difference matrices
beta, Ex, Ey = ps.solve.solve(P, beta_in, neigs = neigs)

# Print out the effective index for each calculated mode
for bout in beta: 
    print("Effective index: {}".format(bout * lam / (2. * cst.pi) ))

# Create plot data array for x-polarised light
E_plot = np.reshape(Ex, (points, points, neigs))

# For each eigenmode, create a real field plot for x-polarised light
for i in range(neigs):
   plt.figure()
   s = str(beta[i]* lam / (2. * cst.pi) ) + '.png'
   plt.pcolor(x*1.E6, y*1.E6, np.real(E_plot[:,:,i]))
   plt.show()
