"""
Created on Thu Oct 12 12:12:29 2017

@author: pbm24

In this example the Ex, Ey feilds are solved for the fundimental mode of a 
square silica waveguide. 

"""

import philsol as ps
import numpy as np
import scipy.constants as cst
import matplotlib.pyplot as plt
from PIL import Image

def bmp_to_array(file):
    img = Image.open(file).convert('L')
    img_array = np.asarray(img.getdata(),dtype=np.float64).reshape((img.size[1],img.size[0]))
    img_array /= 255
    
    return img_array

# Lets generate some scale
points = 200

#lets import an sem image
index_array = bmp_to_array('hollow_core3.bmp')

print(np.max(index_array))
print(np.min(index_array))

# Now we do some scaling so the waveguide guides light at telecom frequecies. 
lam = 1.55E-6
k = 2 * cst.pi / lam 
dx = 20.E-6 / 90
dy = 20.E-6 / 90

 
x = np.array(range(points)) * dx 
y = np.array(range(points)) * dy 


n = np.array([index_array,index_array,index_array]) 
n = np.transpose(n)

n = np.ones_like(n) + 0.45 * (n - np.min(index_array)) / (np.max(index_array) - np.min(index_array)) 

# plot the index profile 
plt.figure()
plt.pcolor(x*1.E6,y*1.E6, n[:,:,1], cmap="gray")
plt.show()

#%%

# Assemble finite difference matrices 
P, _ = ps.eigen_build(k, n, dx, dy, False)

#%% Now we solve
neigs = 5
neff = 0.99806
beta_in = 2 * cst.pi * neff / lam

beta, Ex, Ey = ps.solve.solve_Et(P, beta_in, neigs)
for bout in beta: 
    print(bout * lam / (2. * cst.pi) )

E_plot = np.reshape(Ex, (points, points, neigs))

for i in range(neigs):
   s = str(beta[i]* lam / (2. * cst.pi) ) + '.png'
   plt.pcolor(x*1.E6,y*1.E6, np.real(E_plot[:,:,i]))
   plt.show()
