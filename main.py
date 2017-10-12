"""
Created on Thu Oct 12 12:12:29 2017

@author: pbm24

In this example the Ex, Ey feilds are solved for the fundimental mode of a 
square silica waveguide. 

"""


import philsol


# first we need to build a refractive index profile
points = 60
n = np.ones((points, points, 3))
n[19:39, 19:39, :] = 1.45


# Now we do some scaling so the waveguide guides light at telecom frequecies. 
lam = 1.55E-6
k = 2 * cst.pi / lam 
neff = 1.25
beta_in = 2*cst.pi * neff / lam

# we want structure to be scale of wavelength 
dx = 1.E-6 / 20
dy = 1.E-6 / 20

x = np.array(range(points)) * dx 
y = np.array(range(points)) * dy 


# plot the index profile 
plt.pcolor(x*1.E6,y*1.E6,n[:,:,1])

#%%

# Assemble finite difference matrices 
P = eigen_build(k, n, dx, dy)

# Now we solve
beta, Ex, Ey = SolveE(P, beta_in) 
Ex = np.reshape(Ex, (points, points))    
Ey = np.reshape(Ey, (points, points))

print(beta * lam / (2.*cst.pi) )

plt.figure()
plt.pcolor(x*1.E6,y*1.E6,Ex)