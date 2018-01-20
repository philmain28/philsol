"""
Created on Wed Oct 11 09:36:36 2017

@author: pbm24

fd mode solver based on:
Full-vectorial finite-difference anaylysis of microstructured oprical fibres
Zhaoming Zhu and Thomous G Brown

Please break open the functions and figure out how they work.
"""

from philsol.core import eigen_build

"""
Contains the main maths, lifted from the Zhu paper to assemble Maxwells equations into a linear eigenvalue problem:
- The synatax is as close to verbatim as we could make it
- Returns scipy.sparse matrix which defines the eigen problem for a given refractive index profile. 
- Other operator can be return to help construct other coupled fields (such as E_z, H_x, H_y, H_z)
"""

from philsol import solve

"""
- Handles the eigen solving bit 
- Nothing too complicated here at the moment it just uses standard library scipy sparse libary.
- We hope to add some more fancy eigen solving tricks at a later date (e.g SLEPc )
"""

from philsol import construct

"""
- This function is empty. 
"""
