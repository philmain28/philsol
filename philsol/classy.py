
"""
Created on Mon Oct  1 14:36:49 2018

Sweet lets wrap up all our functions into a higher level class to save us the 
effort of actualy explicitly passing everything around

The n_build function will need to be manually changed


@author: pbm24
"""
import philsol as ps
import numpy as np

class phil_class:
    def __init__(self, n, k0, x_max = None, y_max = None, dx = None, dy = None):
        self.k0 = k0
        self.n = n
        self.num_x, self.num_y, _ = np.shape(n)
        self.Eigs = 1
        self.E_trial = None
        self.P = None
        
        if x_max == None and y_max == None: 
            if dx != None and dy != None:
               self.dx = dx
               self.dy = dy
               self.x = np.array(range(self.num_x)) * self.dx
               self.y = np.array(range(self.num_y)) * self.dy
            else: 
               raise Exception('Gonna need some dimensions yo!')
               
        if dx == None and dy == None: 
            if x_max != None and y_max != None:
               self.dx = x_max / float(self.num_x - 1)
               self.dy = y_max / float(self.num_y - 1)
            else: 
               raise Exception('Gonna need some dimensions yo!')


    def build_stuff(self, x_bound = None, y_bound = None, matrices = None):
        
        
        if matrices == None:
            self.P, _ = ps.core.eigen_build(    self.k0, 
                                                self.n, 
                                                self.dx, 
                                                self.dy, 
                                                x_boundary = x_bound, 
                                                y_boundary = y_bound      )
        else:
            self.P, self.mats  =  ps.core.eigen_build(    self.k0, 
                                                          self.n, 
                                                          self.dx, 
                                                          self.dy, 
                                                          x_boundary = x_bound, 
                                                          y_boundary = y_bound)
    
    def solve_stuff(self, neigs, beta_trial, extra_fields = None):
        '''
        Function to pass object to solver with the option of constructing the 
        ez and h feilds
        '''
        
        # some insults if everything isn't set up
        if self.P == None: 
            print("Build your eigenproblem you idiot")
            
        
        self.Eigs = neigs
        self.beta, self.Ex, self.Ey = ps.solve.solve(   self.P, 
                                                        beta_trial, 
                                                        E_trial = self.E_trial, 
                                                        neigs = self.Eigs    )

        if extra_fields == True:
            self.E = np.empty(self.neigs, 2 * self.num_x * self.num_y, 3)
            self.H = np.empty_like(self.E)
            for i in self.neigs:
                self.E[i, :, 0] = self.Ex[i,:]
                self.E[i, :, 1] = self.Ey[i,:]
                ez, hx, hy, hz = ps.construct()
                self.E[i, :, 2] = ez
                self.H[i, :, 0] = hx
                self.H[i, :, 1] = hy        
                self.H[i, :, 2] = hz        

        
        
