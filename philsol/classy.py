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
        self.mats = None
        
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
               '''
               self.dx = dx
               self.dy = dy
               self.x = np.array(range(self.num_x)) * self.dx
               self.y = np.array(range(self.num_y)) * self.dy
               '''
               self.dx = x_max / float(self.num_x - 1)
               self.dy = y_max / float(self.num_x - 1)
            else: 
               raise Exception('Gonna need some dimensions yo!')


    def build_stuff(self, x_bound = None, y_bound = None, kx_bloch = 0, ky_bloch = 0, matrices = None):
        
        
        if matrices == None:
            self.P, _ = ps.core.eigen_build(    self.k0, 
                                                self.n, 
                                                self.dx, 
                                                self.dy, 
                                                x_boundary = x_bound, 
                                                y_boundary = y_bound     
                                                )

        else:
            self.P, self.mats  =  ps.core.eigen_build(    self.k0, 
                                                          self.n, 
                                                          self.dx, 
                                                          self.dy, 
                                                          x_boundary = x_bound, 
                                                          y_boundary = y_bound
                                                          )
    
    def solve_stuff(self, neigs, beta_trial, extra_fields = False, poynting_vector = False):
        '''
        Function to pass object to solver with the option of constructing the 
        ez and h feilds
        '''
        
        # some insults if everything isn't set up
        if self.P == None: 
            raise Exception("Build your eigenproblem you idiot")
            
        if self.mats == None and extra_fields == True:
            raise Exception("If full vector fields are required the build_stuff step needs matrices = True")
        

        self.Eigs = neigs
        self.beta, self.Ex, self.Ey = ps.solve.solve(   self.P, 
                                                        beta_trial, 
                                                        E_trial = self.E_trial, 
                                                        neigs = self.Eigs    
                                                        )

        if extra_fields == True:

            self.E = np.empty((self.Eigs, self.num_x * self.num_y, 3), dtype = complex)
            self.H = np.empty_like(self.E , dtype = complex )

            for i in range(self.Eigs):
                self.E[i, :, 0] = self.Ex[i,:]
                self.E[i, :, 1] = self.Ey[i,:]
                ez, hx, hy, hz = ps.construct.extra_feilds(self.k0, self.beta[i], self.Ex[i,:], self.Ey[i,:], self.mats)
                self.E[i, :, 2] = ez
                self.H[i, :, 0] = hx
                self.H[i, :, 1] = hy        
                self.H[i, :, 2] = hz  
        
        if poynting_vector == True:
            if not extra_fields:
                raise Exception("Poynting vector requires full E and H fields")
        
            self.S = np.empty((self.Eigs, self.num_x * self.num_y, 3), dtype = complex)
            
            # Calculate Poynting vector components: S = E × H*
            self.S[:, :, 0] = self.E[:, :, 1] * np.conj(self.H[:, :, 2]) - self.E[:, :, 2] * np.conj(self.H[:, :, 1])
            self.S[:, :, 1] = self.E[:, :, 2] * np.conj(self.H[:, :, 0]) - self.E[:, :, 0] * np.conj(self.H[:, :, 2])
            self.S[:, :, 2] = self.E[:, :, 0] * np.conj(self.H[:, :, 1]) - self.E[:, :, 1] * np.conj(self.H[:, :, 0])
                

    
    def destroy_crap(self, fields = False):
        
        self.n = None
        self.P = None
        self.mats = None
        self.Ex = None 
        self.Ey = None
        
        if fields == True: 
            self.E = None
            self.H = None
        