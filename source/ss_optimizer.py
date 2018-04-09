'''
Created on Apr 9, 2018

@author: wfg
'''

import scipy.integrate
import scipy.optimize
import numpy as np
from numpy import Inf
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def fmin(ss, ibg, isample):
    return pow(ibg-ss*isample,2)

def ferror(ss, ibg, isample):
    return ibg-ss*isample
    

class ss_optimizer(object):

    def __init__(self, background_file, sample_file, initial_x, final_x, 
                 method="default"):
        
        self.background_file = background_file
        self.sample_file = sample_file
        self.x0 = initial_x
        self.x1 = final_x
        self.method = method
        
        background_data = np.loadtxt(self.background_file)
        sample_data = np.loadtxt(self.sample_file)
        
        self.background_x = background_data[:,0]
        self.background_y = background_data[:,1]
        
        self.sample_x = sample_data[:,0]
        self.sample_y = sample_data[:,1]
        
        self.ibg = self.integrate(self.background_y, self.background_x, self.x0, self.x1)
        self.isample = self.integrate(self.sample_y, self.sample_x, self.x0, self.x1)
    
        self.nfigures = 0
        self.fig = plt.figure()
    
    # public
    def brent(self):
        ss = scipy.optimize.brentq(ferror, 0.5, 1, args=(self.ibg,self.isample))
        return ss
    
    def bisect(self):
        ss = scipy.optimize.bisect(ferror, 0.5, 1, args=(self.ibg,self.isample))
        return ss
    
    def minimize(self):
        ss = scipy.optimize.minimize(fmin, 0.5, args=(self.ibg,self.isample))
        return ss.x
    
    def integrate(self, y, x, initial_x=0, final_x=Inf, plot=False, rule="trapz", ss=1):
          
        if(final_x == Inf ):
            final_x = x[len(x)-1]
            
        li = self.__index_lower_bound(x, initial_x)
        ui = self.__index_upper_bound(x, final_x)
        
        integral = 0   
        
        if(rule == "trapz"):
            integral = scipy.integrate.trapz(y[li:ui], x[li:ui], 1e-10)
        elif(rule == "simps"):
            integral = scipy.integrate.simps(y[li:ui], x[li:ui], 1e-10)
        
        if(plot):
            plt.plot(self.background_x, self.background_y, "black")
            plt.plot(self.sample_x, self.sample_y, "blue")
            ax = self.fig.add_subplot(111)
            self.nfigures = self.nfigures + 1
            plt.plot(x, y, "red")
            plt.yscale("log")
            plt.fill_between(x[li:ui], 0.001, y[li:ui], color="red")
            plt.text(0.2, 0.85, "integral = " + str(integral), ha='center', va='center', 
                     transform = ax.transAxes )
            plt.text(0.2, 0.75, "ss = " + str(ss[0]), ha='center', va='center', 
                     transform = ax.transAxes )
            
            background_legend = mpatches.Patch(color='black', label='Background')
            sample_legend = mpatches.Patch(color='blue', label='Sample')
            scaled_legend = mpatches.Patch(color='red', label='ss * Sample')
            plt.legend(handles=[background_legend, sample_legend, scaled_legend])
            
            plt.show()
            
        return integral
    
    
    def plot(self, ss):

        ssxsample_y = float(ss) * self.sample_y
        self.integrate(ssxsample_y, self.sample_x, self.x0, self.x1, plot=True, ss=ss)
                       
        
    
    
    #private functions
    def __index_lower_bound(self, array, xref):     
        
        i = 0
        for x in np.nditer(array):
            if(x >= xref):
                if(i == 0):
                    return 0    
                else:
                    return i
            else:
                i = i + 1
                    
        return i
    
    
    def __index_upper_bound(self, array, xref):     
        
        i = 0
        for x in np.nditer(array):
            if(x >= xref):
                return i+1
            else:
                i = i + 1
                    
        return i
       
        
