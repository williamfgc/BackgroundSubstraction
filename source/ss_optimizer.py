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
from numpy.f2py.auxfuncs import throw_error

def fmin(ss, fbackground, fsample):
    return pow(ss*fbackground-fsample,2)
    

class ss_optimizer(object):

    def __init__(self, background_file, sample_file, xmin, xmax, xmin_peak, xmax_peak, 
                 method="default"):

        if(xmin_peak > xmax or xmin_peak < xmin):
            raise ValueError('xmin_peak is out of bounds [xmin, xmax]')
        if(xmax_peak > xmax or xmax_peak < xmin):
            raise ValueError('xmax_peak is out of bounds [xmin, xmax]')
        
        self.background_file = background_file
        self.sample_file = sample_file
        self.xmin = xmin
        self.xmax = xmax
        self.xmin_peak = xmin_peak
        self.xmax_peak = xmax_peak
        self.method = method
        
        background_data = np.loadtxt(self.background_file)
        sample_data = np.loadtxt(self.sample_file)
        
        self.background_x = background_data[:,0]
        self.background_y = background_data[:,1]
        
        self.sample_x = sample_data[:,0]
        self.sample_y = sample_data[:,1]
        
        self.fig = plt.figure()
    
    # public
    def run(self):
        
        fIbg = self.fI(self.background_x,self.background_y,
                       self.xmin, self.xmax, 
                       self.xmin_peak, self.xmax_peak)
        
        fIs = self.fI(self.sample_x,self.sample_y,
                       self.xmin, self.xmax, 
                       self.xmin_peak, self.xmax_peak)
        
        ss = scipy.optimize.minimize(fmin, 0.5, args=(fIbg, fIs))
        return ss
    
    # private
    def fI(self, Ix, Iy, xmin, xmax, xmin_peak, xmax_peak):
        
        dx_peak  = xmax_peak - xmin_peak
        dx_left  = xmin_peak - xmin
        dx_right = xmax - xmax_peak
        
        AI_peak  = self.integrate(Iy, Ix, xmin_peak, xmax_peak)
        AI_left  = self.integrate(Iy, Ix, xmin, xmin_peak)
        AI_right = self.integrate(Iy, Ix, xmax_peak, xmax)
        
        fI = AI_peak - 0.5 * dx_peak * (AI_left/dx_left + AI_right/dx_right)
        return fI
    
    
    def integrate(self, y, x, initial_x=0, final_x=Inf, plot=False, rule="trapz", ss=1,
                  log_scale=False):
          
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
            plt.axvline(x=self.xmin_peak, linestyle = 'dashed', color='purple')
            plt.axvline(x=self.xmax_peak, linestyle = 'dashed', color='purple')
            
            if(log_scale):
                plt.yscale("log")
            
            plt.plot(self.background_x[li:ui], self.background_y[li:ui], "black", 
                    linestyle = 'dashed')
            plt.plot(self.sample_x[li:ui], self.sample_y[li:ui], "red")
            ax = self.fig.add_subplot(111)
            
            plt.plot(x[li:ui], y[li:ui], "green")
            plt.fill_between(x[li:ui], 0.001, y[li:ui], color="green")
            plt.text(0.2, 0.85, "integral = " + str(integral), ha='center', va='center', 
                     transform = ax.transAxes )
            plt.text(0.2, 0.75, "ss = " + str(ss[0]), ha='center', va='center', 
                     transform = ax.transAxes )
            
            background_legend = mpatches.Patch(color='black', label='Background')
            sample_legend = mpatches.Patch(color='red', label='Sample')
            scaled_legend = mpatches.Patch(color='green', label='ss * Background')
            plt.legend(handles=[background_legend, sample_legend, scaled_legend])
            
            plt.show()
            
        return integral
    
    
    def plot(self, ss):

        ssxbackground_y = float(ss) * self.background_y
        self.integrate(ssxbackground_y, self.background_x, self.xmin, self.xmax, plot=True, ss=ss)
    
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
       
        
