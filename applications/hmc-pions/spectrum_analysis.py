import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import jackknife

class Spectrum():
    def __init__(self, Nt, block_size):
        self.block_size = block_size
        self.Nt = Nt
        
        self.corrs = {}
        self.energies = {}
        self.energy_errors = {}
        self.A = {}
        self.matrix = [[]]
        
        self.cosh_model = lambda nt,A0,m : [A0*np.cosh((self.Nt/2.0-nt[i])*m)/np.cosh((self.Nt/2.0-1)*m) for i in range(len(nt))]
    
    def add_corrs(self, name, corrs):
        self.corrs[name] = corrs
    
    def find_E_from_fit(self, corr_avgs, fit_range=[], E_guess=None):
        if(len(fit_range)==0):
            fit_range=[0,int(self.Nt)]
        nt = range(fit_range[0],fit_range[1])
        A0_guess = corr_avgs[1]
        if(E_guess==None):
            E_guess = -np.log(corr_avgs[1]/corr_avgs[0])
        a = corr_avgs[fit_range[0]:fit_range[1]]
        pi_opt, pi_cov = curve_fit(self.cosh_model, nt, a, p0=[A0_guess,E_guess])
        return [np.abs(pi_opt[1]), pi_opt[0]]
    
    def get_energy(self, name):
        if(not (name in self.corrs)):
            print(f"Error: Correlators for {name} have not been saved.")
            return
        if(name in self.energies):
            print(f"Already Calculated. {name} energy is {self.energies[name]}/a +- {self.energy_errors[name]}/a")
            return
        blocks = jackknife.get_jackknife_blocks(self.corrs[name], self.block_size, self.find_E_from_fit)
        [[self.energies[name], self.A[name]],[self.energy_errors[name],A_err]] = jackknife.get_errors_from_blocks(self.find_E_from_fit(np.mean(self.corrs[name],axis=0)), blocks)
        print(f"Pion mass is {self.energies[name]}/a +- {self.energy_errors[name]}/a")
    
    def plot_corrs(self, name):
        blocks = jackknife.get_jackknife_blocks(self.corrs[name], self.block_size)
        corr_avgs, corr_errs = jackknife.get_errors_from_blocks(np.mean(self.corrs[name], axis=0), blocks)
        plt.errorbar(range(len(corr_avgs)), corr_avgs, yerr=corr_errs, label=f"{name} Correlators")
    
    def plot_fit(self, name):
        self.get_energy(name)
        self.plot_corrs(name)
        nt_plt = np.arange(0,len(self.corrs[name][0])-1+0.1,0.1)
        plt.plot(nt_plt, self.cosh_model(nt_plt, self.A[name], self.energies[name]), label=f"{name} Fit")
    
    def use_matrix(self, matrix):
        self.matrix = matrix
    
    def get_matrix_elems(self, t):
        matrix = np.zeros((len(self.matrix),len(self.matrix[0])))
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                matrix[i][j] = np.mean(self.corrs[self.matrix[i][j]], axis=0)[t]
        return matrix
