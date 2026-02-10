import numpy as np
import scipy as sp
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
        self.pion_mass_blocks = []
        
        self.cosh_model = lambda nt,A0,m : [A0*np.cosh((self.Nt/2.0-nt[i])*m)/np.cosh((self.Nt/2.0-1)*m) for i in range(len(nt))]
    
    def add_corrs(self, name, corrs):
        self.corrs[name] = corrs
    
    def find_E_from_fit(self, corr_avgs, fit_range=[], E_guess=None):
        if(len(fit_range)==0):
            fit_range=[0,len(corr_avgs)]
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
        print(f"{name} energy is {self.energies[name]}/a +- {self.energy_errors[name]}/a")

    def get_energy_div_mpi(self, name):
        if(not (name in self.corrs)):
            print(f"Error: Correlators for {name} have not been saved.")
            return
        if not self.pion_mass_blocks:
            self.pion_mass_blocks = jackknife.get_jackknife_blocks(self.corrs["pipi"], self.block_size, lambda x: self.find_E_from_fit(x)[0])
        blocks = jackknife.get_jackknife_blocks(self.corrs[name], self.block_size, lambda x: self.find_E_from_fit(x)[0])
        blocks = np.divide(blocks, self.pion_mass_blocks)
        mpi = self.find_E_from_fit(np.mean(self.corrs["pipi"],axis=0))[0]
        energies, energy_errors = jackknife.get_errors_from_blocks(self.find_E_from_fit(np.mean(self.corrs[name],axis=0))[0]/mpi, blocks)
        print(f"{name} energy is {energies}*m_pi +- {energy_errors}*m_pi")
        return energies, energy_errors
    
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
        
    def get_fit_matrix_elems(self, t):
        matrix = np.zeros((len(self.matrix),len(self.matrix[0])))
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                if(self.matrix[i][j]=="0"):
                    matrix[i][j] = 0
                else:
                    self.get_energy(self.matrix[i][j])
                    matrix[i][j] = self.cosh_model([t], self.A[self.matrix[i][j]], self.energies[self.matrix[i][j]])[0]
        return matrix
    
    def get_matrix_elems(self, t, matrix, corrs):
        elems = np.zeros((len(matrix),len(matrix[0])))
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if(matrix[i][j]=="0"):
                    elems[i][j] = 0
                else:
                    elems[i][j] = corrs[matrix[i][j]][t]
        return elems
    
    def get_spectrum_from_corrs(self, t, matrix, corrs):
        T1=self.get_matrix_elems(t,matrix,corrs)
        T2=self.get_matrix_elems(t+1,matrix,corrs)
        w,v = sp.linalg.eigh(T1,T2)
        return np.concatenate([[np.log(w)], v / np.sum(v**2,axis=0)**0.5])
    
    def get_spectrum(self, t, divide_pion_mass=False):
        names = list(self.corrs)
        N = len(self.corrs[names[0]])
        corrs = []
        for i in range(N):
            j=0
            corrs.append([])
            for name in names:
                corrs[i].append(self.corrs[name][i])
                j+=1
        matrix=[]
        for i in range(len(self.matrix)):
            matrix.append([])
            for j in range(len(self.matrix[0])):
                if self.matrix[i][j]=="0":
                    matrix[i].append("0")
                else:
                    matrix[i].append(names.index(self.matrix[i][j]))
        if divide_pion_mass:
            if not self.pion_mass_blocks:
                self.pion_mass_blocks = jackknife.get_jackknife_blocks(self.corrs["pipi"], self.block_size, lambda x: self.find_E_from_fit(x)[0])
            av = self.get_spectrum_from_corrs(t,matrix,np.mean(corrs,axis=0))
            blocks = jackknife.get_jackknife_blocks(corrs, self.block_size, lambda x: self.get_spectrum_from_corrs(t,matrix,x))
            for b in range(len(blocks)):
                blocks[b][0] = np.divide(blocks[b][0], self.pion_mass_blocks[b])
            mpi = self.find_E_from_fit(np.mean(self.corrs["pipi"],axis=0))[0]
            av[0] = np.divide(av[0], mpi)
            return jackknife.get_errors_from_blocks(av, blocks)
        else:
            blocks = jackknife.get_jackknife_blocks(corrs, self.block_size, lambda x: self.get_spectrum_from_corrs(t,matrix,x))
            return jackknife.get_errors_from_blocks(self.get_spectrum_from_corrs(t,matrix,np.mean(corrs,axis=0)), blocks)
