import glob
import pickle
import matplotlib.pyplot as plt
import numpy as np
import jackknife as jk
from scipy.integrate import quad_vec
from scipy.optimize import curve_fit

import ratios_fit

class Analysis:
    def __init__(self, data):
        self.data = data
    
    # Functions for calculating decay rate =======================
    
    # Calculating M-L ratios ---------------
    
    def get_M_L_blocks(self, Ms, Ls, params, der=False):
        sf_ML = self.data.get_indices(params)[0]
        sfs_M = self.data.replace_params(sf_ML, ["M", "L"], [[M, 1.0] for M in Ms])
        sfs_L = self.data.replace_params(sf_ML, ["M", "L"], [[1.0, L] for L in Ls])
        delta_actions_M = [jk.get_jackknife_blocks(np.exp(self.data.delta_actions_M[sfs_M[i]][str(Ms[i+1])][self.data.cutoff:]), self.data.block_size)
                           for i in range(len(Ms)-1)]
        delta_actions_L = [jk.get_jackknife_blocks(np.exp(self.data.delta_actions_L[sfs_L[i]][str(Ls[i+1])][self.data.cutoff:]), self.data.block_size)
                           for i in range(len(Ls)-1)]
        if(der):
            assert(Ls[0] == 0.0)
            delta_actions_L[0] = jk.get_jackknife_blocks(np.multiply(self.data.delta_actions_t_FV[sfs_L[0]][f"{int(self.data.params[sfs_L[0]]["tFV"])+1}"][self.data.cutoff:],np.exp(self.data.delta_actions_L[sfs_L[0]][str(Ls[1])][self.data.cutoff:])), self.data.block_size)
        return delta_actions_M + delta_actions_L
    
    def calc_ratio(self, delta_actions, N_Ms):
        ratio = 1.0
        for i in range(len(delta_actions)):
            if(i<N_Ms):
                ratio *= np.mean(delta_actions[i])
            else:
                ratio /= np.mean(delta_actions[i])
        return ratio

    def get_R_blocks(self, Ms, Ls, params, der = False):
        ml_blocks = self.get_M_L_blocks(Ms, Ls, params, der)
        return jk.super_jackknife_combine_blocks(ml_blocks, lambda x: self.calc_ratio(x, len(Ms)-1))
    
    # Calculating fit ----------------------
    
    def fit_ratios(self, ratios, ratio_errs, t_TVs, fit_time, dt, start=0, stop=100, fitobject=ratios_fit.GaussianFit):
        fit = fitobject(dt)
        return fit.fit_correction(fit_time, t_TVs, ratios, ratio_errs, filter_x=lambda t: t<start or t>stop)
    
    def get_fit_ratios_blocks(self, params, start=0, stop=100, fitobject=ratios_fit.GaussianFit):
        t_TV = int(params["tTV"])*float(params["dt"])
        params_tTV = params.copy()
        params_tTV["M"] = 1.0
        params_tTV["L"] = 0.0
        del params_tTV["tTV"]
        #
        sfs = self.data.get_indices(params_tTV)
        sfs.sort(key=lambda x: float(self.data.params[x]["tFV"]))
        dt = float(self.data.params[sfs[0]]["dt"])
        #
        dS_blocks = []
        errs = []
        t_TVs = []
        for sf in sfs:
            t = int(self.data.params[sf]["tTV"])*dt
            if(t>start):
                dS_blocks.append(self.get_exp_Ebar_blocks(sf))
                errs.append(jk.get_errors_from_blocks(np.mean(dS_blocks[-1]),dS_blocks[-1])[1])
                t_TVs.append(t)
        
        mean_fit = self.fit_ratios(np.mean(dS_blocks,axis=1),errs,t_TVs,t_TV,dt,start,stop,fitobject)
        print(f"Based on mean, ratios fit is: {mean_fit}")
        
        return mean_fit, jk.super_jackknife_combine_blocks(dS_blocks, lambda x: self.fit_ratios(x,errs,t_TVs,t_TV,dt,start,stop,fitobject))
    
    # Calculating decay rate using fit------
    
    def calc_gamma(self, R, correction_factor, t_full, dt):
        return 2*np.pi*R*correction_factor/(t_full*dt)**2
    
    def calc_gamma_blocks(self, Ms, Ls, fit_start, fit_stop, params, der=False, fitobject=ratios_fit.GaussianFit):
        R_blocks = self.get_R_blocks(Ms, Ls, params, der)
        fit_from_mean, fit_blocks = self.get_fit_ratios_blocks(params, start=fit_start, stop=fit_stop, fitobject=fitobject)
        #
        print(f"Correction factor estimated: {fit_from_mean}")
        print(f"R ratio estimated: {np.mean(R_blocks)}")
        #
        t_full = int(params["tfull"])
        dt = float(params["dt"])
        gamma_blocks = jk.super_jackknife_combine_blocks([R_blocks, fit_blocks], lambda x: self.calc_gamma(x[0], x[1], t_full, dt))
        gamma_mean = self.calc_gamma(np.mean(R_blocks), fit_from_mean, t_full, dt)
        return gamma_mean, gamma_blocks
    
    def calc_gamma_w_errors(self, Ms, Ls, fit_start, fit_stop, params):
        gamma_mean, gamma_blocks = self.calc_gamma_blocks(Ms, Ls, fit_start, fit_stop, params)
        return jk.get_errors_from_blocks(gamma_mean, gamma_blocks)
    
    # Calculating subtracted decay rate (should only work for 1D)-----
    def calc_gamma_sub(self, gamma1, t_full1, gamma2, t_full2):
        return (gamma2**0.5*t_full2-gamma1**0.5*t_full1)**2/(t_full2-t_full1)**2
    
    def calc_gamma_sub_blocks(self, gamma1_blocks, t_full1, gamma2_blocks, t_full2):
        return jk.super_jackknife_combine_blocks([gamma1_blocks, gamma2_blocks], lambda x: self.calc_gamma_sub(x[0], t_full1, x[1], t_full2))
    
    def calc_gamma_sub_w_errors(self, gamma1_blocks, t_full1, gamma2_blocks, t_full2):
        gamma_blocks = self.calc_gamma_sub_blocks(gamma1_blocks, t_full1, gamma2_blocks, t_full2)
        return jk.get_errors_from_blocks(np.mean(gamma_blocks), gamma_blocks)
    
    def calc_gamma_div(self, gamma1, gamma_der, t_full, dt):
        return (gamma_der*t_full*dt/gamma1**0.5)**2
    
    def calc_gamma_div_blocks(self, gamma1_blocks, gamma2_blocks, t_full, dt):
        return jk.super_jackknife_combine_blocks([gamma1_blocks, gamma2_blocks], lambda x: self.calc_gamma_div(x[0], x[1], t_full, dt))
    
    # Estimating systematic errors----------
    
    def calc_gamma_M_L_errors(self, Ms, Ls, fit_start, fit_stop, params):
        gammas = []
        for i in range(1,len(Ms)-1):
            M = Ms.pop(i)
            gammas.append(self.calc_gamma_blocks(Ms, Ls, fit_start, fit_stop, params))
            Ms.insert(i,M)
        for i in range(1,len(Ls)-1):
            L = Ls.pop(i)
            gammas.append(self.calc_gamma_blocks(Ms, Ls, fit_start, fit_stop, params))
            Ls.insert(i,L)
        return gammas
    
    # Retrieving data ============================================
    
    def get_exp_Ebar_blocks(self, sf, delta_t=1):
        t_TV = int(self.data.params[sf]["tTV"])
        t_FV = int(self.data.params[sf]["tFV"])
        blocks_TV = jk.get_jackknife_blocks(np.exp(self.data.delta_actions_t_TV[sf][f"{t_TV+delta_t}"][self.data.cutoff:]), self.data.block_size)
        blocks_FV = jk.get_jackknife_blocks(np.exp(self.data.delta_actions_t_FV[sf][f"{t_FV+delta_t}"][self.data.cutoff:]), self.data.block_size)
        return np.divide(blocks_FV,blocks_TV)
    
    # Functions for plotting data ================================
    
    def plot_exp_Ebar_blocks(self, param, delta_t=1, params={}, get_x=None, filter_x=lambda x: False, sort=None, label=None):
        sfs = self.data.get_indices(params)
        if(sort==None):
            sort = lambda x: float(self.data.params[x][param])
        if(get_x==None):
            get_x=sort
        sfs.sort(key=sort)
        expS = []
        expS_errs = []
        xs = []
        for sf in sfs:
            x=get_x(sf)
            if(filter_x(x)): continue
            blocks = self.get_exp_Ebar_blocks(sf, delta_t)
            dS, err = jk.get_errors_from_blocks(np.mean(blocks), blocks)
            expS.append(dS)
            expS_errs.append(err)
            xs.append(x)
        plt.errorbar(xs, expS, yerr=expS_errs, label=label)
        plt.title("exp(-$\\Delta S_{t_\\text{FV}\\to t_{\\text{FV}+1}})$")
        plt.xlabel(param)
        return xs, expS, expS_errs
    
    def plot_expS(self, delta_action, get_x=float, fact=1.0, label="p", filter_x=lambda x: False):
        expS = []
        expS_errs = []
        x = []
        for k in delta_action:
            if(filter_x(k)): continue
            x.append(get_x(k))
            blocks = jk.get_jackknife_blocks(np.exp(delta_action[k][self.data.cutoff:]), self.data.block_size)
            [eS, err] = jk.get_errors_from_blocks(np.mean(blocks), blocks)
            expS.append(eS*fact)
            expS_errs.append(err*fact)
        plt.errorbar(x, expS, yerr=expS_errs, label=label)
        print(x)
        print(expS)
        print(expS_errs)
        return x, expS, expS_errs
    
    def plot_expS_extend(self, delta_actions, sfs, param, get_x=float, filter_x=lambda x,y: False, sort=None):
        fact = 1.0
        last_params = []
        last_expS = []
        if(sort==None):
            sort = lambda x: float(self.data.params[x][param])
        sfs.sort(key=sort)
        for sf in sfs:
            print(sf)
            p = float(self.data.params[sf][param])
            if(p in last_params):
                fact = last_expS[last_params.index(p)]
            else:
                print(f"No previous factor found for {param}={p}")
            last_params, last_expS, errs = self.plot_expS(delta_actions[sf], get_x, fact, f"{param}={p}", lambda x: filter_x(x,p))
    
    def plot_expS_vs_M(self):
        sfs = list(filter(lambda x: self.data.params[x]["M"]!="1.0", list(self.data.delta_actions_M)))
        self.plot_expS_extend(self.data.delta_actions_M, sfs, "M")
    
    def plot_expS_vs_L(self):
        sfs = list(filter(lambda x: self.data.params[x]["L"]!="1.0", list(self.data.delta_actions_L)))
        self.plot_expS_extend(self.data.delta_actions_L, sfs, "L")
    
    def plot_expS_vs_t_TV(self, t_limit=[-100,100], sf=""):
        if(sf==""):
            sfs = list(self.data.delta_actions_t_TV)
        else:
            sfs = [sf]
        self.plot_expS_extend(self.data.delta_actions_t_TV, sfs, "tTV", get_x=int, filter_x=lambda x,x0: (int(x)-x0)<t_limit[0] or (int(x)-x0)>t_limit[1])
    
    def plot_expS_vs_t_FV(self, t_limit=[-100,100], sf=""):
        if(sf==""):
            sfs = list(self.data.delta_actions_t_FV)
        else:
            sfs = [sf]
        self.plot_expS_extend(self.data.delta_actions_t_FV, sfs, "tFV", get_x=int, filter_x=lambda x,x0: (int(x)-x0)<t_limit[0] or (int(x)-x0)>t_limit[1])
