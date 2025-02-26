import glob
import pickle
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import jackknife as jk
from scipy.integrate import quad_vec
from scipy.optimize import curve_fit

import qlat as q

import ratios_fit

class Analysis:
    def __init__(self, data):
        self.data = data
    
    # Functions for calculating decay rate =======================
    
    # Calculating M-L ratios ---------------
    
    def get_M_L_blocks(self, Ms, Ls, params, der=False):
        sf_ML = self.data.get_indices(params)[0]
        sfs_M = self.data.replace_params(sf_ML, ["M", "L", "proj2"], [[M, 1.0, False] for M in Ms])
        sfs_L = self.data.replace_params(sf_ML, ["M", "L", "proj2"], [[1.0, L, False] for L in Ls])
        delta_actions_M = [jk.get_jackknife_blocks(np.exp(self.data.delta_actions[sfs_M[i]]["M"][str(Ms[i+1])][self.data.cutoff:]), self.data.block_size)
                           for i in range(len(Ms)-1)]
        delta_actions_L = [jk.get_jackknife_blocks(np.exp(self.data.delta_actions[sfs_L[i]]["L"][str(Ls[i+1])][self.data.cutoff:]), self.data.block_size)
                           for i in range(len(Ls)-1)]
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

    def get_P_blocks(self, Ps, params, der=False):
        try:
            sf_P = self.data.get_indices(params)[0]
        except IndexError as e:
            print(params)
            raise(e)
        sfs_P = self.data.replace_params(sf_P, ["P"], [[P] for P in Ps])
        delta_actions_P = [jk.get_jackknife_blocks(np.exp(self.data.delta_actions[sfs_P[i]]["P"][str(Ps[i+1])][self.data.cutoff:]), self.data.block_size)
                           for i in range(len(Ps)-1)]
        return delta_actions_P
    
    def calc_P_ratio(self, delta_actions):
        ratio = 1.0
        for i in range(len(delta_actions)):
            ratio *= np.mean(delta_actions[i])
        return ratio
        
    def get_ddR_blocks(self, Ms, Ls, Ps, params):
        R_blocks = self.get_R_blocks(Ms, Ls, params)
        P_blocks = self.get_P_blocks(Ps, params)
        ddR_over_R_blocks = jk.super_jackknife_combine_blocks(P_blocks, lambda x: self.calc_P_ratio(x))
        return jk.super_jackknife_combine_blocks([R_blocks, ddR_over_R_blocks], lambda x: x[0]*x[1])
    
    # Calculating fit ----------------------
    
    def fit_ratios(self, ratios, ratio_errs, t_TVs, fit_time, dt, start=0, stop=100, fitobject=ratios_fit.GaussianFit):
        fit = fitobject(dt)
        return fit.fit_correction(fit_time, t_TVs, ratios, ratio_errs, filter_x=lambda t: t<start or t>stop)
    
    def get_fit_ratios_blocks(self, params, start=0, stop=100, fitobject=ratios_fit.GaussianFit):
        t_TV = int(params["tTV"])*float(params["dt"])
        params_tTV = params.copy()
        params_tTV["M"] = 1.0
        params_tTV["L"] = 0.0
        params_tTV["P"] = 1.0
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
                dS_blocks.append(self.get_ddR_div_ddR_blocks(sf))
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
    
    #def get_R_div_R_blocks(self, sf, delta_t=1):
        #t_TV = int(self.data.params[sf]["tTV"])
        #t_FV = int(self.data.params[sf]["tFV"])
        #blocks_TV = jk.get_jackknife_blocks(np.exp(self.data.delta_actions_t_TV[sf][f"{t_TV+delta_t}"][self.data.cutoff:]), self.data.block_size)
        #blocks_FV = jk.get_jackknife_blocks(np.exp(self.data.delta_actions_t_FV[sf][f"{t_FV+delta_t}"][self.data.cutoff:]), self.data.block_size)
        #return np.divide(blocks_FV,blocks_TV)

    def get_ddR_div_ddR_blocks(self, sf):
        t_TV = int(self.data.params[sf]["tTV"])
        t_FV = int(self.data.params[sf]["tFV"])
        t_full2 = int(self.data.params[sf]["tfull2"])
        blocks_A = jk.get_jackknife_blocks(np.exp(self.data.delta_actions[sf]["DD"]["A"][self.data.cutoff:]), self.data.block_size)
        blocks_B = jk.get_jackknife_blocks(np.exp(self.data.delta_actions[sf]["DD"]["B"][self.data.cutoff:]), self.data.block_size)
        sf2 = self.data.replace_params(sf, ["disp", "tfull2", "tTV"], [["False", t_full2-1, t_TV+1]])[0]
        blocks_C = jk.get_jackknife_blocks(np.exp(self.data.delta_actions[sf2]["DD"]["C"][self.data.cutoff:]), self.data.block_size)
        blocks_D = jk.get_jackknife_blocks(np.exp(self.data.delta_actions[sf2]["DD"]["D"][self.data.cutoff:]), self.data.block_size)
        return np.multiply(blocks_A, blocks_C), np.multiply(blocks_B*blocks_D), np.multiply(blocks_A, blocks_C) / np.multiply(blocks_B*blocks_D)

    # Functions for plotting data ================================
    
    def plot_exp_Ebar_blocks(self, param, params={}, get_x=None, filter_x=lambda x: False, sort=None, label=None):
        sfs = self.data.get_indices(params)
        if(sort==None):
            sort = lambda x: float(self.data.params[x][param])
        if(get_x==None):
            get_x=sort
        sfs.sort(key=sort)
        expS1 = []
        expS_errs1 = []
        xs1 = []
        expS2 = []
        expS_errs2 = []
        xs2 = []
        expS3 = []
        expS_errs3 = []
        xs3 = []
        for sf in sfs:
            x=get_x(sf)
            if(filter_x(x)): continue
            blocks1, blocks2, blocks3 = self.get_ddR_div_ddR_blocks(sf)
            dS1, err1 = jk.get_errors_from_blocks(np.mean(blocks1), blocks1)
            dS2, err2 = jk.get_errors_from_blocks(np.mean(blocks2), blocks2)
            dS3, err3 = jk.get_errors_from_blocks(np.mean(blocks3), blocks3)
            expS1.append(dS1)
            expS_errs1.append(err1)
            xs1.append(x)
            expS2.append(dS2)
            expS_errs2.append(err2)
            xs2.append(x)
            expS3.append(dS3)
            expS_errs3.append(err3)
            xs3.append(x)
        #plt.errorbar(xs1, expS1, yerr=expS_errs1, label=label)
        #plt.errorbar(xs2, expS2, yerr=expS_errs2, label=label)
        plt.errorbar(xs3, expS3, yerr=expS_errs3, label=label)
        plt.title("exp(-$\\Delta S_{t_\\text{FV}\\to t_{\\text{FV}+1}})$")
        plt.xlabel(param)
        return xs1, expS1, expS_errs1
    
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
    
    def plot_expS_extend(self, delta_actions, key, sfs, param, get_x=float, filter_x=lambda x,y: False, sort=None):
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
            last_params, last_expS, errs = self.plot_expS(delta_actions[sf][key], get_x, fact, f"{param}={p}", lambda x: filter_x(x,p))
    
    def plot_expS_vs_M(self):
        sfs = list(filter(lambda x: self.data.params[x]["M"]!="1.0", list(self.data.delta_actions)))
        self.plot_expS_extend(self.data.delta_actions, "M", sfs, "M")
    
    def plot_expS_vs_L(self):
        sfs = list(filter(lambda x: self.data.params[x]["L"]!="1.0", list(self.data.delta_actions)))
        self.plot_expS_extend(self.data.delta_actions, "L", sfs, "L")
    
    def plot_expS_vs_P(self, params={"proj2": False}):
        sfs = list(filter(lambda x: self.data.params[x]["P"]!="1.0" and self.data.params[x]["M"]=="1.0" and self.data.params[x]["L"]=="0.0", list(self.data.get_indices(params))))
        self.plot_expS_extend(self.data.delta_actions, "P", sfs, "P")
    
    #def plot_expS_vs_t_TV(self, t_limit=[-100,100], sf=""):
        #if(sf==""):
        #    sfs = list(self.data.delta_actions_t_TV)
        #else:
        #    sfs = [sf]
        #self.plot_expS_extend(self.data.delta_actions_t_TV, sfs, "tTV", get_x=int, filter_x=lambda x,x0: (int(x)-x0)<t_limit[0] or (int(x)-x0)>t_limit[1])
    
    #def plot_expS_vs_t_FV(self, t_limit=[-100,100], sf=""):
        #if(sf==""):
        #    sfs = list(self.data.delta_actions_t_FV)
        #else:
        #    sfs = [sf]
        #self.plot_expS_extend(self.data.delta_actions_t_FV, sfs, "tFV", get_x=int, filter_x=lambda x,x0: (int(x)-x0)<t_limit[0] or (int(x)-x0)>t_limit[1])
    
    # Functions for plotting data ================================
    def plot_mean_path(self, params={}, label="tTV", ax=None, color="blue"):
        sfs = self.data.get_indices(params)
        for sf in sfs:
            if(ax==None):
                plt.plot(np.mean(self.data.timeslices[sf][self.data.cutoff:],axis=0), label=f"{label} = {self.data.params[sf][label]}", color=color)
            else:
                ax.plot(np.mean(self.data.timeslices[sf][self.data.cutoff:],axis=0), label=f"{label} = {self.data.params[sf][label]}", color=color)
    
    def plot_paths(self, params={}, sampling_freq=100, new_plot=10000, cutoff=0, ax=None, alpha=0.7, color="red"):
        sfs = self.data.get_indices(params)
        for sf in sfs:
            i=0
            for ts in self.data.timeslices[sf][cutoff:]:
                if(ax==None):
                    if (i+1)%sampling_freq==0: plt.plot(ts, alpha=alpha, color=color)
                    if (i+1)%new_plot==0: plt.show()
                else:
                    if (i+1)%sampling_freq==0: ax.plot(ts, alpha=alpha, color=color)
                i+=1

    def plot_potential(self, params, xmin=-1, xmax=2, fig=None, ax=None, vmin=-1, vmax=2, cmap="grey"):
        sf = self.data.get_indices(params)[0]
        action = q.QMAction(float(self.data.params[sf]["alpha"]), float(self.data.params[sf]["beta"]), float(self.data.params[sf]["FVoff"]), float(self.data.params[sf]["TVoff"]), float(self.data.params[sf]["bar"]), float(self.data.params[sf]["M"]), float(self.data.params[sf]["L"]), float(self.data.params[sf]["P"]), float(self.data.params[sf]["eps"]), int(self.data.params[sf]["tfull1"]), int(self.data.params[sf]["tfull2"]), int(self.data.params[sf]["tFVout"]), int(self.data.params[sf]["tFVmid"]), 0, float(self.data.params[sf]["dt"]), bool(self.data.params[sf]["Hlow"]), bool(self.data.params[sf]["disp"]))
        xs = np.arange(xmin,xmax,0.01)
        ts = np.arange(0, params["Nt"], 1)
        V_data = np.array([[action.V(x,t)-action.V(0,0) for t in ts[:-1]] for x in xs[:-1]])
        if(fig==None or ax==None):
            fig, ax = plt.subplots()
        pcm = ax.pcolormesh(ts, xs, V_data, cmap=matplotlib.colormaps[cmap], vmin=vmin, vmax=vmax)
        fig.colorbar(pcm, ax=ax)

    def check_data(self, n_traj = 50000):
        #self.plot_mean_path()
        for i in self.data.accept_rates:
            if(np.mean(self.data.accept_rates[i]) < 0.7 or len(self.data.trajs[i])<n_traj):
                print(i)
                print(len(self.data.trajs[i]))
                print(np.mean(self.data.accept_rates[i][self.data.cutoff:]))
