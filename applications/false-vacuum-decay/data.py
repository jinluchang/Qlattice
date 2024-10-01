import glob
import pickle
import matplotlib.pyplot as plt
import numpy as np
import jackknife as jk
from scipy.integrate import quad_vec
from scipy.optimize import curve_fit

import ratios_fit

class Data:
    def __init__(self, cutoff, block_size):
        self.cutoff = cutoff
        self.block_size = block_size
        # Stores a dictionary of parameters
        self.params = {}
        # Stores the trajectory number for debugging purposes
        self.trajs = {}
        # Save the acceptance rates
        self.accept_rates = {}
        # Stores the average phi^2 for each trajectory
        self.psq_list={}
        # Stores the average values of each phi_i for each trajectory
        self.phi_list={}
        #
        self.timeslices={}
        #
        self.fields={}
        self.momentums={}
        self.forces={}
        #
        self.delta_actions_M = {}
        self.actions_M = {}
        self.delta_actions_L = {}
        self.actions_L = {}
        self.delta_actions_t_FV = {}
        self.actions_t_FV = {}
        self.delta_actions_t_TV = {}
        self.actions_t_TV = {}
    
    def load(self, save_file):
        files = glob.glob(save_file)
        if len(files):
            for sf in files:
                with open(sf,"rb") as input:
                    print(f"Loading {sf} ...")
                    data = pickle.load(input)
                    sf = self.remove_date(sf)
                    if(sf in list(self.trajs)):
                        print(f"Already loaded ensemble with same parameters as {sf}")
                        continue
                    self.params[sf] = self.parse(sf)
                    self.trajs[sf] = data["trajs"]
                    self.accept_rates[sf] = data["accept_rates"]
                    self.psq_list[sf] = data["psq_list"]
                    self.phi_list[sf] = data["phi_list"]
                    self.timeslices[sf] = data["timeslices"]
                    self.fields[sf] = data["fields"]
                    self.momentums[sf] = data["momentums"]
                    self.forces[sf] = data["forces"]
                    self.delta_actions_M[sf] = data["delta_actions_M"]
                    self.delta_actions_L[sf] = data["delta_actions_L"]
                    self.delta_actions_t_FV[sf] = data["delta_actions_t_FV"]
                    self.delta_actions_t_TV[sf] = data["delta_actions_t_TV"]
                    print(f"# traj: {len(data['trajs'])}")
                    print(f"Accept rate: {np.mean(data['accept_rates'])}")
                    print(f"... Loaded {sf}")
    
    # Functions for working with file names ======================

    def parse(self, sf):
        params = {}
        p = False
        i = 0
        a = sf.split("_")
        while i < len(a):
            if p and len(a[i].split("-"))==1:
                params[a[i]] = a[i+1]
                i += 1
            elif len(a[i].split("x"))==2:
                params["Nt"] = int(a[i].split("x")[1])
                p = True
            i += 1
        return params
    
    def remove_date(self, sf):
        a = sf.split("_")
        for i in range(len(a)):
            if len(a[i].split("-"))==3:
                a[i] = "*"
        return "_".join(a)
    
    def replace_param(self, sf, param, val):
        a = sf.split("_")
        a[a.index(param)+1] = str(val)
        return "_".join(a)
    
    def replace_params(self, sf, params, values):
        rtn = []
        for v in values:
            s = sf
            for p in range(len(params)):
                s = self.replace_param(s, params[p], v[p])
            rtn.append(s)
        return rtn
    
    def get_indices(self, params):
        if(params=={}):
            return list(self.params)
        rtn = []
        for sf in self.params:
            append = True
            for p in params:
                if(str(params[p])!=str(self.params[sf][p])):
                    append = False
                    break
            if(append):
                rtn.append(sf)
        return rtn
    
    # Functions for calculating decay rate =======================
    
    # Calculating M-L ratios ---------------
    
    def get_M_L_blocks(self, Ms, Ls, params):
        sf_ML = self.get_indices(params)[0]
        sfs_M = self.replace_params(sf_ML, ["M", "L"], [[M, 1.0] for M in Ms])
        sfs_L = self.replace_params(sf_ML, ["M", "L"], [[1.0, L] for L in Ls])
        delta_actions_M = [jk.get_jackknife_blocks(np.exp(self.delta_actions_M[sfs_M[i]][str(Ms[i+1])][self.cutoff:]), self.block_size)
                           for i in range(len(Ms)-1)]
        delta_actions_L = [jk.get_jackknife_blocks(np.exp(self.delta_actions_L[sfs_L[i]][str(Ls[i+1])][self.cutoff:]), self.block_size)
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
    
    # Calculating fit ----------------------
    
    def fit_ratios(self, ratios, ratio_errs, t_TVs, fit_time, dt, start=0, stop=100, fitobject=ratios_fit.GaussianFit):
        fit = fitobject(dt)
        return fit.fit_correction(fit_time, t_TVs, ratios, ratio_errs, filter_x=lambda t: t<start or t>stop)
    
    def get_fit_ratios_blocks(self, params_tFV, t_TV, start=0, stop=100, fitobject=ratios_fit.GaussianFit):
        sfs = self.get_indices(params_tFV)
        sfs.sort(key=lambda x: float(self.params[x]["tFV"]))
        dt = float(self.params[sfs[0]]["dt"])
        
        dS_blocks = []
        errs = []
        t_TVs = []
        for sf in sfs:
            t = int(self.params[sf]["tTV"])*dt
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
    
    def calc_gamma_blocks(self, Ms, Ls, fit_start, fit_stop, params):
        params_tFV = params.copy()
        params_tFV["M"] = 1.0
        params_tFV["L"] = 0.0
        params_tFV.remove("tFV")
        #
        t_full = int(self.params[params]["tfull"])
        dt = float(self.params[params]["dt"])
        fit_time = int(self.params[params]["tTV"])*dt
        #
        ml_blocks = self.get_M_L_blocks(Ms, Ls, params)
        R_blocks = jk.super_jackknife_combine_blocks(ml_blocks, lambda x: self.calc_ratio(x, len(Ms)-1))
        fit_from_mean, fit_blocks = self.get_fit_ratios_blocks(params_tFV, fit_time, start=fit_start, stop=fit_stop)
        #
        print(f"Correction factor estimated: {fit_from_mean}")
        #
        gamma_blocks = jk.super_jackknife_combine_blocks([R_blocks, fit_blocks], lambda x: self.calc_gamma(x[0], x[1], t_full, dt))
        gamma_mean = self.calc_gamma(np.mean(R_blocks), fit_from_mean, t_full, dt)
        return gamma_mean, gamma_blocks
    
    def calc_gamma_w_errors(self, Ms, Ls, fit_start, fit_stop, params):
        gamma_mean, gamma_blocks = self.calc_gamma_blocks(Ms, Ls, fit_start, fit_stop, params)
        return jk.get_errors_from_blocks(gamma_mean, gamma_blocks)
    
    # Calculating subtracted decay rate-----
    def calc_gamma_sub(self, gamma1, t_full1, gamma2, t_full2):
        return (gamma2**0.5*t_full2-gamma1**0.5*t_full1)**2/(t_full2-t_full1)**2
    
    def calc_gamma_sub_blocks(self, gamma1_blocks, t_full1, gamma2_blocks, t_full2):
        return jk.super_jackknife_combine_blocks([gamma1_blocks, gamma2_blocks], lambda x: self.calc_gamma_sub(x[0], t_full1, x[1], t_full2))
    
    def calc_gamma_sub_w_errors(self, gamma1_blocks, t_full1, gamma2_blocks, t_full2):
        gamma_blocks = self.calc_gamma_sub_blocks(gamma1_blocks, t_full1, gamma2_blocks, t_full2)
        return jk.get_errors_from_blocks(np.mean(gamma_blocks), gamma_blocks)
    
    # Estimating systematic errors----------
    
    def calc_gamma_M_L_errors(self, Ms, Ls, fit_start, fit_stop, params):
        gammas = []
        for i in range(1,len(Ms)-1):
            M = Ms.pop(i)
            gammas.append(self.calc_gamma_blocks(Ms, Ls, fit_start, fit_stop, params)[0])
            Ms.insert(i,M)
        for i in range(1,len(Ls)-1):
            L = Ls.pop(i)
            gammas.append(self.calc_gamma_blocks(Ms, Ls, fit_start, fit_stop, params)[0])
            Ls.insert(i,L)
        return gammas
    
    # Retrieving data ============================================
    
    def get_exp_Ebar_blocks(self, sf, delta_t=1):
        t_TV = int(self.params[sf]["tTV"])
        t_FV = int(self.params[sf]["tFV"])
        blocks_TV = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_TV[sf][f"{t_TV+delta_t}"][self.cutoff:]), self.block_size)
        blocks_FV = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_FV[sf][f"{t_FV+delta_t}"][self.cutoff:]), self.block_size)
        return np.divide(blocks_FV,blocks_TV)
    
    # Functions for plotting data ================================
    
    def plot_mean_path(self, params={}):
        sfs = self.get_indices(params)
        for sf in sfs:
            plt.plot(np.mean(self.timeslices[sf][self.cutoff:],axis=0), label=sf)
    
    def plot_paths(self):
        for sf in self.timeslices:
            i=0
            #plt.plot(self.timeslices[sf][0])
            #plt.show()
            for ts in self.timeslices[sf][:]:
                if (i+1)%100==0: plt.plot(ts)
                if (i+1)%10000==0: plt.show()
                i+=1
    
    def plot_exp_Ebar_blocks(self, param, delta_t=1, sfs=None, get_x=None, filter_x=lambda x: False, sort=None):
        if(sfs==None):
            sfs = list(self.params)
        if(sort==None):
            sort = lambda x: float(self.params[x][param])
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
        plt.errorbar(xs, expS, yerr=expS_errs)
        return xs, expS, expS_errs
    
    def plot_expS(self, delta_action, get_x=float, fact=1.0, label="p", filter_x=lambda x: False):
        expS = []
        expS_errs = []
        x = []
        for k in delta_action:
            if(filter_x(k)): continue
            x.append(get_x(k))
            blocks = jk.get_jackknife_blocks(np.exp(delta_action[k][self.cutoff:]), self.block_size)
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
            sort = lambda x: float(self.params[x][param])
        sfs.sort(key=sort)
        for sf in sfs:
            print(sf)
            p = float(self.params[sf][param])
            if(p in last_params):
                fact = last_expS[last_params.index(p)]
            else:
                print(f"No previous factor found for {param}={p}")
            last_params, last_expS, errs = self.plot_expS(delta_actions[sf], get_x, fact, f"{param}={p}", lambda x: filter_x(x,p))
    
    def plot_expS_vs_M(self):
        sfs = list(filter(lambda x: self.params[x]["M"]!="1.0", list(self.delta_actions_M)))
        self.plot_expS_extend(self.delta_actions_M, sfs, "M")
    
    def plot_expS_vs_L(self):
        sfs = list(filter(lambda x: self.params[x]["L"]!="1.0", list(self.delta_actions_L)))
        self.plot_expS_extend(self.delta_actions_L, sfs, "L")
    
    def plot_expS_vs_t_TV(self, t_limit=[-100,100], sf=""):
        if(sf==""):
            sfs = list(self.delta_actions_t_TV)
        else:
            sfs = [sf]
        self.plot_expS_extend(self.delta_actions_t_TV, sfs, "tTV", get_x=int, filter_x=lambda x,x0: (int(x)-x0)<t_limit[0] or (int(x)-x0)>t_limit[1])
    
    def plot_expS_vs_t_FV(self, t_limit=[-100,100], sf=""):
        if(sf==""):
            sfs = list(self.delta_actions_t_FV)
        else:
            sfs = [sf]
        self.plot_expS_extend(self.delta_actions_t_FV, sfs, "tFV", get_x=int, filter_x=lambda x,x0: (int(x)-x0)<t_limit[0] or (int(x)-x0)>t_limit[1])
    
    # Obsolete ===================================================
    
    # Calculating decay rate using derivatives at a single point instead of a fit to all data
    
    def calc_gamma_Ebar_delta(self, R, Ebar, delta_E, t_full, dt):
        return R*(2*np.pi)**0.5/delta_E * np.exp(-Ebar**2/2/delta_E**2)/(t_full*dt)**2
    
    def calc_gamma_blocks_Ebar_delta(self, Ms, Ls, params_ML, params_tFV, der=1):
        t_full = int(params_ML["tfull"])
        dt = float(params_ML["dt"])
        t_FV = int(params_ML["tFV"])
        ml_blocks = self.get_M_L_blocks(Ms, Ls, params_ML)
        R_blocks = jk.super_jackknife_combine_blocks(ml_blocks, lambda x: self.calc_ratio(x, len(Ms)-1))
        params_ML["M"] = 1.0
        params_ML["L"] = 0.0
        Ebar_blocks = self.get_Ebar_blocks(self.get_indices(params_ML)[0])
        #
        sfs = self.get_indices(params_tFV)
        sfs.sort(key=lambda x: float(self.params[x]["tFV"]))
        if der==0:
            params_tFV["tFV"] = t_FV
            sf1 = self.get_indices(params_tFV)[0]
            sf2 = list(filter(lambda x: int(self.params[x]["tFV"])>t_FV, sfs))[0]
        elif der==1:
            sf1 = list(filter(lambda x: int(self.params[x]["tFV"])<t_FV, sfs))[-1]
            sf2 = list(filter(lambda x: int(self.params[x]["tFV"])>t_FV, sfs))[0]
        else:
            sf1 = list(filter(lambda x: int(self.params[x]["tFV"])<t_FV, sfs))[-1]
            params_tFV["tFV"] = t_FV
            sf2 = self.get_indices(params_tFV)[0]
        dE_blocks = self.get_Ebar_slope_blocks(sf1, sf2)
        print(f"Calculating dE with t_FV={self.params[sf1]['tFV']} and t_FV={self.params[sf2]['tFV']}")
        print(f"Correction factor estimated: {1 / (2*np.pi)**0.5 / np.mean(dE_blocks) * np.exp(-np.mean(Ebar_blocks)**2/2/np.mean(dE_blocks)**2)}")
        #
        gamma_blocks = jk.super_jackknife_combine_blocks([R_blocks, Ebar_blocks, dE_blocks], lambda x: self.calc_gamma_Ebar_delta(x[0], x[1], x[2], t_full, dt))
        gamma_mean = self.calc_gamma_Ebar_delta(np.mean(R_blocks), np.mean(Ebar_blocks), np.mean(dE_blocks), t_full, dt)
        return gamma_mean, gamma_blocks
    
    def calc_gamma_w_errors_Ebar_delta(self, Ms, Ls, params_ML, params_tFV, der=1):
        gamma_mean, gamma_blocks = self.calc_gamma_blocks_Ebar_delta(Ms, Ls, params_ML, params_tFV, der)
        return jk.get_errors_from_blocks(gamma_mean, gamma_blocks)
    
    def calc_gamma_dis_errors_Ebar_delta(self, Ms, Ls, params_ML, params_tFV):
        gamma = self.calc_gamma_blocks_Ebar_delta(Ms,Ls,params_ML, params_tFV, der=1)[0]
        fd = self.calc_gamma_blocks_Ebar_delta(Ms,Ls,params_ML, params_tFV, der=2)[0]
        bd = self.calc_gamma_blocks_Ebar_delta(Ms,Ls,params_ML, params_tFV, der=0)[0]
        lerr = min([fd,bd]) - gamma
        uerr = max([fd,bd]) - gamma
        if lerr>0 or uerr<0:
            lerr = max([abs(lerr),abs(uerr)])
            uerr = lerr
        return [abs(lerr), abs(uerr)]
    
    def get_Ebar_blocks(self, sf, delta_t=1):
        dt = float(self.params[sf]["dt"])
        bdiv = np.log(self.get_exp_Ebar_blocks(sf, delta_t))/(dt*delta_t)
        return bdiv
    
    def get_Ebar_E_FV(self, sf, delta_t=1):
        bdiv = self.get_Ebar_blocks(sf, delta_t)
        return jk.get_errors_from_blocks(np.mean(bdiv), bdiv)
    
    def get_Ebar_slope_blocks(self, sf1, sf2, delta_t=1):
        dt_TV = self.get_t_TV(sf2)*float(self.params[sf2]["dt"]) - self.get_t_TV(sf1)*float(self.params[sf2]["dt"])
        bdiv1 = self.get_Ebar_blocks(sf1, delta_t)
        bdiv2 = self.get_Ebar_blocks(sf2, delta_t)
        return ((bdiv1 - bdiv2) / dt_TV)**0.5
    
    def get_Ebar_slope(self, sf1, sf2, delta_t=1):
        bdiv = self.get_Ebar_slope_blocks(sf1, sf2, delta_t)
        return jk.get_errors_from_blocks(np.mean(bdiv), bdiv)
    
    def get_delta_E(self, sf):
        delta_t=1
        delta_t2=2
        t_TV = self.get_t_TV(sf)
        t_FV = int(self.params[sf]["tFV"])
        dt = float(self.params[sf]["dt"])
        
        blocks_TVa = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_TV[sf][f"{t_TV+delta_t}"][self.cutoff:]), self.block_size)
        blocks_FVa = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_FV[sf][f"{t_FV+delta_t}"][self.cutoff:]), self.block_size)
        
        blocks_TVa2 = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_TV[sf][f"{t_TV+delta_t2}"][self.cutoff:]), self.block_size)
        blocks_FVa2 = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_FV[sf][f"{t_FV+delta_t2}"][self.cutoff:]), self.block_size)
        
        bdiv = np.log(np.divide(blocks_FVa2,blocks_TVa2)/np.divide(blocks_FVa,blocks_TVa)**2.0)**0.5/(dt*delta_t2)
        return jk.get_errors_from_blocks(np.mean(bdiv), bdiv)
    
    # Plotting------------------------------
    
    def plot_Ebar_E_FV(self, params_tFV, delta_t=1):
        sfs = self.get_indices(params_tFV)
        sfs.sort(key=self.get_t_TV)
        t_TVs = []
        Es = []
        E_errs = []
        for sf in sfs:
            t_TVs.append(self.get_t_TV(sf)*float(self.params[sf]["dt"]))
            E, err = self.get_Ebar_E_FV(sf, delta_t)
            Es.append(E)
            E_errs.append(err)
        plt.errorbar(t_TVs, Es, yerr=E_errs)
    
    def plot_delta_E(self, params_tFV):
        sfs = self.get_indices(params_tFV)
        sfs.sort(key=self.get_t_TV)
        t_TVs = []
        dEs = []
        dE_errs = []
        for sf in sfs:
            t_TVs.append(self.get_t_TV(sf)*float(self.params[sf]["dt"]))
            dE, err = self.get_delta_E(sf)
            dEs.append(dE)
            dE_errs.append(err)
        plt.errorbar(t_TVs, dEs, yerr=dE_errs)
    
    def plot_Ebar_slope(self, params_tFV, delta_t=1):
        sfs = self.get_indices(params_tFV)
        sfs.sort(key=self.get_t_TV)
        t_TVs = []
        dEs = []
        dE_errs = []
        for sf in range(len(sfs)-1):
            t_TVs.append(self.get_t_TV(sfs[sf])*float(self.params[sfs[sf]]["dt"]))
            dE, err = self.get_Ebar_slope(sfs[sf], sfs[sf+1], delta_t)
            dEs.append(dE)
            dE_errs.append(err)
        plt.errorbar(t_TVs, dEs, yerr=dE_errs)
    
    def plot_change_over_mdtime(self, obs, block_size, t_limit=100):
        y = []
        y_err = []
        x = []
        for i in range(int(len(obs)/block_size)):
            x.append(i*block_size)
            blocks = jk.get_jackknife_blocks(obs[i*block_size:(i+1)*block_size], self.block_size)
            [meas, err] = jk.get_errors_from_blocks(np.mean(blocks), blocks)
            y.append(meas)
            y_err.append(err)
        plt.errorbar(x,y,yerr=y_err)