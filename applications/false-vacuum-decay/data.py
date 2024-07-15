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
    
    def remove_date(self, sf):
        a = sf.split("_")
        for i in range(len(a)):
            if len(a[i].split("-"))==3:
                a[i] = "*"
        return "_".join(a)
    
    def get_t_TV(self, sf):
        return int(self.get_param(sf, "measurements").split("x")[1]) - 2*int(self.get_param(sf, "tfull")) - int(self.get_param(sf, "tFV"))
    
    def get_param(self, sf, param):
        if(param=="tTV"):
            return self.get_t_TV(sf)
        return sf.split(param)[1].split("_")[1]
    
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
    
    def get_sfs_list(self, sfs, profile):
        if(profile=="*"):
            return sfs
        rtn = []
        profile = profile.split("_")
        for sf in sfs:
            a = sf.split("_")
            append = True
            for i in range(len(a)):
                if(profile[i]=="*"):
                    continue
                elif(a[i]!=profile[i]):
                    append = False
                    break
            if(append):
                rtn.append(sf)
        return rtn
    
    # Functions for calculating decay rate =======================
    
    # Calculating M-L ratios ---------------
    
    def get_M_L_blocks(self, Ms, Ls, profile):
        sfs_M = self.replace_params(profile, ["M", "L"], [[M, 1.0] for M in Ms])
        sfs_L = self.replace_params(profile, ["M", "L"], [[1.0, L] for L in Ls])
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
    
    def get_fit_ratios_blocks(self, profile_tFV, t_TV, start=0, stop=100, fitobject=ratios_fit.GaussianFit):
        sfs = self.get_sfs_list(list(self.delta_actions_t_FV), profile_tFV)
        sfs.sort(key=lambda x: float(self.get_param(x, "tFV")))
        dt = float(self.get_param(sfs[0], "dt"))
        
        dS_blocks = []
        errs = []
        t_TVs = []
        for sf in sfs:
            t = self.get_t_TV(sf)*dt
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
    
    def calc_gamma_blocks(self, Ms, Ls, fit_start, fit_stop, profile):
        profile_tFV = self.replace_params(profile,["M","L","tFV"],[[1.0,0.0,"*"]])[0]
        #
        t_full = int(self.get_param(profile,"tfull"))
        dt = float(self.get_param(profile, "dt"))
        fit_time = int(self.get_param(profile,"tTV"))*dt
        #
        ml_blocks = self.get_M_L_blocks(Ms, Ls, profile)
        R_blocks = jk.super_jackknife_combine_blocks(ml_blocks, lambda x: self.calc_ratio(x, len(Ms)-1))
        fit_from_mean, fit_blocks = self.get_fit_ratios_blocks(profile_tFV, fit_time, start=fit_start, stop=fit_stop)
        #
        print(f"Correction factor estimated: {fit_from_mean}")
        #
        gamma_blocks = jk.super_jackknife_combine_blocks([R_blocks, fit_blocks], lambda x: self.calc_gamma(x[0], x[1], t_full, dt))
        gamma_mean = self.calc_gamma(np.mean(R_blocks), fit_from_mean, t_full, dt)
        return gamma_mean, gamma_blocks
    
    def calc_gamma_w_errors(self, Ms, Ls, fit_start, fit_stop, profile):
        gamma_mean, gamma_blocks = self.calc_gamma_blocks(Ms, Ls, fit_start, fit_stop, profile)
        return jk.get_errors_from_blocks(gamma_mean, gamma_blocks)
    
    # Estimating systematic errors----------
    
    def calc_gamma_M_L_errors(self, Ms, Ls, fit_start, fit_stop, profile):
        gammas = []
        for i in range(1,len(Ms)-1):
            M = Ms.pop(i)
            gammas.append(self.calc_gamma_blocks(Ms, Ls, fit_start, fit_stop, profile)[0])
            Ms.insert(i,M)
        for i in range(1,len(Ls)-1):
            L = Ls.pop(i)
            gammas.append(self.calc_gamma_blocks(Ms, Ls, fit_start, fit_stop, profile)[0])
            Ls.insert(i,L)
        return gammas
    
    # Retrieving data ============================================
    
    def get_exp_Ebar_blocks(self, sf, delta_t=1):
        t_TV = self.get_t_TV(sf)
        t_FV = int(self.get_param(sf,"tFV"))
        blocks_TV = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_TV[sf][f"{t_TV+delta_t}"][self.cutoff:]), self.block_size)
        blocks_FV = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_FV[sf][f"{t_FV+delta_t}"][self.cutoff:]), self.block_size)
        return np.divide(blocks_FV,blocks_TV)
    
    # Functions for plotting data ================================
    
    def plot_mean_path(self, profile="*"):
        sfs = self.get_sfs_list(list(self.timeslices), profile)
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
            sfs = list(self.delta_actions_t_TV)
        if(sort==None):
            sort = lambda x: float(self.get_param(x, param))
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
            sort = lambda x: float(self.get_param(x, param))
        sfs.sort(key=sort)
        for sf in sfs:
            print(sf)
            p = float(self.get_param(sf,param))
            if(p in last_params):
                fact = last_expS[last_params.index(p)]
            else:
                print(f"No previous factor found for {param}={p}")
            last_params, last_expS, errs = self.plot_expS(delta_actions[sf], get_x, fact, f"{param}={p}", lambda x: filter_x(x,p))
    
    def plot_expS_vs_M(self):
        sfs = list(filter(lambda x: self.get_param(x,"M")!="1.0", list(self.delta_actions_M)))
        self.plot_expS_extend(self.delta_actions_M, sfs, "M")
    
    def plot_expS_vs_L(self):
        sfs = list(filter(lambda x: self.get_param(x,"L")!="1.0", list(self.delta_actions_L)))
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
    
    def calc_gamma_blocks_Ebar_delta(self, Ms, Ls, profile_ML, profile_tFV, der=1):
        t_full = int(self.get_param(profile_ML,"tfull"))
        dt = float(self.get_param(profile_ML, "dt"))
        t_FV = int(self.get_param(profile_ML,"tFV"))
        ml_blocks = self.get_M_L_blocks(Ms, Ls, profile_ML)
        R_blocks = jk.super_jackknife_combine_blocks(ml_blocks, lambda x: self.calc_ratio(x, len(Ms)-1))
        Ebar_blocks = self.get_Ebar_blocks(self.replace_params(profile_ML,["M","L"],[[1.0,0.0]])[0])
        #
        sfs = self.get_sfs_list(list(self.delta_actions_t_FV), profile_tFV)
        sfs.sort(key=lambda x: float(self.get_param(x, "tFV")))
        if der==0:
            sf1 = self.replace_param(profile_tFV,"tFV",t_FV)
            sf2 = list(filter(lambda x: int(self.get_param(x,"tFV"))>t_FV, sfs))[0]
        elif der==1:
            sf1 = list(filter(lambda x: int(self.get_param(x,"tFV"))<t_FV, sfs))[-1]
            sf2 = list(filter(lambda x: int(self.get_param(x,"tFV"))>t_FV, sfs))[0]
        else:
            sf1 = list(filter(lambda x: int(self.get_param(x,"tFV"))<t_FV, sfs))[-1]
            sf2 = self.replace_param(profile_tFV,"tFV",t_FV)
        dE_blocks = self.get_Ebar_slope_blocks(sf1, sf2)
        print(f"Calculating dE with t_FV={self.get_param(sf1,'tFV')} and t_FV={self.get_param(sf2,'tFV')}")
        print(f"Correction factor estimated: {1 / (2*np.pi)**0.5 / np.mean(dE_blocks) * np.exp(-np.mean(Ebar_blocks)**2/2/np.mean(dE_blocks)**2)}")
        #
        gamma_blocks = jk.super_jackknife_combine_blocks([R_blocks, Ebar_blocks, dE_blocks], lambda x: self.calc_gamma_Ebar_delta(x[0], x[1], x[2], t_full, dt))
        gamma_mean = self.calc_gamma_Ebar_delta(np.mean(R_blocks), np.mean(Ebar_blocks), np.mean(dE_blocks), t_full, dt)
        return gamma_mean, gamma_blocks
    
    def calc_gamma_w_errors_Ebar_delta(self, Ms, Ls, profile_ML, profile_tFV, der=1):
        gamma_mean, gamma_blocks = self.calc_gamma_blocks_Ebar_delta(Ms, Ls, profile_ML, profile_tFV, der)
        return jk.get_errors_from_blocks(gamma_mean, gamma_blocks)
    
    def calc_gamma_dis_errors_Ebar_delta(self, Ms, Ls, profile_ML, profile_tFV):
        gamma = self.calc_gamma_blocks_Ebar_delta(Ms,Ls,profile_ML, profile_tFV, der=1)[0]
        fd = self.calc_gamma_blocks_Ebar_delta(Ms,Ls,profile_ML, profile_tFV, der=2)[0]
        bd = self.calc_gamma_blocks_Ebar_delta(Ms,Ls,profile_ML, profile_tFV, der=0)[0]
        lerr = min([fd,bd]) - gamma
        uerr = max([fd,bd]) - gamma
        if lerr>0 or uerr<0:
            lerr = max([abs(lerr),abs(uerr)])
            uerr = lerr
        return [abs(lerr), abs(uerr)]
    
    def get_Ebar_blocks(self, sf, delta_t=1):
        dt = float(self.get_param(sf, "dt"))
        bdiv = np.log(self.get_exp_Ebar_blocks(sf, delta_t))/(dt*delta_t)
        return bdiv
    
    def get_Ebar_E_FV(self, sf, delta_t=1):
        bdiv = self.get_Ebar_blocks(sf, delta_t)
        return jk.get_errors_from_blocks(np.mean(bdiv), bdiv)
    
    def get_Ebar_slope_blocks(self, sf1, sf2, delta_t=1):
        dt_TV = self.get_t_TV(sf2)*float(self.get_param(sf2, "dt")) - self.get_t_TV(sf1)*float(self.get_param(sf2, "dt"))
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
        t_FV = int(self.get_param(sf,"tFV"))
        dt = float(self.get_param(sf, "dt"))
        
        blocks_TVa = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_TV[sf][f"{t_TV+delta_t}"][self.cutoff:]), self.block_size)
        blocks_FVa = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_FV[sf][f"{t_FV+delta_t}"][self.cutoff:]), self.block_size)
        
        blocks_TVa2 = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_TV[sf][f"{t_TV+delta_t2}"][self.cutoff:]), self.block_size)
        blocks_FVa2 = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_FV[sf][f"{t_FV+delta_t2}"][self.cutoff:]), self.block_size)
        
        bdiv = np.log(np.divide(blocks_FVa2,blocks_TVa2)/np.divide(blocks_FVa,blocks_TVa)**2.0)**0.5/(dt*delta_t2)
        return jk.get_errors_from_blocks(np.mean(bdiv), bdiv)
    
    # Plotting------------------------------
    
    def plot_Ebar_E_FV(self, profile_tFV, delta_t=1):
        sfs = self.get_sfs_list(list(self.delta_actions_t_TV), profile_tFV)
        sfs.sort(key=self.get_t_TV)
        t_TVs = []
        Es = []
        E_errs = []
        for sf in sfs:
            t_TVs.append(self.get_t_TV(sf)*float(self.get_param(sf,"dt")))
            E, err = self.get_Ebar_E_FV(sf, delta_t)
            Es.append(E)
            E_errs.append(err)
        plt.errorbar(t_TVs, Es, yerr=E_errs)
    
    def plot_delta_E(self, profile_tFV):
        sfs = self.get_sfs_list(list(self.delta_actions_t_TV), profile_tFV)
        sfs.sort(key=self.get_t_TV)
        t_TVs = []
        dEs = []
        dE_errs = []
        for sf in sfs:
            t_TVs.append(self.get_t_TV(sf)*float(self.get_param(sf,"dt")))
            dE, err = self.get_delta_E(sf)
            dEs.append(dE)
            dE_errs.append(err)
        plt.errorbar(t_TVs, dEs, yerr=dE_errs)
    
    def plot_Ebar_slope(self, profile_tFV, delta_t=1):
        sfs = self.get_sfs_list(list(self.delta_actions_t_TV), profile_tFV)
        sfs.sort(key=self.get_t_TV)
        t_TVs = []
        dEs = []
        dE_errs = []
        for sf in range(len(sfs)-1):
            t_TVs.append(self.get_t_TV(sfs[sf])*float(self.get_param(sfs[sf],"dt")))
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