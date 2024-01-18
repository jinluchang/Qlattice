import glob
import pickle
import matplotlib.pyplot as plt
import numpy as np
import jackknife as jk

class Data:
    def __init__(self, Nt, cutoff, block_size):
        self.Nt = Nt
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
    
    def get_sfs_list(self, sfs, params, values):
        rtn = []
        for sf in sfs:
            append = True
            for p in range(len(params)):
                if(self.get_param(sf,params[p])!=str(values[p])):
                    append = False
                    continue
            if(append):
                rtn.append(sf)
        return rtn
    
    def get_M_L_list(self, Ms, Ls):
        sf = list(self.delta_actions_M)[0]
        sfs_M = self.replace_params(sf, ["M", "L"], [[M, 1.0] for M in Ms])
        sfs_L = self.replace_params(sf, ["M", "L"], [[1.0, L] for L in Ls])
        delta_actions_M = [np.exp(self.delta_actions_M[sfs_M[i]][str(Ms[i+1])][self.cutoff:]) for i in range(len(Ms)-1)]
        delta_actions_L = [np.exp(self.delta_actions_L[sfs_L[i]][str(Ls[i+1])][self.cutoff:]) for i in range(len(Ls)-1)]
        return delta_actions_M + delta_actions_L
    
    def calc_ratio(self, delta_actions, N_Ms):
        ratio = 1.0
        for i in range(len(delta_actions)):
            if(i<N_Ms):
                ratio *= np.mean(delta_actions[i])
            else:
                ratio /= np.mean(delta_actions[i])
        return ratio
    
    def plot_mean_path(self):
        #x = np.arange(-5,5,0.1)
        #for t in range(0,self.Nt,int(self.Nt/20)):
        #    plt.plot([min(self.action.V(i,t)*self.Nt/20.0, 200.0) + t for i in x],x)
        for sf in self.timeslices:
            plt.plot(np.mean(self.timeslices[sf][self.cutoff:],axis=0))
    
    def plot_paths(self):
        #x = np.arange(-5,5,0.1)
        #for t in range(0,self.Nt,int(self.Nt/20)):
        #    plt.plot([min(self.action.V(i,t)*self.Nt/20.0, 200.0) + t for i in x],x)
        for sf in self.timeslices:
            i=0
            #plt.plot(self.timeslices[sf][0])
            #plt.show()
            for ts in self.timeslices[sf][:]:
                if (i+1)%100==0: plt.plot(ts)
                if (i+1)%10000==0: plt.show()
                i+=1
    
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
    
    def plot_all_expS(self, delta_actions, sfs, param, get_x=float, filter_x=lambda x,y: False, sort=None):
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
        self.plot_all_expS(self.delta_actions_M, sfs, "M")
    
    def plot_expS_vs_L(self):
        sfs = list(filter(lambda x: self.get_param(x,"L")!="1.0", list(self.delta_actions_L)))
        self.plot_all_expS(self.delta_actions_L, sfs, "L")
    
    def plot_expS_vs_t_TV(self, t_limit=[-100,100], sf=""):
        if(sf==""):
            sfs = list(self.delta_actions_t_TV)
        else:
            sfs = [sf]
        self.plot_all_expS(self.delta_actions_t_TV, sfs, "tTV", get_x=int, filter_x=lambda x,x0: (int(x)-x0)<t_limit[0] or (int(x)-x0)>t_limit[1])
    
    def plot_delta_actions_vs_t_FV(self, t_limit=[-100,100], sf=""):
        if(sf==""):
            sfs = list(self.delta_actions_t_FV)
        else:
            sfs = [sf]
        self.plot_all_expS(self.delta_actions_t_FV, sfs, "tFV", get_x=int, filter_x=lambda x,x0: (int(x)-x0)<t_limit[0] or (int(x)-x0)>t_limit[1])
    
    def get_Ebar_blocks(self, sf, delta_t=1):
        t_TV = self.get_t_TV(sf)
        t_FV = int(self.get_param(sf,"tFV"))
        dt = float(self.get_param(sf, "dt"))
        blocks_TV = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_TV[sf][f"{t_TV+delta_t}"][self.cutoff:]), self.block_size)
        blocks_FV = jk.get_jackknife_blocks(np.exp(self.delta_actions_t_FV[sf][f"{t_FV+delta_t}"][self.cutoff:]), self.block_size)
        bdiv = np.log(np.divide(blocks_FV,blocks_TV))/(dt*delta_t)
        return bdiv
    
    def get_Ebar_E_FV(self, sf, delta_t=1):
        bdiv = self.get_Ebar_blocks(sf, delta_t)
        return jk.get_errors_from_blocks(np.mean(bdiv), bdiv)
    
    def get_Ebar_slope(self, sf1, sf2, delta_t=1):
        dt_TV = self.get_t_TV(sf2)*float(self.get_param(sf2, "dt")) - self.get_t_TV(sf1)*float(self.get_param(sf2, "dt"))
        bdiv1 = self.get_Ebar_blocks(sf1, delta_t)
        bdiv2 = self.get_Ebar_blocks(sf2, delta_t)
        bdiv = (bdiv1 - bdiv2)**0.5 / dt_TV
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
        #print(bdiv)
        return jk.get_errors_from_blocks(np.mean(bdiv), bdiv)
    
    def plot_Ebar_E_FV(self, delta_t=1):
        sfs = self.get_sfs_list(list(self.delta_actions_t_TV), ["M", "L"], [1.0, 0.0])
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
    
    def plot_delta_E(self):
        sfs = self.get_sfs_list(list(self.delta_actions_t_TV), ["M", "L"], [1.0, 0.0])
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
    
    def plot_Ebar_slope(self, delta_t=1):
        sfs = self.get_sfs_list(list(self.delta_actions_t_TV), ["M", "L"], [1.0, 0.0])
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
    
    def plot_ratio1_vs_t(self, delta_action, dt, t_limit=[-100,100], fact=1.0, label="X"):
        delta_actions = []
        delta_action_errs = []
        t = []
        ks = list(delta_action)
        blocks_prev = jk.get_jackknife_blocks(np.exp(delta_action[ks[0]][self.cutoff:]), self.block_size)
        for k in ks[1:]:
            try:
                blocks = jk.get_jackknife_blocks(np.exp(delta_action[k][self.cutoff:]), self.block_size)
                bdiv = np.log(np.divide(blocks,blocks_prev))/dt
                #print(f"dS: {delta_actions[-1]}")
                #print(np.mean(blocks))
                #print(np.mean(blocks_prev))
                blocks_prev = blocks
                #print(f"{int(k)}, {t_limit}")
                if(int(k)<t_limit[0] or int(k)>t_limit[1]): continue
                [da, err] = jk.get_errors_from_blocks(np.mean(bdiv), bdiv)
                t.append(int(k))
                delta_actions.append(da*fact)
                delta_action_errs.append(err*fact)
            except ValueError:
                continue
        plt.errorbar(t, delta_actions, yerr=delta_action_errs, label=label)
        print(t)
        print(delta_actions)
    
    def plot_ratio1_vs_t(self, delta_action, dt, t_limit=[-100,100], fact=1.0, label="X"):
        delta_actions = []
        delta_action_errs = []
        t = []
        ks = list(delta_action)
        blocks_prev = jk.get_jackknife_blocks(np.exp(delta_action[ks[0]][self.cutoff:]), self.block_size)
        for k in ks[1:]:
            try:
                blocks = jk.get_jackknife_blocks(np.exp(delta_action[k][self.cutoff:]), self.block_size)
                bdiv = np.log(np.divide(blocks,blocks_prev))/dt
                #print(f"dS: {delta_actions[-1]}")
                #print(np.mean(blocks))
                #print(np.mean(blocks_prev))
                blocks_prev = blocks
                #print(f"{int(k)}, {t_limit}")
                if(int(k)<t_limit[0] or int(k)>t_limit[1]): continue
                [da, err] = jk.get_errors_from_blocks(np.mean(bdiv), bdiv)
                t.append(int(k))
                delta_actions.append(da*fact)
                delta_action_errs.append(err*fact)
            except ValueError:
                continue
        plt.errorbar(t, delta_actions, yerr=delta_action_errs, label=label)
        print(t)
        print(delta_actions)
    
    def plot_ratio1_vs_t_TV(self, t_limit=100, sf=""):
        if(sf==""):
            sfs = self.delta_actions_t_TV
        else:
            sfs = [sf]
        fact = 1.0
        for sf in sfs:
            print(sf)
            if(self.get_param(sf,"M")!="1.0" or self.get_param(sf,'L')!="0.0"):
                continue
            t_TV = self.get_t_TV(sf)
            label = f"tTV: {t_TV}, tFV: {self.get_param(sf,'tFV')}, tfull: {self.get_param(sf,'tfull')}"
            self.plot_ratio1_vs_t(self.delta_actions_t_TV[sf], float(self.get_param(sf,"dt")), [t_TV-t_limit, t_TV+t_limit], fact, label)
    
    def plot_ratio2_vs_t(self, delta_action, dt, t_limit=[-100,100], fact=1.0, label="X"):
        delta_actions = []
        delta_action_errs = []
        t = []
        ks = list(delta_action)
        blocks_prev2 = jk.get_jackknife_blocks(np.exp(delta_action[ks[0]][self.cutoff:]), self.block_size)
        blocks_prev = jk.get_jackknife_blocks(np.exp(delta_action[ks[1]][self.cutoff:]), self.block_size)
        for k in ks[2:]:
            try:
                blocks = jk.get_jackknife_blocks(np.exp(delta_action[k][self.cutoff:]), self.block_size)
                bdiv = np.log(np.divide(blocks,blocks_prev)*np.divide(blocks_prev2,blocks_prev))/dt/dt
                #print(f"dS: {delta_actions[-1]}")
                #print(np.mean(blocks))
                #print(np.mean(blocks_prev))
                blocks_prev2 = blocks_prev
                blocks_prev = blocks
                #print(f"{int(k)}, {t_limit}")
                if(int(k)<t_limit[0] or int(k)>t_limit[1]): continue
                [da, err] = jk.get_errors_from_blocks(np.mean(bdiv), bdiv)
                t.append(int(k))
                delta_actions.append(da*fact)
                delta_action_errs.append(err*fact)
            except ValueError:
                continue
        plt.errorbar(t, delta_actions, yerr=delta_action_errs, label=label)
        print(t)
        print(delta_actions)
    
    def plot_ratio2_vs_t_TV(self, t_limit=100, sf=""):
        if(sf==""):
            sfs = self.delta_actions_t_TV
        else:
            sfs = [sf]
        fact = 1.0
        for sf in sfs:
            print(sf)
            if(self.get_param(sf,"M")!="1.0" or self.get_param(sf,'L')!="0.0"):
                continue
            t_TV = self.get_t_TV(sf)
            label = f"tTV: {t_TV}, tFV: {self.get_param(sf,'tFV')}, tfull: {self.get_param(sf,'tfull')}"
            self.plot_ratio2_vs_t(self.delta_actions_t_TV[sf], float(self.get_param(sf,"dt")), [t_TV-t_limit, t_TV+t_limit], fact, label)
    
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
    
    def load(self, save_file):
        files = glob.glob(save_file)
        if len(files):
            for sf in files:
                print(f"Loading {sf}")
                with open(sf,"rb") as input:
                    data = pickle.load(input)
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