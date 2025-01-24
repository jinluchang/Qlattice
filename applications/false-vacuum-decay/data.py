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
        self.timeslices={}
        #
        self.delta_actions = {}
    
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
                    self.delta_actions[sf] = data["delta_actions"]
        
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
    
    # Functions for plotting data ================================
    
    def plot_mean_path(self, params={}, label="tTV"):
        sfs = self.get_indices(params)
        for sf in sfs:
            plt.plot(np.mean(self.timeslices[sf][int(self.cutoff/100):],axis=0), label=f"{label} = {self.params[sf][label]}")
    
    def plot_paths(self, params={}, sampling_freq=100, new_plot=10000):
        sfs = self.get_indices(params)
        for sf in sfs:
            i=0
            for ts in self.timeslices[sf][:]:
                if (i+1)%sampling_freq==0: plt.plot(ts)
                if (i+1)%new_plot==0: plt.show()
                i+=1

    def check_data(self, n_traj = 50000):
        #self.plot_mean_path()
        for i in self.accept_rates:
            if(np.mean(self.accept_rates[i]) < 0.7 or len(self.trajs[i])<n_traj):
                print(i)
                print(len(self.trajs[i]))
                print(np.mean(self.accept_rates[i][self.cutoff:]))
