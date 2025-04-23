import glob
import pickle
import numpy as np

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
                        #self.trajs[sf].extend(data["trajs"])
                        #self.accept_rates[sf].extend(data["accept_rates"])
                        #self.psq_list[sf].extend(data["psq_list"])
                        #self.phi_list[sf].extend(data["phi_list"])
                        #self.timeslices[sf].extend(data["timeslices"])
                        #for key1 in self.delta_actions[sf]:
                        #    for key2 in self.delta_actions[sf][key1]:
                        #        self.delta_actions[sf][key1][key2].extend(data["delta_actions"][key1][key2])
                    else:
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
