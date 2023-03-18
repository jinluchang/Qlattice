import numpy as np
import pickle
import glob
import datetime
import sys

class Correlators():
    def __init__(self, Nx, Nt, msq, lmbd, alpha, version, cutoff, block_size):
        self.psq_list=[]
        self.phi_list=[]
        self.timeslices=[]
        self.timeslices_m=[]
        self.hm_timeslices=[]
        self.A_timeslices=[]
        self.psq_pred_list=[]
        self.phi_pred_list=[]
        self.timeslices_pred=[]
        self.hm_timeslices_pred=[]
        self.A_timeslices_pred=[]
        self.loaded_files=[]
        self.data_len=0
        # Extent of the lattice
        self.Nx = Nx
        self. Nt = Nt
        # Spacial volume of the lattice
        self.Vx = self.Nx**3
        # The parameters in the action
        self.msq = msq
        self.lmbd = lmbd
        self.alpha = alpha
        # Version of the data generation code
        self.version = version
        # The number of trajectories until thermalization
        self.cutoff = cutoff
        # The size of blocks to use for error estimation
        self.block_size = block_size
        # Initial "done" variables
        self.reset_calcs()
    
    def reset_calcs(self):
        # These variables keep track of what's already been calculated
        self.done = []
        #
        self.vev={"sigma": 0, "psq": 0, "psqm": []}
        self.vev_err={"sigma": 0, "psq": 0, "psqm": []}
        self.corrs={}
        self.corr_avgs={}
    
    def load_data(self, date, day, filename=None):
        self.date = date
        self.day = day
        if(filename==None):
            filename = f"sigma_pion_corrs_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{date}-{day}_{self.version}.bin"
        print(f"Loading {filename}")
        self.loaded_files.append(filename)
        with open(filename,"rb") as input:
            data = pickle.load(input)
            self.psq_list.extend(data["psq_list"])
            self.phi_list.extend(data["phi_list"])
            self.timeslices.extend(data["timeslices"])
            self.timeslices_m.extend(data["timeslices_m"])
            self.hm_timeslices.extend(data["hm_timeslices"])
            self.A_timeslices.extend(data["ax_cur_timeslices"])
            self.psq_pred_list.extend(data["psq_pred_list"])
            self.phi_pred_list.extend(data["phi_pred_list"])
            self.timeslices_pred.extend(data["timeslices_pred"])
            self.hm_timeslices_pred.extend(data["hm_timeslices_pred"])
            self.A_timeslices_pred.extend(data["ax_cur_timeslices_pred"])
        self.data_len = len(self.phi_list)
        self.reset_calcs()
    
    def find_latest_date(self, files):
        if not len(files):
            return ""
        datenos=[]
        for f in files:
            info=f.split("_")
            for i in info:
                date = i.split("-")
                if(len(date)==3 and len(date[0])==4 and len(date[1])==2 and len(date[2])==2):
                    try:
                        datenos.append(float(date[0])+float(date[1])/12.0+float(date[2])/365.25)
                        break
                    except ValueError:
                        pass
        return datenos
    
    def load_latest_data(self):
        files = glob.glob(f"sigma_pion_corrs_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_*_{self.version}.bin")
        datenos = self.find_latest_date(files)
        self.load_data("-", "-", files[datenos.index(max(datenos))])
    
    def load_all_data(self):
        files = glob.glob(f"sigma_pion_corrs_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_*_{self.version}.bin")
        datenos = self.find_latest_date(files)
        files = sorted(files, key=lambda f: datenos[files.index(f)])
        for f in files:
            self.load_data("-", "-", f)
    
    def get_jackknife_blocks(self, data, f=lambda x:x):
        N = int(len(data)/self.block_size)
        data_mean = np.mean(data,axis=0)*N*self.block_size
        block_avgs = []
        for i in range(N):
            block_av=np.copy(data_mean)
            for j in range(self.block_size):
                block_av -= data[i*self.block_size+j]
            block_av /= (N-1)*self.block_size
            block_avgs.append(f(block_av))
        return block_avgs
    
    def get_errors_from_blocks(self, est_value, blocks):
        N = len(blocks)
        err = 0
        bias = 0
        for i in range(N):
            err = np.add(err, (N-1)/N*np.power(np.subtract(est_value,blocks[i]),2))
            bias = np.add(bias,np.divide(blocks[i],N))
        err = np.power(err,0.5)
        return [np.add(est_value,np.multiply(N-1,np.subtract(est_value,bias))), err]
    
    def apply_to_timeslices(self, ts, f):
        return [[f(ts[i][t]) for t in range(len(self.timeslices[0]))] 
            for i in range(self.cutoff,self.data_len)]
    
    def apply_to_timeslices_m(self, ts, f):
        return [[[f(ts[i][m][t]) for t in range(len(ts[0][0]))] 
            for m in range(len(ts[0]))]
            for i in range(self.cutoff,self.data_len)]
    
    def apply_to_obs(self, obs, f):
        return [f(obs[i]) for i in range(self.cutoff,self.data_len)]
    
    def apply_to_obs_m(self, obs, f):
        return [[f(obs[i][m]) for m in range(len(obs[0]))]
            for i in range(self.cutoff,self.data_len)]
    
    def calc_vev(self, name, values, m=""):
        if(f"{name}{m}" in self.done):
            return
        vev_blocks = self.get_jackknife_blocks(values)
        a,b = self.get_errors_from_blocks(np.mean(vev_blocks),vev_blocks)
        print(f"Vacuum Expectation Value of {name} = {a} +- {b}")
        if(m==""):
            self.vev[name], self.vev_err[name] = a,b
        else:
            self.vev[name][m], self.vev_err[name][m] = a,b
        self.done.append(f"{name}{m}")
    
    def calc_vev_m(self, name, values):
        M = len(values[0])
        self.vev[name] = [0.0]*M
        self.vev_err[name] = [0.0]*M
        values = np.array(values)
        for m in range(M):
            self.calc_vev(name, values[:,m], m)
    
    def calc_sigma_vev(self):
        self.calc_vev("sigma", self.apply_to_obs(self.phi_list, lambda x: x[0]))
    
    def calc_psq_vev(self):
        self.calc_vev("psq", self.apply_to_obs(np.array(self.timeslices),
            lambda ts: np.mean((ts[:,1]+ts[:,2]+ts[:,3])**2)/3.0/self.Vx/self.Nt))
    
    def calc_psqm_vev(self):
        self.calc_vev_m("psqm", self.apply_to_obs_m(np.array(self.timeslices_m),
            lambda ts: np.mean((ts[:,1]+ts[:,2]+ts[:,3])*np.conj(ts[:,1]+ts[:,2]+ts[:,3]))/3/self.Vx/self.Nt))
    
    def correlator(self,tslices1,tslices2,delta_t):
        rtn = 0
        for t0 in range(self.Nt):
            rtn += tslices1[t0%self.Nt]*tslices2[(t0+delta_t)%self.Nt]
        return rtn/self.Nt
    
    def calc_corrs(self, name, tslices1, tslices2, m=""):
        if(f"{name}{m}corrs" in self.done):
            return
        if(m==""):
            self.corrs[name]=[]
            for i in range(len(tslices1)):
                self.corrs[name].append([self.correlator(tslices1[i],tslices2[i],dt) for dt in range(self.Nt)])
            self.corr_avgs[name] = np.mean(self.corrs[name],axis=0)
        else:
            self.corrs[name][m]=[]
            try:
                for i in range(len(tslices1)):
                    self.corrs[name][m].append([self.correlator(tslices1[i][m],tslices2[i][m],dt) for dt in range(self.Nt)])
            except:
                for i in range(len(tslices1)):
                    self.corrs[name][m].append([self.correlator(tslices1[i][m],tslices2[i],dt) for dt in range(self.Nt)])
            self.corr_avgs[name][m] = np.mean(self.corrs[name][m],axis=0)
        self.done.append(f"{name}{m}corrs")
    
    def calc_corrs_m(self, name, tslices1, tslices2):
        M = len(tslices1[0])
        self.corrs[name] = [0.0]*M
        self.corr_avgs[name] = [0.0]*M
        for m in range(M):
            self.calc_corrs(name, tslices1, tslices2, m)
    
    def calc_pipi_corrs(self):
        pion_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: (ts[1]+ts[2]+ts[3])/3**0.5)
        self.calc_corrs("pipi", pion_ts, pion_ts)
    
    def calc_ss_corrs(self):
        self.calc_sigma_vev()
        sigma_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: ts[0] - self.vev["sigma"]*self.Vx)
        self.calc_corrs("ss", sigma_ts, sigma_ts)
    
    def calc_pipi_s_corrs(self):
        self.calc_sigma_vev()
        self.calc_psq_vev()
        pipi_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: (ts[1]+ts[2]+ts[3])**2/3 - self.vev["psq"]*self.Vx)
        sigma_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: ts[0] - self.vev["sigma"]*self.Vx)
        self.calc_corrs("pipi_s", pipi_ts, sigma_ts)
    
    def calc_pipim_corrs(self):
        pion_ts_m = self.apply_to_timeslices_m(self.timeslices_m,
            lambda ts: (ts[1]+ts[2]+ts[3])/3**0.5)
        self.calc_corrs_m("pipim", pion_ts_m, np.conj(pion_ts_m))
    
    def calc_ssm_corrs(self):
        sigma_ts_m = self.apply_to_timeslices_m(self.timeslices_m,
            lambda ts: ts[0])
        self.calc_corrs_m("ssm", sigma_ts_m, np.conj(sigma_ts_m))
    
    def calc_pipim_pipim_corrs(self):
        self.calc_psqm_vev()
        pipim_ts = self.apply_to_timeslices_m(np.array(self.timeslices_m),
            lambda ts: (ts[1]+ts[2]+ts[3])*np.conj(ts[1]+ts[2]+ts[3])/3)
        pipim_ts -= np.array(self.vev["psqm"])[None,:,None]*self.Vx
        self.calc_corrs_m("pipim_pipim", pipim_ts, np.conj(pipim_ts))
    
    def calc_pipim_s_corrs(self):
        self.calc_sigma_vev()
        self.calc_psqm_vev()
        pipim_ts = self.apply_to_timeslices_m(np.array(self.timeslices_m),
            lambda ts: (ts[1]+ts[2]+ts[3])*np.conj(ts[1]+ts[2]+ts[3])/3)
        pipim_ts -= np.array(self.vev["psqm"])[None,:,None]*self.Vx
        sigma_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: ts[0] - self.vev["sigma"]*self.Vx)
        self.calc_corrs_m("pipim_s", pipim_ts, sigma_ts)
    
    def save(self):
        date = datetime.datetime.now().date()
        with open(f"output_data/corrs/corrs_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{date}_corrs.bin", "wb") as file:
            pickle.dump([self.loaded_files,
                         self.done,
                         self.vev,
                         self.vev_err,
                         self.corrs,
                         self.corr_avgs], file)
    
    def load(self, date):
        with open(f"output_data/corrs/corrsmeasurements_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{date}_corrs.bin","rb") as file:
            data = pickle.load(file)
            for f in data[0]:
                self.load_data("-","-",f)
            [self.loaded_files,
             self.done,
             self.vev,
             self.vev_err,
             self.corrs,
             self.corr_avgs]=data

def main():
    # The lattice dimensions
    total_site = [4,4,4,8]
    # The multiplicity of the scalar field
    mult = 4
    # Use action for a Euclidean scalar field. The Lagrangian will be:
    # (1/2)*[sum i]|dphi_i|^2 + (1/2)*m_sq*[sum i]|phi_i|^2
    #     + (1/24)*lmbd*([sum i]|phi_i|^2)^2
    m_sq = -8.
    lmbd = 32.0
    alpha = 0.1
    #
    version = "1-9"
    #
    cutoff = 500
    block_size = 100

    for i in range(1,len(sys.argv)):
        try:
            if(sys.argv[i]=="-d"):
                a = sys.argv[i+1].split("x")
                total_site = [int(a[j]) for j in range(4)]
            elif(sys.argv[i]=="-m"):
                m_sq = float(sys.argv[i+1])
            elif(sys.argv[i]=="-l"):
                lmbd = float(sys.argv[i+1])
            elif(sys.argv[i]=="-a"):
                alpha = float(sys.argv[i+1])
            elif(sys.argv[i]=="-c"):
                cutoff = float(sys.argv[i+1])
            elif(sys.argv[i]=="-b"):
                block_size = float(sys.argv[i+1])
        except:
            raise Exception("Invalid arguments")
    
    corrs = Correlators(total_site[0], total_site[3], m_sq, lmbd, alpha, version, cutoff, block_size)
    corrs.load_all_data()
    print("pipi corrs==============================")
    corrs.calc_pipi_corrs()
    corrs.save()
    print("ss corrs==============================")
    corrs.calc_ss_corrs()
    corrs.save()
    print("pipim corrs==============================")
    corrs.calc_pipim_corrs()
    corrs.save()
    print("ssm corrs==============================")
    corrs.calc_ssm_corrs()
    corrs.save()
    print("pipim pipim corrs==============================")
    corrs.calc_pipim_pipim_corrs()
    corrs.save()  
    print("pipim s corrs==============================")
    corrs.calc_pipim_s_corrs()
    corrs.save()
