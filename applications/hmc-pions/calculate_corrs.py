import numpy as np
import pickle
import glob
import datetime
import sys
import jackknife

class Correlators():
    def __init__(self, Nx, Nt, msq, lmbd, alpha, version, t_max, cutoff, block_size):
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
        # The maximum time extent for which to save the correlators
        self.t_max = t_max
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
            filename = f"output_data/sigma_pion_corrs_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{date}-{day}_{self.version}.bin"
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
        files = glob.glob(f"output_data/*_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_*_{self.version}.bin")
        datenos = self.find_latest_date(files)
        self.load_data("-", "-", files[datenos.index(max(datenos))])
    
    def load_all_data(self):
        files = glob.glob(f"output_data/*_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_*_{self.version}.bin")
        datenos = self.find_latest_date(files)
        files = sorted(files, key=lambda f: datenos[files.index(f)])
        for f in files:
            self.load_data("-", "-", f)
    
    def get_jackknife_blocks(self, data, f=lambda x:x):
        return jackknife.get_jackknife_blocks(data, self.block_size, f)
    
    def get_errors_from_blocks(self, est_value, blocks):
        return jackknife.get_errors_from_blocks(est_value, blocks)
    
    def apply_to_timeslices(self, ts, f):
        return [[f(ts[i][t]) for t in range(len(self.timeslices[0]))] 
            for i in range(self.cutoff,self.data_len)]
    
    def apply_to_timeslices_m(self, ts, f):
        return [[[f(ts[i][m][t]) for t in range(len(ts[0][0]))] 
            for m in range(len(ts[0]))]
            for i in range(self.cutoff,self.data_len)]

    def proj_s(self, timeslices_m):
        tslices_s = []
        for m_list in [[0,2],[3,5],[6,9],[10,12],[13,15]]:
            # Note the projection onto the s-wave is only correct for the lowest relative momentum
            tslices_s.append(np.mean(np.array(timeslices_m)[:,m_list[0]:m_list[1],:],axis=1))
        return np.swapaxes(np.array(tslices_s), 0, 1)
    
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
        if((not name in self.vev) or len(self.vev[name])==0):
            self.vev[name] = [0.0]*M
            self.vev_err[name] = [0.0]*M
        values = np.array(values)
        for m in range(M):
            self.calc_vev(name, values[:,m], m)
    
    def calc_sigma_vev(self):
        self.calc_vev("sigma", self.apply_to_obs(self.phi_list, self.sigma_interp_no_vevsub))
    
    def calc_psq_I0_vev(self):
        self.calc_vev("psq_I0", np.mean(self.apply_to_timeslices(np.array(self.timeslices),
            lambda ts: self.pipi_I0_interp_no_vevsub(ts)/self.Vx), axis=1))
    
    def calc_psqm_I0_vev(self):
        self.calc_vev_m("psqm_I0", np.mean(self.proj_s(self.apply_to_timeslices_m(np.array(self.timeslices_m),
            lambda ts: self.pipi_I0_interp_no_vevsub(ts)/self.Vx)), axis=2))
    
    def correlator(self,tslices1,tslices2,delta_t):
        rtn = 0
        for t0 in range(self.Nt):
            rtn += tslices1[t0%self.Nt]*tslices2[(t0+delta_t)%self.Nt]
        return rtn/self.Nt
    
    def calc_corrs(self, name, tslices1, tslices2, m=""):
        if(f"{name}{m}corrs" in self.done):
            return
        if(m==""):
            corrs=[]
            for i in range(len(tslices1)):
                #print(i)
                corrs.append([self.correlator(tslices1[i],tslices2[i],dt) for dt in range(self.t_max)])
            self.corr_avgs[name] = np.mean(corrs,axis=0)
            self.corrs[name] = corrs #self.get_jackknife_blocks(corrs)
        else:
            corrs=[]
            try:
                for i in range(len(tslices1)):
                    corrs.append([self.correlator(tslices1[i][m],tslices2[i][m],dt) for dt in range(self.t_max)])
            except:
                for i in range(len(tslices1)):
                    corrs.append([self.correlator(tslices1[i][m],tslices2[i],dt) for dt in range(self.t_max)])
            self.corr_avgs[name][m] = np.mean(corrs,axis=0)
            self.corrs[name][m] = corrs #self.get_jackknife_blocks(corrs)
        self.done.append(f"{name}{m}corrs")
    
    def calc_corrs_m(self, name, tslices1, tslices2):
        M = len(tslices1[0])
        if((not name in self.corrs) or len(self.corrs[name])==0):
            self.corrs[name] = [0.0]*M
            self.corr_avgs[name] = [0.0]*M
        for m in range(M):
            self.calc_corrs(name, tslices1, tslices2, m)
    
    def pion_interp(self, ts):
        return (ts[1]+ts[2]+ts[3])/3**0.5
    
    def sigma_interp(self, ts):
        return ts[0] - self.vev["sigma"]*self.Vx
    
    def sigma_interp_no_vevsub(self, ts):
        return ts[0]
    
    def pipi_I0_interp_no_vevsub(self, ts):
        return (-ts[1]*np.conj(ts[1])-ts[2]*np.conj(ts[2])-ts[3]*np.conj(ts[3]))/3**0.5
    
    def pipi_I0_interp(self, ts):
        return (-ts[1]*np.conj(ts[1])-ts[2]*np.conj(ts[2])-ts[3]*np.conj(ts[3]))/3**0.5 - self.vev["psq_I0"]*self.Vx
    
    def pipi_I2_interp(self, ts):
        return (-ts[1]*np.conj(ts[1])-ts[2]*np.conj(ts[2])+2*ts[3]*np.conj(ts[3]))/6**0.5
    
    def calc_pipi_corrs(self):
        pion_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: self.pion_interp(ts))
        self.calc_corrs("pipi", pion_ts, pion_ts)
    
    def calc_ss_corrs(self):
        self.calc_sigma_vev()
        sigma_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: self.sigma_interp(ts))
        self.calc_corrs("ss", sigma_ts, sigma_ts)
    
    def calc_pipi_s_corrs(self):
        self.calc_sigma_vev()
        self.calc_psq_I0_vev()
        pipi_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: self.pipi_I0_interp(ts))
        sigma_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: self.sigma_interp(ts))
        self.calc_corrs("pipi_s", pipi_ts, sigma_ts)
    
    def calc_pipi_pipi_I0_corrs(self):
        self.calc_psq_I0_vev()
        pipi_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: self.pipi_I0_interp(ts))
        self.calc_corrs("pipi_pipi_I0", pipi_ts, pipi_ts)
    
    def calc_pipi_pipi_I2_corrs(self):
        pipi_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: self.pipi_I2_interp(ts))
        self.calc_corrs("pipi_pipi_I2", pipi_ts, pipi_ts)
    
    def calc_pipim_corrs(self):
        pion_ts_m = self.apply_to_timeslices_m(self.timeslices_m,
            lambda ts: self.pion_interp(ts))
        self.calc_corrs_m("pipim", pion_ts_m, np.conj(pion_ts_m))
    
    def calc_ssm_corrs(self):
        sigma_ts_m = self.apply_to_timeslices_m(self.timeslices_m,
            lambda ts: self.sigma_interp_no_vevsub(ts))
        self.calc_corrs_m("ssm", sigma_ts_m, np.conj(sigma_ts_m))
    
    def calc_pipim_pipi_I0_corrs(self):
        self.calc_psq_I0_vev()
        self.calc_psqm_I0_vev()
        pipim_ts = self.proj_s(self.apply_to_timeslices_m(np.array(self.timeslices_m),
            lambda ts: self.pipi_I0_interp_no_vevsub(ts)))
        pipim_ts -= np.array(self.vev["psqm_I0"])[None,:,None]*self.Vx
        pipi_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: self.pipi_I0_interp(ts))
        self.calc_corrs_m("pipim_pipi_I0", pipim_ts, pipi_ts)
    
    def calc_pipim_pipi_I2_corrs(self):
        pipim_ts = self.proj_s(self.apply_to_timeslices_m(np.array(self.timeslices_m),
            lambda ts: self.pipi_I2_interp(ts)))
        pipi_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: self.pipi_I2_interp(ts))
        self.calc_corrs_m("pipim_pipi_I2", pipim_ts, pipi_ts)
    
    def calc_pipim_pipim_I0_corrs(self):
        self.calc_psqm_I0_vev()
        pipim_ts = self.proj_s(self.apply_to_timeslices_m(np.array(self.timeslices_m),
            lambda ts: self.pipi_I0_interp_no_vevsub(ts)))
        pipim_ts -= np.array(self.vev["psqm_I0"])[None,:,None]*self.Vx
        self.calc_corrs_m("pipim_pipim_I0", pipim_ts, np.conj(pipim_ts))
    
    def calc_pipim_pipim_I2_corrs(self):
        pipim_ts = self.proj_s(self.apply_to_timeslices_m(np.array(self.timeslices_m),
            lambda ts: self.pipi_I2_interp(ts)))
        self.calc_corrs_m("pipim_pipim_I2", pipim_ts, np.conj(pipim_ts))
    
    def calc_pipim_s_corrs(self):
        self.calc_sigma_vev()
        self.calc_psqm_I0_vev()
        pipim_ts = self.proj_s(self.apply_to_timeslices_m(np.array(self.timeslices_m),
            lambda ts: self.pipi_I0_interp_no_vevsub(ts)))
        pipim_ts -= np.array(self.vev["psqm_I0"])[None,:,None]*self.Vx
        sigma_ts = self.apply_to_timeslices(self.timeslices,
            lambda ts: self.sigma_interp(ts))
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
        with open(f"output_data/corrs/corrs_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{date}_corrs.bin","rb") as file:
            data = pickle.load(file)
            for f in data[0]:
                try:
                    self.load_data("-","-",f)
                except FileNotFoundError:
                    print("Warning: One of the files on which these correlators are based could not be found.")
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
    version = "3-1"
    #
    t_max = 32
    cutoff = 500
    block_size = 10

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
                cutoff = int(sys.argv[i+1])
            elif(sys.argv[i]=="-b"):
                block_size = int(sys.argv[i+1])
            elif(sys.argv[i]=="-t"):
                t_max = int(sys.argv[i+1])
        except:
            raise Exception("Invalid arguments")
    
    corrs = Correlators(total_site[0], total_site[3], m_sq, lmbd, alpha, version, t_max, cutoff, block_size)
    print("Loading data...")
    corrs.load_all_data()
    print("pipi corrs==============================")
    corrs.calc_pipi_corrs()
    corrs.save()
    print("ss corrs==============================")
    corrs.calc_ss_corrs()
    corrs.save()
    print("pipi pipi I0 corrs==============================")
    corrs.calc_pipi_pipi_I0_corrs()
    corrs.save()
    print("pipi pipi I2 corrs==============================")
    corrs.calc_pipi_pipi_I2_corrs()
    corrs.save()
    print("pipi s corrs==============================")
    corrs.calc_pipi_s_corrs()
    corrs.save()
    print("pipim corrs==============================")
    corrs.calc_pipim_corrs()
    corrs.save()
    print("ssm corrs==============================")
    corrs.calc_ssm_corrs()
    corrs.save()
    print("pipim pipi I0 corrs==============================")
    corrs.calc_pipim_pipi_I0_corrs()
    corrs.save()
    print("pipim pipi I2 corrs==============================")
    corrs.calc_pipim_pipi_I2_corrs()
    corrs.save()
    print("pipim pipim I0 corrs==============================")
    corrs.calc_pipim_pipim_I0_corrs()
    corrs.save()
    print("pipim pipim I2 corrs==============================")
    corrs.calc_pipim_pipim_I2_corrs()
    corrs.save()
    print("pipim s corrs==============================")
    corrs.calc_pipim_s_corrs()
    corrs.save()

if __name__ == '__main__':
    print("Hello")
    main()
