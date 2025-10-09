import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import pickle
import glob
import datetime

def autocorr(data):
    N = len(data)
    mean = np.mean(data)
    variance = np.var(data)
    data = np.subtract(data,mean)
    r = np.correlate(data, data, mode = 'full')[-N:]
    assert np.allclose(r, np.array([np.sum(np.multiply(data[:N-k],data[-(N-k):])) for k in range(N)]))
    result = r/(variance*(np.arange(N, 0, -1)))
    return result

class LatticeData():
    def __init__(self, Nx, Nt, msq, lmbd, alpha, version, cutoff, block_size):
        self.trajs=[]
        self.accept_rates=[]
        self.psq_list=[]
        self.phi_list=[]
        self.timeslices=[]
        self.timeslices_m=[]
        self.hm_timeslices=[]
        self.A_timeslices=[]
        self.phi_sq_dist=[]
        self.phi_i_dist=[]
        self.theta_dist=[]
        self.psq_dist_center = 0.0
        self.phi_dist_center = 0.0
        self.psq_pred_list=[]
        self.phi_pred_list=[]
        self.timeslices_pred=[]
        self.hm_timeslices_pred=[]
        self.A_timeslices_pred=[]
        self.fields=[]
        self.momentums=[]
        self.fields_pred=[]
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
        # The size of the blocks needed to get uncorrelated block averages
        self.block_size = block_size
        # The model to use for calculating mass
        self.cosh_model = lambda nt,A0,m : [A0*np.cosh((self.Nt/2.0-nt[i])*m)/np.cosh((self.Nt/2.0-1)*m) for i in range(len(nt))]
        # The model to use for calculating Fpi
        self.fpi_model = lambda nt,f_pi,m_pi : [-f_pi**2*m_pi*self.Vx*np.exp(-m_pi*self.Nt/2)*np.sinh((self.Nt/2-nt[i]+1/2)*m_pi)*(np.sinh((self.Nt/2-nt[i]-1/2)*m_pi)-np.sinh((self.Nt/2-nt[i]+0.5)*m_pi))/(self.alpha*np.cosh((self.Nt/2-nt[i])*m_pi)) for i in range(len(nt))]
        self.fpi2_model = lambda nt,f_pi,m_pi : [-f_pi/self.alpha*(m_pi*self.Vx*np.exp(-m_pi*self.Nt/2))**0.5*(np.sinh((self.Nt/2.0-nt[i]-1/2)*m_pi)-np.sinh((self.Nt/2.0-nt[i]+1/2)*m_pi))/(np.cosh((self.Nt/2.0-nt[i])*m_pi))**0.5 for i in range(len(nt))]
        # Initial "done" variables
        self.reset_calcs()
    
    def reset_calcs(self):
        # These variables keep track of what's already been calculated
        self.sigma_vev_done = False
        self.sigma_mass_eff_done = False
        self.pion_corrs_done = False
        self.sigma_corrs_done = False
        self.pion_mass_done = False
        self.sigma_mass_done = False
        self.A_corrs_done = False
        self.fpi_done = False
        self.fpi2_done = False
        self.scan_pion_mass_done = False
        self.scan_phi0_done = False
        self.pipi_sigma_corrs_done = False
        self.pipi_pipi_corrs_done = False
        self.phisq_done = False
        try:
            self.pipim_sigma_corrs_done = [False for m in range(len(self.timeslices_m[0]))]
            self.pipim_sigma_energy_done = [False for m in range(len(self.timeslices_m[0]))]
        except:
            self.pipim_sigma_corrs_done = [False]
            self.pipim_sigma_energy_done = [False]
        self.pipi_sigma_energy_done = False
        self.psqm_done = False
        #
        self.vev_sigma=0
        self.vev_sigma_err=0
        self.sigma_mass_eff=0
        self.sigma_mass_eff_err=0
        self.psq=0
        self.psq_err=0
        self.phisq_av=0
        self.phisq_av_err=0
        self.psqm_av=0
        self.psqm_av_err=0
        self.pi_corrs=0
        self.pi_corr_avgs=0
        self.s_corrs=0
        self.s_corr_avgs=0
        self.pipi_s_corrs=0
        self.pipi_s_corr_avgs=0
        self.pipi_pipi_corrs=0
        self.pipi_pipi_corr_avgs=0
        self.pipim_s_corrs=0
        self.pipim_s_corr_avgs=0
        self.pion_mass_est=0
        self.sigma_mass_est=0
        self.fpi_est=0
        self.pion_mass=0
        self.A_pion=0
        self.pion_mass_err=0
        self.A_pion_err=0
        self.sigma_mass=0
        self.A_sigma=0
        self.sigma_mass_err=0
        self.A_sigma_err=0
        self.pipim_sigma_E=0
        self.pipim_sigma_E_err=0
        self.pipi_sigma_E=0
        self.pipi_sigma_E_err=0
        self.A_corrs=0
        self.A_corr_avgs=0
        self.fpi=0
        self.ainv=0
        self.fpi2=0
        self.mpi_from_fpi=0
        self.mpi_from_fpi2=0
        self.fpi_err=0
        self.ainv_err=0
        self.fpi2_err=0
        self.mpi_from_fpi_err=0
        self.mpi_from_fpi2_err=0
        self.pion_masses=0
        self.pion_mass_errors=0
        self.k_pion_mass=0
        self.phi0s=0
        self.phi0_errors=0
        self.k_phi0=0
    
    def load_data(self, date, day, filename=None):
        self.date = date
        self.day = day
        if(filename==None):
            filename = f"output_data/measurements_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{date}-{day}_{self.version}.bin"
        print(f"Loading {filename}")
        self.loaded_files.append(filename)
        with open(filename,"rb") as input:
            data = pickle.load(input)
            self.trajs.extend(data["trajs"])
            self.accept_rates.extend(data["accept_rates"])
            self.psq_list.extend(data["psq_list"])
            self.phi_list.extend(data["phi_list"])
            self.timeslices.extend(data["timeslices"])
            self.timeslices_m.extend(data["timeslices_m"])
            self.hm_timeslices.extend(data["hm_timeslices"])
            self.A_timeslices.extend(data["ax_cur_timeslices"])
            self.phi_sq_dist.extend(data["phi_sq_dist"])
            self.phi_i_dist.extend(data["phi_i_dist"])
            self.theta_dist.extend(data["theta_dist"])
            self.psq_dist_center = data["psq_dist_center"]
            self.phi_dist_center = data["phi_dist_center"]
            self.psq_pred_list.extend(data["psq_pred_list"])
            self.phi_pred_list.extend(data["phi_pred_list"])
            self.timeslices_pred.extend(data["timeslices_pred"])
            self.hm_timeslices_pred.extend(data["hm_timeslices_pred"])
            self.A_timeslices_pred.extend(data["ax_cur_timeslices_pred"])
            self.fields.extend(data["fields"])
            self.momentums.extend(data["momentums"])
            self.fields_pred.extend(data["field_pred"])
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
    
    def cut_old_data(self):
        # If data from an earlier version was accidently included, this finds the most 
        # recent location where the trajectory count restarted
        try:
            start = next(i for i in reversed(range(len(self.trajs))) if self.trajs[i] == 2)
        except StopIteration:
            return
        self.trajs=self.trajs[start:]
        self.accept_rates=self.accept_rates[start:]
        self.psq_list=self.psq_list[start:]
        self.phi_list=self.phi_list[start:]
        self.timeslices=self.timeslices[start:]
        self.timeslices_m=self.timeslices_m[start:]
        self.hm_timeslices=self.hm_timeslices[start:]
        self.A_timeslices=self.A_timeslices[start:]
        self.phi_sq_dist=self.phi_sq_dist[start:]
        self.phi_i_dist=self.phi_i_dist[start:]
        self.theta_dist=self.theta_dist[start:]
        self.psq_pred_list=self.psq_pred_list[start:]
        self.phi_pred_list=self.phi_pred_list[start:]
        self.timeslices_pred=self.timeslices_pred[start:]
        self.hm_timeslices_pred=self.hm_timeslices_pred[start:]
        self.A_timeslices_pred=self.A_timeslices_pred[start:]
        self.fields=self.fields[start:]
        self.momentums=self.momentums[start:]
        self.fields_pred=self.fields_pred[start:]
        self.data_len = len(self.phi_list)
        self.reset_calcs()
    
    def plot_phi(self, zoom_win=300):
        # Plot an observable
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12,3))
        axes[0].plot([self.phi_list[i][1] for i in range(0,len(self.phi_list))])
        axes[1].plot([self.phi_list[i][1] for i in range(self.cutoff,self.cutoff+zoom_win)])
        fig.suptitle(f"Average phi_i \n {self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{self.date}-{self.day}_{self.version}")
        plt.tight_layout()
        print(f"Average acceptance rate after cutoff: {np.mean(self.accept_rates[self.cutoff:])}")
     
    def plot_sigma(self, zoom_win=300):
        # Plot an observable
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12,3))
        axes[0].plot([self.phi_list[i][0] for i in range(0,len(self.phi_list))])
        axes[1].plot([self.phi_list[i][0] for i in range(self.cutoff,self.cutoff+zoom_win)])
        fig.suptitle(f"Average sigma \n {self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{self.date}-{self.day}_{self.version}")
        plt.tight_layout()
        print(f"Average acceptance rate after cutoff: {np.mean(self.accept_rates[self.cutoff:])}")
        
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
    
    def calc_sigma_vev(self):
        if(self.sigma_vev_done):
            return
        vev_sigma_blocks = self.get_jackknife_blocks([self.phi_list[i][0] for i in range(self.cutoff,self.data_len)])
        self.vev_sigma, self.vev_sigma_err = self.get_errors_from_blocks(np.mean(vev_sigma_blocks),vev_sigma_blocks)
        print(f"Vacuum Expectation Value of sigma = {self.vev_sigma} +- {self.vev_sigma_err}")
        self.sigma_vev_done = True
    
    def calc_phisq(self):
        psq_blocks = self.get_jackknife_blocks([self.psq_list[i] for i in range(self.cutoff,self.data_len)])
        self.psq, self.psq_err = self.get_errors_from_blocks(np.mean(psq_blocks),psq_blocks)
        print(f"Average phi^2 = {self.psq} +- {self.psq_err}")
        self.sigma_vev_done = True
    
    def calc_phisq_av(self):
        # This is the phisq appropriate for correlation function calculations involving pi*pi
        if(self.phisq_done):
            return
        corrs = [self.pion_correlator(self.timeslices[i],0) for i in range(self.cutoff, self.data_len)]
        phisq_blocks = self.get_jackknife_blocks(corrs)
        self.phisq_av, self.phisq_av_err = self.get_errors_from_blocks(np.mean(corrs),phisq_blocks)
        print(f"Vacuum Expectation Value of phi_i^2 = {self.phisq_av} +- ?")
        self.phisq_done = True
    
    def calc_sigma_mass_eff(self):
        if(self.sigma_mass_eff_done):
            return
        self.calc_sigma_corrs()
        sigma_mass_eff_blocks = self.get_jackknife_blocks(self.s_corrs, lambda x : -np.log(x[1]/x[0]))
        self.sigma_mass_eff, self.sigma_mass_eff_err = self.get_errors_from_blocks(np.mean(sigma_mass_eff_blocks),sigma_mass_eff_blocks)
        print(f"Effective mass of sigma = {self.sigma_mass_eff} +- {self.sigma_mass_eff_err}")
        self.sigma_mass_eff_done = True
    
    def calc_psqm_av(self):
        if(self.psqm_done):
            return
        self.psqm_av = []
        self.psqm_av_err = []
        for m in range(len(self.timeslices_m[0])):
            corrs = [self.pipim_w_vev_correlator(self.timeslices_m[i][m],0) for i in range(self.cutoff, self.data_len)]
            psqm_blocks = self.get_jackknife_blocks(corrs)
            a, b = self.get_errors_from_blocks(np.mean(corrs),psqm_blocks)
            self.psqm_av.append(a)
            self.psqm_av_err.append(b)
            print(f"Vacuum Expectation Value of pipim = {self.psqm_av[m]} +- {self.psqm_av_err[m]}")
        self.psqm_done = True
    
    def pion_correlator(self,tslices,delta_t):
        rtn = 0
        for t0 in range(self.Nt):
            rtn += (tslices[t0%self.Nt][1]*tslices[(t0+delta_t)%self.Nt][1] + 
                    tslices[t0%self.Nt][2]*tslices[(t0+delta_t)%self.Nt][2] + 
                    tslices[t0%self.Nt][3]*tslices[(t0+delta_t)%self.Nt][3])/3.0
        return rtn/self.Nt
    
    def sigma_correlator(self,tslices,vev_sigma,delta_t):
        rtn = 0
        for t0 in range(self.Nt):
            rtn += (tslices[t0%self.Nt][0]/self.Vx-vev_sigma)*(tslices[(t0+delta_t)%self.Nt][0]/self.Vx-vev_sigma)
        return rtn/self.Nt
    
    def pipim_w_vev_correlator(self,tslices_m,delta_t):
        rtn = 0
        for t0 in range(self.Nt):
            rtn += (tslices_m[t0%self.Nt][1]*np.conj(tslices_m[(t0+delta_t)%self.Nt][1]) + 
                    tslices_m[t0%self.Nt][2]*np.conj(tslices_m[(t0+delta_t)%self.Nt][2]) + 
                    tslices_m[t0%self.Nt][3]*np.conj(tslices_m[(t0+delta_t)%self.Nt][3]))/3.0
        return rtn/self.Nt
    
    def pipi_sigma_correlator(self,tslices,vev_sigma,phisq_av,delta_t):
        rtn = 0
        for t0 in range(self.Nt):
            rtn += ((tslices[t0%self.Nt][1]*tslices[t0%self.Nt][1] + 
                    tslices[t0%self.Nt][2]*tslices[t0%self.Nt][2] + 
                    tslices[t0%self.Nt][3]*tslices[t0%self.Nt][3])/3.0-phisq_av)*(tslices[(t0+delta_t)%self.Nt][0]-vev_sigma*self.Vx)
        return rtn/self.Nt
    
    def pipi_pipi_correlator(self,tslices,phisq_av,delta_t):
        rtn = 0
        for t0 in range(self.Nt):
            rtn += ((tslices[t0%self.Nt][1]*tslices[t0%self.Nt][1] + 
                    tslices[t0%self.Nt][2]*tslices[t0%self.Nt][2] + 
                    tslices[t0%self.Nt][3]*tslices[t0%self.Nt][3])/3.0-phisq_av)*(
                    (tslices[(t0+delta_t)%self.Nt][1]*tslices[(t0+delta_t)%self.Nt][1] + 
                    tslices[(t0+delta_t)%self.Nt][2]*tslices[(t0+delta_t)%self.Nt][2] + 
                    tslices[(t0+delta_t)%self.Nt][3]*tslices[(t0+delta_t)%self.Nt][3])/3.0-phisq_av)
        return rtn/self.Nt
    
    def pipim_sigma_correlator(self,tslices,tslices_m,vev_sigma,psqm_av,delta_t):
        rtn = 0
        for t0 in range(self.Nt):
            rtn += ((tslices_m[t0%self.Nt][1]*np.conj(tslices_m[t0%self.Nt][1]) + 
                    tslices_m[t0%self.Nt][2]*np.conj(tslices_m[t0%self.Nt][2]) + 
                    tslices_m[t0%self.Nt][3]*np.conj(tslices_m[t0%self.Nt][3]))/3.0-psqm_av)*(tslices[(t0+delta_t)%self.Nt][0]-vev_sigma*self.Vx)
        return rtn/self.Nt
    
    def calc_pion_corrs(self):
        if(self.pion_corrs_done):
            return
        self.pi_corrs = [[self.pion_correlator(self.timeslices[i],dt) for dt in range(self.Nt)] for i in range(self.cutoff,self.data_len)]
        self.pi_corr_avgs = np.mean(self.pi_corrs,axis=0)
        self.pion_corrs_done = True
    
    def calc_sigma_corrs(self):
        if(self.sigma_corrs_done):
            return
        self.calc_sigma_vev()
        self.s_corrs = [[self.sigma_correlator(self.timeslices[i],self.vev_sigma,dt) for dt in range(self.Nt)] for i in range(self.cutoff,self.data_len)]
        self.s_corr_avgs = np.mean(self.s_corrs,axis=0)
        self.sigma_corrs_done = True
    
    def calc_pipi_sigma_corrs(self):
        if(self.pipi_sigma_corrs_done):
            return
        self.calc_sigma_vev()
        self.calc_phisq_av()
        self.pipi_s_corrs = [[self.pipi_sigma_correlator(self.timeslices[i],self.vev_sigma,self.phisq_av,dt) for dt in range(self.Nt)] for i in range(self.cutoff,self.data_len)]
        self.pipi_s_corr_avgs = np.mean(self.pipi_s_corrs,axis=0)
        self.pipi_sigma_corrs_done = True
    
    def calc_pipi_pipi_corrs(self):
        if(self.pipi_pipi_corrs_done):
            return
        self.calc_phisq()
        self.pipi_pipi_corrs = [[self.pipi_pipi_correlator(self.timeslices[i],self.phisq_av,dt) for dt in range(self.Nt)] for i in range(self.cutoff,self.data_len)]
        self.pipi_pipi_corr_avgs = np.mean(self.pipi_pipi_corrs,axis=0)
        self.pipi_pipi_corrs_done = True
    
    def calc_pipim_sigma_corrs(self,m):
        if(self.pipim_sigma_corrs_done[m]):
            return
        if(not any(self.pipim_sigma_corrs_done)):
            self.pipim_s_corrs = [[] for m in range(len(self.timeslices_m[0]))]
            self.pipim_s_corr_avgs = [[] for m in range(len(self.timeslices_m[0]))]
        self.calc_sigma_vev()
        self.calc_psqm_av()
        self.pipim_s_corrs[m] = []
        for i in range(self.cutoff,self.data_len):
            self.pipim_s_corrs[m].append([self.pipim_sigma_correlator(self.timeslices[i],self.timeslices_m[i][m],self.vev_sigma,self.psqm_av[m],dt) for dt in range(self.Nt)])
            print(f"{m} - {i}/{self.data_len}")
        self.pipim_s_corr_avgs[m] = np.mean(self.pipim_s_corrs[m],axis=0)
        self.pipim_sigma_corrs_done[m] = True
    
    def calc_quick_ests(self):
        self.calc_pion_corrs()
        self.calc_sigma_corrs()
        self.pion_mass_est = -np.log(self.pi_corr_avgs[1]/self.pi_corr_avgs[0])
        self.sigma_mass_est = -np.log(self.s_corr_avgs[1]/self.s_corr_avgs[0])
        self.fpi_est = 2**0.5*self.alpha*self.pi_corr_avgs[0]**0.5/((-np.log(self.pi_corr_avgs[1]/self.pi_corr_avgs[0]))**(3/2)*self.Vx**.5)
        print(f"The pion mass is approximately: {self.pion_mass_est}")
        print(f"The sigma mass is approximately: {self.sigma_mass_est}")
        print(f"F_pi is approximately: {self.fpi_est}")
    
    def find_mass_from_fit(self, corr_avgs, fit_range=[], mass_guess=None):
        if(len(fit_range)==0):
            fit_range=[0,int(self.Nt)]
        nt = range(fit_range[0],fit_range[1])
        A0_guess = corr_avgs[1]
        if(mass_guess==None):
            mass_guess = -np.log(corr_avgs[1]/corr_avgs[0])
        a = corr_avgs[fit_range[0]:fit_range[1]]
        pi_opt, pi_cov = curve_fit(self.cosh_model, nt, a, p0=[A0_guess,mass_guess])
        return [np.abs(pi_opt[1]), pi_opt[0]]
    
    def calc_pion_mass(self):
        if(self.pion_mass_done):
            return
        self.calc_pion_corrs()
        pion_mass_blocks = self.get_jackknife_blocks(self.pi_corrs, self.find_mass_from_fit)
        [[self.pion_mass, self.A_pion],[self.pion_mass_err,self.A_pion_err]] = self.get_errors_from_blocks(self.find_mass_from_fit(self.pi_corr_avgs), pion_mass_blocks)
        print(f"Pion mass is {self.pion_mass}/a +- {self.pion_mass_err}/a")
        self.pion_mass_done=True
        
    def calc_sigma_mass(self):
        if(self.sigma_mass_done):
            return
        self.calc_sigma_corrs()
        sigma_mass_blocks = self.get_jackknife_blocks(self.s_corrs, self.find_mass_from_fit)
        [[self.sigma_mass, self.A_sigma],[self.sigma_mass_err,self.A_sigma_err]] = self.get_errors_from_blocks(self.find_mass_from_fit(self.s_corr_avgs), sigma_mass_blocks)
        print(f"Sigma mass is {self.sigma_mass}/a +- {self.sigma_mass_err}/a")
        self.sigma_mass_done=True
    
    def calc_pipi_sigma_energy(self):
        if(self.pipi_sigma_energy_done):
            return
        self.calc_pipi_sigma_corrs()
        blocks = self.get_jackknife_blocks(np.real(self.pipi_s_corrs), self.find_mass_from_fit)
        [[E, A],[E_err,A_err]] = self.get_errors_from_blocks(self.find_mass_from_fit(np.real(self.pipi_s_corr_avgs)), blocks)
        print(f"<pipi sigma> energy is {E}/a +- {E_err}/a")
        self.pipi_sigma_E = E
        self.pipi_sigma_E_err = E_err
        self.pipi_sigma_energy_done=True
    
    def calc_pipim_sigma_energy(self,m):
        if(self.pipim_sigma_energy_done[m]):
            return
        if(not any(self.pipim_sigma_energy_done)):
            self.pipim_sigma_E = [[] for m in range(len(self.timeslices_m[0]))]
            self.pipim_sigma_E_err = [[] for m in range(len(self.timeslices_m[0]))]
        self.calc_pipim_sigma_corrs(m)
        blocks = self.get_jackknife_blocks(np.real(self.pipim_s_corrs[m]), self.find_mass_from_fit)
        [[E, A],[E_err,A_err]] = self.get_errors_from_blocks(self.find_mass_from_fit(np.real(self.pipim_s_corr_avgs[m])), blocks)
        print(f"<pipim sigma> energy is {E}/a +- {E_err}/a")
        self.pipim_sigma_E[m] = E
        self.pipim_sigma_E_err[m] = E_err
        self.pipim_sigma_energy_done[m]=True
    
    def axpi_correlator(self,tslices,Atslices,delta_t):
        rtn = 0
        for t0 in range(self.Nt):
            rtn += (Atslices[(t0+delta_t)%self.Nt][0]*tslices[t0%self.Nt][1] + 
                    Atslices[(t0+delta_t)%self.Nt][1]*tslices[t0%self.Nt][2] + 
                    Atslices[(t0+delta_t)%self.Nt][2]*tslices[t0%self.Nt][3])/3.0
        return rtn/self.Nt
    
    def calc_axpi_corrs(self):
        if(self.A_corrs_done):
            return
        #self.A_corrs = [[self.axpi_correlator(self.timeslices[i],self.A_timeslices[i],dt) for dt in range(self.Nt)] for i in range(self.cutoff,self.data_len)]
        self.A_corrs = []
        for i in range(self.cutoff,self.data_len):
            self.A_corrs.append([self.axpi_correlator(self.timeslices[i],self.A_timeslices[i],dt) for dt in range(self.Nt)])
            print(f"{i} / {self.data_len}")
        self.A_corr_avgs = np.mean(self.A_corrs,axis=0)
        self.A_corrs_done = True
    
    def fpi_fit(self, A_corr_avgs, fit_range=[], fpi_guess=0.2, mpi_guess=0.2):
        if(len(fit_range)==0):
            fit_range=[1,int(self.Nt)]
        nt = range(fit_range[0],fit_range[1])
        pi_opt, pi_cov = curve_fit(lambda nt, fpi, mpi : self.fpi_model(nt,fpi,mpi), nt, 
                                   A_corr_avgs[fit_range[0]:fit_range[1]],
                                   p0=[fpi_guess, mpi_guess])
        return pi_opt
    
    def calc_fpi(self):
        if(self.fpi_done):
            return
        self.calc_axpi_corrs()
        fpi_blocks = self.get_jackknife_blocks(self.A_corrs, self.fpi_fit)
        [[self.fpi, self.mpi_from_fpi],[self.fpi_err,self.mpi_from_fpi_err]] = self.get_errors_from_blocks(self.fpi_fit(self.A_corr_avgs), fpi_blocks)
        self.ainv, self.ainv_err = self.get_errors_from_blocks(0.092/self.fpi_fit(self.A_corr_avgs), np.divide(0.092,fpi_blocks))
        print(f"F_pi is {self.fpi}/a +- {self.fpi_err}/a")
        print(f"m_pi is {self.mpi_from_fpi}/a +- {self.mpi_from_fpi_err}/a")
        print(f"L is {self.fpi / 0.092 * 0.19732698044404103 * self.Nx} +- {self.fpi_err / 0.092 * 0.19732698044404103 * self.Nx} fm")
        print(f"a inverse is {self.ainv[0]} +- {self.ainv_err[0]} GeV")
        self.calc_pion_corrs()
        pion_mass_blocks = self.get_jackknife_blocks(self.pi_corrs, self.find_mass_from_fit)
        self.mpi_over_fpi, self.mpi_over_fpi_err = self.get_errors_from_blocks(self.find_mass_from_fit(self.pi_corr_avgs)[0]/self.fpi_fit(self.A_corr_avgs)[0], np.divide(pion_mass_blocks,fpi_blocks)[:,0])
        print(f"m_pi over f_pi is {self.mpi_over_fpi} +- {self.mpi_over_fpi_err}")
        self.fpi_done=True
    
    def fpi2_fit(self, pi_corr_avgs, fit_range=[], fpi_guess=0.2, mpi_guess=0.2):
        if(len(fit_range)==0):
            fit_range=[1,int(self.Nt)]
        a = np.power(pi_corr_avgs[fit_range[0]:fit_range[1]], 0.5)
        a = np.nan_to_num(a)
        nt = range(fit_range[0],fit_range[1])
        pi_opt, pi_cov = curve_fit(lambda nt, fpi, mpi : self.fpi2_model(nt,fpi,mpi), 
                                   nt, a, p0=[fpi_guess, mpi_guess])
        return pi_opt
    
    def calc_fpi2(self):
        if(self.fpi2_done):
            return
        self.calc_pion_corrs()
        fpi2_blocks = self.get_jackknife_blocks(self.pi_corrs, self.fpi2_fit)
        [[self.fpi2, self.mpi_from_fpi2],[self.fpi2_err,self.mpi_from_fpi2_err]] = self.get_errors_from_blocks(self.fpi2_fit(self.pi_corr_avgs), fpi2_blocks)
        print(f"F_pi (from pion correlators) is {self.fpi2}/a +- {self.fpi2_err}/a")
        print(f"m_pi is {self.mpi_from_fpi2}/a +- {self.mpi_from_fpi2_err}/a")
        self.fpi_done=True
    
    def plot_pion_fit(self):
        self.calc_pion_corrs()
        self.calc_pion_mass()
        pi_corr_avg_blocks = self.get_jackknife_blocks(self.pi_corrs)
        pi_corr_avgs_bc, pi_corr_errs = self.get_errors_from_blocks(self.pi_corr_avgs, pi_corr_avg_blocks)
        plt.errorbar(range(0,len(self.pi_corr_avgs)), pi_corr_avgs_bc, yerr=pi_corr_errs)
        nt_plt = np.arange(0,len(self.pi_corr_avgs)-1+0.1,0.1)
        plt.plot(nt_plt, self.cosh_model(nt_plt, self.A_pion, self.pion_mass))
        plt.title(f"Pion correlator fit \n {self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{self.date}-{self.day}_{self.version}")
        
    def plot_sigma_fit(self):
        self.calc_sigma_corrs()
        self.calc_sigma_mass()
        s_corr_avg_blocks = self.get_jackknife_blocks(self.s_corrs)
        s_corr_avgs_bc, s_corr_errs = self.get_errors_from_blocks(self.s_corr_avgs, s_corr_avg_blocks)
        plt.errorbar(range(0,len(self.s_corr_avgs)), s_corr_avgs_bc, yerr=s_corr_errs)
        nt_plt = np.arange(0,len(self.s_corr_avgs)-1+0.1,0.1)
        plt.plot(nt_plt, self.cosh_model(nt_plt, self.A_sigma, self.sigma_mass))
        plt.title(f"Sigma correlator fit \n {self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{self.date}-{self.day}_{self.version}")
    
    def plot_fpi_fit(self):
        self.calc_axpi_corrs()
        self.calc_fpi()
        A_blocks = self.get_jackknife_blocks(self.A_corrs, lambda d : d[1:int(self.Nt)])
        a,b = self.get_errors_from_blocks(self.A_corr_avgs[1:int(self.Nt)],A_blocks)
        plt.errorbar(range(1,self.Nt), a, yerr=b)
        plt.plot(range(1,self.Nt), self.fpi_model(range(1,int(self.Nt)), self.fpi, self.mpi_from_fpi))
        plt.title(f"F_pi fit \n {self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{self.date}-{self.day}_{self.version}")
    
    def scan_pion_mass(self, k_pion_mass=[]):
        if(self.scan_pion_mass_done):
            return
        self.pion_masses=[]
        self.pion_mass_errors=[]
        if(len(k_pion_mass)==0):
            self.k_pion_mass=range(1,100)
        else:
            self.k_pion_mass=k_pion_mass
        self.calc_pion_corrs()
        block_size = self.block_size
        for k in self.k_pion_mass:
            print(k)
            self.block_size = k
            blocks = self.get_jackknife_blocks(self.pi_corrs, self.find_mass_from_fit)
            [[mk,Ak],[ek,eAk]] = self.get_errors_from_blocks(self.find_mass_from_fit(self.pi_corr_avgs), blocks)
            self.pion_masses.append(mk)
            self.pion_mass_errors.append(ek)
        self.block_size = block_size
        self.scan_pion_mass_done = True
        
    def scan_phi0(self, k_phi0=[]):
        if(self.scan_phi0_done):
            return
        self.phi0s=[]
        self.phi0_errors=[]
        if(len(k_phi0)==0):
            self.k_phi0=range(1,100)
        else:
            self.k_phi0=k_phi0
        block_size = self.block_size
        for k in self.k_phi0:
            self.block_size = k
            blocks = self.get_jackknife_blocks([self.phi_list[i][0] for i in range(self.cutoff,self.data_len)])
            [pk, epk] = self.get_errors_from_blocks(np.mean([self.phi_list[i][0] for i in range(self.cutoff,self.data_len)]), blocks)
            self.phi0s.append(pk)
            self.phi0_errors.append(epk)
        self.block_size = block_size
        self.scan_phi0_done = True
    
    def plot_pion_scan(self):
        self.scan_pion_mass()
        plt.errorbar(self.k_pion_mass, self.pion_masses, yerr=self.pion_mass_errors)
        plt.title(f"Pion Mass Error vs Block Size \n {self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{self.date}-{self.day}_{self.version}")
    
    def plot_phi0_scan(self):
        self.scan_phi0()
        plt.errorbar(self.k_phi0, self.phi0s, yerr=self.phi0_errors)
        plt.title(f"Sigma VEV Error vs Block Size \n {self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{self.date}-{self.day}_{self.version}")

    def save(self):
        date = datetime.datetime.now().date()
        with open(f"analysis_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{date}.bin", "wb") as file:
            pickle.dump([self.loaded_files,
                self.sigma_vev_done,
                self.sigma_mass_eff_done,
                self.pion_corrs_done,
                self.sigma_corrs_done,
                self.pion_mass_done,
                self.sigma_mass_done,
                self.A_corrs_done,
                self.fpi_done,
                self.fpi2_done,
                self.scan_pion_mass_done,
                self.scan_phi0_done,
                self.pipi_sigma_corrs_done,
                self.pipi_pipi_corrs_done,
                self.phisq_done,
                self.pipim_sigma_corrs_done,
                self.pipim_sigma_energy_done,
                self.pipi_sigma_energy_done,
                self.psqm_done,
                self.vev_sigma,
                self.vev_sigma_err,
                self.sigma_mass_eff,
                self.sigma_mass_eff_err,
                self.psq,
                self.psq_err,
                self.phisq_av,
                self.phisq_av_err,
                self.psqm_av,
                self.psqm_av_err,
                self.pi_corrs,
                self.pi_corr_avgs,
                self.s_corrs,
                self.s_corr_avgs,
                self.pipi_s_corrs,
                self.pipi_s_corr_avgs,
                self.pipi_pipi_corrs,
                self.pipi_pipi_corr_avgs,
                self.pipim_s_corrs,
                self.pipim_s_corr_avgs,
                self.pion_mass_est,
                self.sigma_mass_est,
                self.fpi_est,
                self.pion_mass,
                self.A_pion,
                self.pion_mass_err,
                self.A_pion_err,
                self.sigma_mass,
                self.A_sigma,
                self.sigma_mass_err,
                self.A_sigma_err,
                self.pipim_sigma_E,
                self.pipim_sigma_E_err,
                self.pipi_sigma_E,
                self.pipi_sigma_E_err,
                self.A_corrs,
                self.A_corr_avgs,
                self.fpi,
                self.mpi_from_fpi,
                self.fpi_err,
                self.mpi_from_fpi_err,
                self.fpi2,
                self.mpi_from_fpi2,
                self.fpi2_err,
                self.mpi_from_fpi2_err,
                self.pion_masses,
                self.pion_mass_errors,
                self.k_pion_mass,
                self.phi0s,
                self.phi0_errors,
                self.k_phi0], file)
    
    def load(self, date):
        with open(f"analysis_{self.Nx}x{self.Nt}_msq_{self.msq}_lmbd_{self.lmbd}_alph_{self.alpha}_{date}.bin","rb") as file:
            data = pickle.load(file)
            for f in data[0]:
                self.load_data("-","-",f)
            [self.loaded_files,
            self.sigma_vev_done,
            self.sigma_mass_eff_done,
            self.pion_corrs_done,
            self.sigma_corrs_done,
            self.pion_mass_done,
            self.sigma_mass_done,
            self.A_corrs_done,
            self.fpi_done,
            self.fpi2_done,
            self.scan_pion_mass_done,
            self.scan_phi0_done,
            self.pipi_sigma_corrs_done,
            self.pipi_pipi_corrs_done,
            self.phisq_done,
            self.pipim_sigma_corrs_done,
            self.pipim_sigma_energy_done,
            self.psqm_done,
            self.vev_sigma,
            self.vev_sigma_err,
            self.sigma_mass_eff,
            self.sigma_mass_eff_err,
            self.psq,
            self.psq_err,
            self.phisq_av,
            self.phisq_av_err,
            self.psqm_av,
            self.psqm_av_err,
            self.pi_corrs,
            self.pi_corr_avgs,
            self.s_corrs,
            self.s_corr_avgs,
            self.pipi_s_corrs,
            self.pipi_s_corr_avgs,
            self.pipi_pipi_corrs,
            self.pipi_pipi_corr_avgs,
            self.pipim_s_corrs,
            self.pipim_s_corr_avgs,
            self.pion_mass_est,
            self.sigma_mass_est,
            self.fpi_est,
            self.pion_mass,
            self.A_pion,
            self.pion_mass_err,
            self.A_pion_err,
            self.sigma_mass,
            self.A_sigma,
            self.sigma_mass_err,
            self.A_sigma_err,
            self.pipim_sigma_E,
            self.pipim_sigma_E_err,
            self.A_corrs,
            self.A_corr_avgs,
            self.fpi,
            self.mpi_from_fpi,
            self.fpi_err,
            self.mpi_from_fpi_err,
            self.fpi2,
            self.mpi_from_fpi2,
            self.fpi2_err,
            self.mpi_from_fpi2_err,
            self.pion_masses,
            self.pion_mass_errors,
            self.k_pion_mass,
            self.phi0s,
            self.phi0_errors,
            self.k_phi0]=data
