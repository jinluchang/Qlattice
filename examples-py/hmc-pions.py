#!/usr/bin/env python3

import sys
import numpy as np
import pickle
import datetime
import glob

import qlat as q


class Field_fft:
    def __init__(self, geo, mult=1):
        self.field = q.Field("Double", geo, mult)
        self.field_ft = q.Field("Complex", geo, mult)
        self.updated = True
        self.updated_ft = True
        self.fft = q.mk_fft(True, is_normalizing=True)
        self.ifft = q.mk_fft(False, is_normalizing=True)
        self.V = geo.total_volume()
        self.mult = mult
    
    def geo(self):
        return self.field.geo()
    
    def get_field(self):
        if(not self.updated):
            q.field_double.set_double_from_complex(self.field, self.ifft*self.field_ft)
            self.updated = True
        return self.field
    
    def get_field_ft(self):
        if(not self.updated_ft):
            q.field_double.set_complex_from_double(self.field_ft,self.field)
            self.field_ft = self.fft*self.field_ft
            self.updated_ft = True
        return self.field_ft
    
    def set_unit(self):
        q.set_unit(self.field)
        self.updated_ft = False
        self.updated = True
    
    def set_zero(self):
        q.set_zero(self.field)
        self.updated_ft = False
        self.updated = True
    
    def set_rand_momentum(self, action, masses, rng_state):
        action.hmc_set_rand_momentum(self.field_ft, masses, rng_state)
        # Performing the inverse Fourier transform this way projects
        # to real momenta
        q.field_double.set_double_from_complex(self.field,self.ifft*self.field_ft)
        self.field *= 2**0.5
        self.updated_ft = False
        self.updated = True
    
    def set_field(self, f):
        self.field @= f
        self.updated_ft = False
        self.updated = True
    
    def set_field_ft(self, f):
        self.field_ft @= f
        self.updated = False
        self.updated_ft = True
    
    def hmc_evolve(self, action, momentum, masses, dt):
        self.get_field_ft()
        action.hmc_field_evolve(self.field_ft, momentum.get_field_ft(), masses, dt)
        self.updated = False
    
    def hmc_set_force(self, action, field):
        action.hmc_set_force(self.field, field)
        self.updated = True
        self.updated_ft = False
    
    def hmc_predict_field(self, action, momentum, masses, vev):
        action.hmc_predict_field(self.field_ft, momentum.get_field_ft(), masses, vev)
        self.updated_ft = True
        self.updated = False
    
    def vacuum_subtract(self, vev):
        self.get_field_ft()
        self.field_ft.set_elem([0,0,0,0],0,np.array([self.field_ft.get_elem([0,0,0,0],0)-vev*self.V**0.5], dtype='c16').tobytes())
        self.updated = False
    
    def load(self, path):
        self.field.load_double(path)
        self.updated_ft = False
        self.updated = True
    
    def add(self, f2):
        self.get_field()
        self.field += f2
        self.updated_ft = False
    
    def multiply(self, factor):
        self.get_field()
        self.field *= factor
        self.updated_ft = False
    
    def add_ft(self, f2):
        self.get_field_ft()
        self.field_ft += f2
        self.updated = False
    
    def glb_sum(self):
        self.get_field()
        return self.field.glb_sum()
    
    def remove_low_modes(self):
        total_site = self.geo().total_site()
        self.get_field_ft()
        for m in range(self.mult):
            for x in [[0,0,0,0],[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[total_site[0]-1,0,0,0],[0,total_site[1]-1,0,0],[0,0,total_site[2]-1,0],[0,0,0,total_site[3]-1]]:
                self.field_ft.set_elem(x,m,np.array([0.0j], dtype='c16').tobytes())
        self.updated = False
    
    def get_representatives_ft(self):
        field_ft = self.get_field_ft()
        return [[field_ft.get_elem([0,0,0,0],0),field_ft.get_elem([1,0,0,0],0),field_ft.get_elem([0,2,0,0],0),field_ft.get_elem([3,0,0,0],0)],
                [field_ft.get_elem([0,0,0,0],1),field_ft.get_elem([1,0,0,0],1),field_ft.get_elem([0,2,0,0],1),field_ft.get_elem([3,0,0,0],1)]]

class HMC:
    def __init__(self, m_sq, lmbd, alpha, total_site, mult, steps, mass_force_coef, recalculate_masses, fresh_start):
        self.m_sq = m_sq
        self.lmbd = lmbd
        self.alpha = alpha
        self.mult = mult
        self.steps = steps
        self.mass_force_coef = mass_force_coef
        self.Vx = total_site[0]*total_site[1]*total_site[2]
        self.V = self.Vx*total_site[3]
        self.total_site = total_site
        
        self.fileid = f"{self.total_site[0]}x{self.total_site[3]}_msq_{self.m_sq}_lmbd_{self.lmbd}_alph_{self.alpha}_{date}_{version}"
        self.fileidwc = f"{self.total_site[0]}x{self.total_site[3]}_msq_{self.m_sq}_lmbd_{self.lmbd}_alph_{self.alpha}_*_{version}"
        
        # The number of trajectories to calculate before taking measurements
        self.start_measurements = 0
        self.init_length = 20
        self.block_init_length = 20
        self.block_length = 95
        self.num_blocks = 4
        self.final_block_length = 200
        # A variable to store the estimated vacuum expectation value of sigma
        self.vev = 0
        self.vevs=[self.vev]
        self.vev_est = 0

        self.traj = 1
        self.estimate_masses = True
        self.safe_estimate_masses = True
        self.masses_loaded = False
        self.perform_metro = False
        
        self.action = q.ScalarAction(m_sq, lmbd, alpha)
        geo = q.Geometry(total_site, mult)
        # Create a random number generator that can be split between 
        # different portions of the lattice
        self.rs = q.RngState("test_hmc_pions-{}x{}x{}x{}".format(total_site[0], total_site[1], total_site[2], total_site[3]))
        
        # Create the scalar field and set all field values to 1
        self.field = Field_fft(geo,mult)
        
        # Create a field to store field configurations predicted based on 
        # the initial momenta (with the assumption of harmonic evolution)
        self.field_predicted = Field_fft(geo,mult)
        self.field_predicted.set_unit()
        
        # Create a field to store the masses used for Fourier acceleration
        self.masses = q.Field("Double",geo,mult)
        # Create a field to store the estimated optimal Fourier accleration
        # masses before lower bounds are applied
        self.masses_est = q.Field("Double",geo,mult)
        
        if(fresh_start):
            self.masses.set_unit()
            self.field.set_unit()
        elif(recalculate_masses):
            self.masses.set_unit()
            self.load_field()
        else:
            self.load_field()
            self.load_masses()
        
        # Create an auxiliary field to store the field as it evolves
        self.f0 = Field_fft(self.field.geo(), self.mult)
        
        # Fields to store everything we need to do linear regression to
        # determine what HMC masses to use for Fourier acceleration
        self.field_av = q.Field("Double",geo,geo.multiplicity())
        self.force_av = q.Field("Double",geo,geo.multiplicity())
        self.field_sq_av = q.Field("Double",geo,geo.multiplicity())
        self.force_mod_av = q.Field("Double",geo,geo.multiplicity())
        self.field_mod_av = q.Field("Double",geo,geo.multiplicity())
        self.field_force_cor = q.Field("Double",geo,geo.multiplicity())
        self.divisor = 0
        self.mask = q.Field("Double",geo,geo.multiplicity())
        self.aux1 = q.Field("Double",geo,geo.multiplicity())
        self.aux2 = q.Field("Double",geo,geo.multiplicity())
        self.reset_fit_variables()
    
    @q.timer_verbose
    def run_traj(self):
        if(not self.masses_loaded and self.traj<self.init_length):
            self.run_hmc(self.rs.split("1hmc-{}".format(self.traj)))
            self.update_masses_w_safe_fit()
        elif(not self.masses_loaded and self.traj<self.init_length+self.num_blocks*self.block_length):
            self.perform_metro = True
            self.safe_estimate_masses = True
            if((self.traj-self.init_length) % self.block_length == 0 and (self.traj-self.init_length) > 0):
                self.vev=np.mean(self.vevs)
                self.vevs=[]
                self.update_masses_w_fit()
            if((self.traj-self.init_length) % self.block_length < self.block_init_length):
                self.estimate_masses = False
                self.run_hmc(self.rs.split("2hmc-{}".format(self.traj)))
            else:
                self.run_hmc_w_mass_est()
        elif(not self.masses_loaded and self.traj<self.init_length+self.num_blocks*self.block_length+self.final_block_length):
            if(self.traj - self.init_length - self.num_blocks*self.block_length == 0):
                self.update_masses_w_fit()
                self.vev=np.mean(self.vevs)
                self.vevs=[]
            if(self.traj-self.init_length-self.num_blocks*self.block_length < self.block_init_length):
                self.estimate_masses = False
                self.run_hmc(self.rs.split("3hmc-{}".format(self.traj)))
            else:
                self.run_hmc_w_mass_est()
        else:
            self.perform_metro = True
            self.estimate_masses = False
            self.safe_estimate_masses = False
            if(not self.masses_loaded and self.traj==self.init_length+self.num_blocks*self.block_length+self.final_block_length):
                self.update_masses_w_fit()
                self.vev=np.mean(self.vevs)
                self.masses.save_double(f"output_data/masses_{self.fileid}.field")
                self.masses_est.save_double(f"output_data/masses_wo_lower_bound_{self.fileid}.field")
            self.run_hmc(self.rs.split("4hmc-{}".format(self.traj)))
        self.traj += 1
    
    def update_masses_w_fit(self):
        # Estimate the masses we should use in order to evolve each field 
        # mode by half of its period
        self.force_av.multiply_double(self.field_av)
        self.force_av*=1/self.divisor
        self.field_force_cor-=self.force_av
        self.field_av.multiply_double(self.field_av)
        self.field_av*=1/self.divisor
        self.field_sq_av-=self.field_av
        q.field_double.set_ratio_double(self.masses,self.field_force_cor,self.field_sq_av)
        # After multiplying the ratio of force_mod_av/field_mod_av by
        # (pi/2)**(-2), we have our estimated masses
        self.masses *= 4/np.pi**2
        self.masses_est@=self.masses
        # A safer method for estimating the masses away from equilibrium
        self.force_mod_av*=1.0/self.divisor/self.mass_force_coef
        # Choose the larger of the two choices
        self.choose_larger(self.masses, self.force_mod_av)
        self.aux2.set_unit()
        self.aux2*=0.01
        self.choose_larger(self.masses, self.aux2)
        
        self.reset_fit_variables()
    
    def choose_larger(self, field1, field2):
        # TODO: Make a more efficient, safer way to do it
        q.field_double.less_than_double(field1, field2, self.mask)
        field2.multiply_double(self.mask)
        self.aux1.set_unit()
        q.field_double.less_than_double(self.mask, self.aux1, self.mask)
        field1.multiply_double(self.mask)
        field1+=field2
    
    def update_masses_w_safe_fit(self):
        self.update_masses_w_fit()
    
    def reset_fit_variables(self):
        # Reset all accumulated variables
        self.field_av.set_zero()
        self.force_av.set_zero()
        self.field_sq_av.set_zero()
        self.force_mod_av.set_zero()
        self.field_mod_av.set_zero()
        self.field_force_cor.set_zero()
        self.divisor = 0
    
    def load_masses(self):
        filename = self.find_latest_date(f"output_data/masses_{self.fileidwc}.field")
        if(not filename==""):
            self.masses.load_double(filename)
            self.masses_loaded=True
        else:
            self.masses.set_unit()
    
    def load_field(self):
        filename, self.traj = self.find_latest_traj(f"output_data/fields/hmc_pions_traj_*_{self.fileidwc}.field")
        self.init_length+=self.traj-1
        if(not filename==""):
            self.field.load(filename)
        else:
            self.field.set_unit()
    
    def save_field(self):
        self.field.get_field().save_double(f"output_data/fields/hmc_pions_traj_{self.traj}_{self.fileid}.field")
    
    def find_latest_date(self, filewc):
        files = glob.glob(filewc)
        if not len(files):
            return ""
        datenos=[]
        for f in files:
            info=f.split("_")
            for i in info:
                date = i.split("-")
                if(len(date)==3 and len(date[0])==4 and len(date[1])==2 and len(date[2])==2):
                    try:
                        datenos.append(float(date[0])+float(date[1])/12.0+float(date[2])/31.0)
                        break
                    except ValueError:
                        pass
        return files[datenos.index(max(datenos))]
    
    def find_latest_traj(self, filewc):
        files = glob.glob(filewc)
        if not len(files):
            return "", 1
        trajnos=[]
        for f in files:
            info=f.split("_")
            for i in range(len(info)):
                if(info[i]=="traj"):
                    try:
                        trajnos.append(int(info[i+1]))
                        break
                    except ValueError:
                        pass
                    except IndexError:
                        pass
        traj = max(trajnos)
        return files[trajnos.index(traj)], traj
    
    def run_hmc_w_mass_est(self):
        self.estimate_masses = True
        # Keep using the same estimated masses for every trajectory 
        # in this block
        self.run_hmc(self.rs.split("hmc-est-mass{}".format(self.traj)))
        self.vevs.append(self.field.glb_sum()[0]/self.V)
    
    @q.timer_verbose
    def run_hmc(self, rs):
        # Create a copy of the scalar field
        self.f0.set_field(self.field.get_field())
        
        momentum = Field_fft(self.field.geo(), self.mult)
        
        # Create a random momentum field, distributed according to a 
        # special distribution that is Gaussian for any given point in
        # the lattice momentum space, but varries in width from point to 
        # point because of the momentum-dependent mass term used for Fourier
        # acceleration
        momentum.set_rand_momentum(self.action, self.masses, rs.split("set_rand_momentum"))
        momentums.append(momentum.get_representatives_ft())
        
        # Predicts the field value at the end of the trajectory based on the
        # assumption that the evolution is a perfect harmonic oscillator
        self.field_predicted.hmc_predict_field(self.action, momentum, self.masses, self.vev)
        
        fields_pred.append(self.field_predicted.get_representatives_ft())
        
        # Evolve the field over time md_time using the given momenta and 
        # the Hamiltonian appropriate for the given action
        delta_h = self.run_hmc_evolve(self.f0, momentum, rs)
        
        # Decide whether to accept or reject the field update using the 
        # metropolis algorithm
        flag, accept_prob = self.metropolis_accept(delta_h, self.traj, rs.split("metropolis_accept"))

        accept_rates.append(accept_prob)
        
        fields.append(self.f0.get_representatives_ft())
        
        # If the field update is accepted or we are within the first few 
        # trajectories, save the field update
        if flag or not self.perform_metro:
            q.displayln_info(f"run_hmc: update field (traj={self.traj})")
            self.field.set_field(self.f0.get_field())
    
    @q.timer_verbose
    def run_hmc_evolve(self, field, momentum, rs):
        # Calculate the value of the molecular dynamics Hamiltonian for the 
        # initial field and momentum configuration
        energy = self.action.hmc_m_hamilton_node(momentum.get_field_ft(), self.masses) + self.action.action_node(field.get_field())
        
        # Evolve the field forward in molecular dynamics time using the 
        # given momenta and the Hamiltonian appropriate for the action
        dt = 1/self.steps
        self.hmc_evolve(field, momentum)
        
        # Calculate the change in the value of the molecular dynamics 
        # Hamilton after the evolution 
        delta_h = self.action.hmc_m_hamilton_node(momentum.get_field_ft(), self.masses) + self.action.action_node(field.get_field()) - energy;
        
        # Sum over delta_h for every parallel node (each node handles part 
        # of the lattice)
        delta_h = q.glb_sum(delta_h)
        
        return delta_h
    
    @q.timer_verbose
    def hmc_evolve(self, field, momentum):
        # Evolve the field according to the given action using the force 
        # gradient algorithm
        lam = 0.5 * (1.0 - 1.0 / 3.0**0.5);
        theta = (2.0 - 3.0**0.5) / 48.0;
        dt = 1/self.steps
        ttheta = theta * dt * dt * dt;
        # The Fourier transformed field is updated, and then the field is 
        # updated based on the new Fourier transformed field
        field.hmc_evolve(self.action, momentum, self.masses, lam*dt)
        for i in range(self.steps):
            self.sm_evolve(field, momentum, 4.0 * ttheta / dt, 0.5 * dt)
            field.hmc_evolve(self.action, momentum, self.masses, (1.0 - 2.0 * lam) * dt)
            force = self.sm_evolve(field, momentum, 4.0 * ttheta / dt, 0.5 * dt)
            if i < steps - 1:
                field.hmc_evolve(self.action, momentum, self.masses, 2.0 * lam * dt)
            else:
                field.hmc_evolve(self.action, momentum, self.masses, lam * dt)
            if(self.safe_estimate_masses):
                q.field_double.set_abs_from_complex(self.aux1,force.get_field_ft())
                self.force_mod_av+=self.aux1
                q.field_double.set_abs_from_complex(self.aux1,field.get_field_ft())
                self.field_mod_av+=self.aux1
                self.divisor+=1
            if(self.estimate_masses):
                q.field_double.set_double_from_complex(self.aux1,field.get_field_ft())
                self.field_av+=self.aux1
                q.field_double.set_double_from_complex(self.aux2,force.get_field_ft())
                self.force_av+=self.aux2
                self.aux2.multiply_double(self.aux1)
                self.field_force_cor+=self.aux2
                self.aux1.multiply_double(self.aux1)
                self.field_sq_av+=self.aux1
                self.divisor+=1
    
    @q.timer_verbose
    def sm_evolve(self, field_init, momentum, fg_dt, dt):
        # Evolve the momentum field according to the given action using the  
        # force gradient algorithm
        field = q.Field("Double",field_init.geo(),self.mult)
        field @= field_init.get_field()
        force = Field_fft(field_init.geo(), self.mult)
        force.hmc_set_force(self.action, field)
        
        force.multiply(fg_dt)
        field += force.get_field()
        force.multiply(1/fg_dt)
        
        force.hmc_set_force(self.action, field)
        
        force.multiply(-dt)
        momentum.add(force.get_field())
        force.multiply(-1/dt)
        return force
    
    @q.timer_verbose
    def metropolis_accept(self, delta_h, traj, rs):
        # This variable will store whether or not we've decided to accept 
        # the trajectory
        flag_d = 0.0
        # This variable will store the acceptance probability
        accept_prob = 0.0
        # Since we only want to decide once whether or not to accept, we 
        # only run this code on node 0 (delta_h has already been summed over
        # all nodes)
        if q.get_id_node() == 0:
            if delta_h <= 0.0:
                accept_prob = 1.0
                flag_d = 1.0
            else:
                accept_prob = np.exp(-delta_h)
                rand_num = rs.u_rand_gen(1.0, 0.0)
                if rand_num <= accept_prob:
                    flag_d = 1.0
        flag_d = q.glb_sum(flag_d)
        accept_prob = q.glb_sum(accept_prob)
        # If we decided to accept, flag_d should be 1.0. Therefore, flag is
        # true if we decided to accept.
        flag = flag_d > 0.5
        q.displayln_info("metropolis_accept: flag={:d} with accept_prob={:.1f}% delta_h={:.16f} traj={:d}".format(
            flag, accept_prob * 100.0, delta_h, traj))
        return flag, accept_prob
    
    def display_masses(self, msg, masses):
        q.displayln_info(msg)
        q.displayln_info([masses.get_elem([0,0,0,0],0),masses.get_elem([1,0,0,0],0),masses.get_elem([4,0,0,0],0)])
        q.displayln_info([masses.get_elem([0,0,0,0],1),masses.get_elem([1,0,0,0],1),masses.get_elem([4,0,0,0],1)])

def phi_squared(field,action):
    # Calculate the average value of phi^2
    phi_sq = action.sum_sq(field) # Returns sum of field^2/2
    phi_sq = q.glb_sum(phi_sq) # Sums over all nodes
    geo = field.geo()
    return phi_sq/geo.total_volume()/geo.multiplicity()

def histogram_bin(val,midpoint,n):
    j = 1
    for i in range(n):
        if(val<j*midpoint/2**i):
            j = j*2 - 1
        else:
            j = j*2 + 1
    return int((j-1)/2.0)

def update_phi_sq_dist(elems,vev_sigma,norm_factor):
    phi_sq = 0.0
    for elem in elems:
        phi_sq+=elem**2
    phi_sq_dist[histogram_bin(phi_sq,np.abs(vev_sigma)*50,6)]+=1.0/norm_factor

def update_phi_i_dist(phi,vev_sigma,norm_factor):
    phi_i_dist[histogram_bin(np.abs(phi),np.abs(vev_sigma)*50,6)]+=1.0/norm_factor

def update_theta_dist(elems,norm_factor):
    phi_sq = 0.0
    for elem in elems:
        phi_sq+=elem**2
    for elem in elems[1:]:
        theta_dist[histogram_bin(np.pi/2+np.arccos(elems[0]/phi_sq**0.5),np.pi/2,6)]+=1.0/norm_factor

def save_observables():
    with open(f"output_data/sigma_pion_corrs_{total_site[0]}x{total_site[3]}_msq_{m_sq}_lmbd_{lmbd}_alph_{alpha}_{date}_{version}.bin", "wb") as output:
        pickle.dump({"trajs": trajs,
                    "accept_rates": accept_rates, 
                    "psq_list": psq_list, 
                    "phi_list": phi_list, 
                    "timeslices": timeslices, 
                    "hm_timeslices": hm_timeslices,
                    "ax_cur_timeslices": ax_cur_timeslices,
                    "polar_timeslices": polar_timeslices,
                    "phi_sq_dist": phi_sq_dist, 
                    "phi_i_dist": phi_i_dist,
                    "theta_dist": theta_dist, 
                    "psq_pred_list": psq_pred_list,
                    "phi_pred_list": phi_pred_list, 
                    "timeslices_pred": timeslices_pred,
                    "hm_timeslices_pred": hm_timeslices_pred,
                    "ax_cur_timeslices_pred": ax_cur_timeslices_pred,
                    "fields": fields,
                    "momentums": momentums,
                    "field_pred": fields_pred},output)

def load_observables():
    filename = f"output_data/sigma_pion_corrs_{total_site[0]}x{total_site[3]}_msq_{m_sq}_lmbd_{lmbd}_alph_{alpha}_{date}_{version}.bin"
    if len(glob.glob(filename)):
        with open(filename,"rb") as input:
            data = pickle.load(input)
            trajs.extend(data["trajs"])
            accept_rates.extend(data["accept_rates"])
            psq_list.extend(data["psq_list"])
            phi_list.extend(data["phi_list"])
            timeslices.extend(data["timeslices"])
            hm_timeslices.extend(data["hm_timeslices"])
            ax_cur_timeslices.extend(data["ax_cur_timeslices"])
            polar_timeslices.extend(data["polar_timeslices"])
            phi_sq_dist.extend(data["phi_sq_dist"])
            phi_i_dist.extend(data["phi_i_dist"])
            theta_dist.extend(data["theta_dist"])
            psq_pred_list.extend(data["psq_pred_list"])
            phi_pred_list.extend(data["phi_pred_list"])
            timeslices_pred.extend(data["timeslices_pred"])
            hm_timeslices_pred.extend(data["hm_timeslices_pred"])
            ax_cur_timeslices_pred.extend(data["ax_cur_timeslices_pred"])
            fields.extend(data["fields"])
            momentums.extend(data["momentums"])
            fields_pred.extend(data["field_pred"])

@q.timer_verbose
def main():
    # If observables have been saved from a previous calculation (on the
    # same day), then load that file first
    load_observables()
    print(len(psq_list))
    
    hmc = HMC(m_sq,lmbd,alpha,total_site,mult,steps,mass_force_coef,recalculate_masses,fresh_start)
    
    # Create the geometry for the axial current field
    geo_cur = q.Geometry(total_site, 3)
    # This field will store the calculated axial currents
    axial_current = q.Field("Double",geo_cur)
    #
    hm_field = Field_fft(hmc.field.geo(),4)
    hm_field_pred = Field_fft(hmc.field.geo(),4)
    #
    polar_field = q.Field("Double",hmc.field.geo())
    
    for traj in range(1,n_traj+1):
        # Run the HMC algorithm to update the field configuration
        trajs.append(hmc.traj)
        hmc.run_traj()

        # Calculate the expectation values of phi and phi^2
        q.displayln_info("Average phi^2:")
        psq = phi_squared(hmc.field.get_field(), hmc.action)
        q.displayln_info(psq)
        q.displayln_info("Average phi^2 predicted:")
        psq_predicted = phi_squared(hmc.field_predicted.get_field(), hmc.action)
        q.displayln_info(psq_predicted)

        q.displayln_info("Average phi:")
        field_sum = hmc.field.get_field().glb_sum()
        phi=[field_sum[i]/hmc.field.V for i in range(mult)]
        q.displayln_info([hmc.field.get_field().get_elem([4,0,0,0],0),hmc.field.get_field().get_elem([4,0,0,0],1),hmc.field.get_field().get_elem([4,0,0,0],2),hmc.field.get_field().get_elem([4,0,0,0],3)])
        #q.displayln_info(phi)
        field_sum = hmc.field_predicted.get_field().glb_sum()
        phi_predicted=[field_sum[i]/hmc.field.V for i in range(mult)]
        
        #q.displayln_info("Analytic masses (free case):")
        #ms=[calc_mass([0,0,0,0]),calc_mass([1,0,0,0]),calc_mass([2,0,0,0]),calc_mass([3,0,0,0]),calc_mass([4,0,0,0]),calc_mass([5,0,0,0])]
        #q.displayln_info(ms)
        
        hmc.display_masses("Masses:", hmc.masses)
        q.displayln_info("vev: ")
        q.displayln_info(hmc.vev)
        
        tslices = hmc.field.get_field().glb_sum_tslice()
        #
        hm_field.set_field_ft(hmc.field.get_field_ft())
        hm_field.remove_low_modes()
        hm_tslices = hm_field.get_field().glb_sum_tslice()
        #
        tslices_predicted = hmc.field_predicted.get_field().glb_sum_tslice()
        #
        hm_field_pred.set_field_ft(hmc.field_predicted.get_field_ft())
        hm_field_pred.remove_low_modes()
        hm_tslices_pred = hm_field_pred.get_field().glb_sum_tslice()
        #
        hmc.action.get_polar_field(polar_field, hmc.field.get_field())
        polar_tslices = polar_field.glb_sum_tslice()
                
        # Calculate the axial current of the current field configuration
        # and save it in axial_current
        hmc.action.axial_current_node(axial_current, hmc.field.get_field())
        tslices_ax_cur = axial_current.glb_sum_tslice()
        hmc.action.axial_current_node(axial_current, hmc.field_predicted.get_field())
        tslices_ax_cur_predicted = axial_current.glb_sum_tslice()
        
        if traj>hmc.start_measurements:
            psq_list.append(psq)
            phi_list.append(phi)
            timeslices.append(tslices.to_numpy())
            ax_cur_timeslices.append(tslices_ax_cur.to_numpy())
            psq_pred_list.append(psq_predicted)
            phi_pred_list.append(phi_predicted)
            timeslices_pred.append(tslices_predicted.to_numpy())
            ax_cur_timeslices_pred.append(tslices_ax_cur_predicted.to_numpy())
            hm_timeslices.append(hm_tslices.to_numpy())
            hm_timeslices_pred.append(hm_tslices_pred.to_numpy())
            polar_timeslices.append(polar_tslices.to_numpy())
        if traj>hmc.init_length+hmc.num_blocks*hmc.block_length+hmc.final_block_length:
            field = hmc.field.get_field()
            norm_factor = hmc.V*(n_traj+1-hmc.init_length-hmc.num_blocks*hmc.block_length-hmc.final_block_length)
            for x in range(hmc.total_site[0]):
                for y in range(hmc.total_site[1]):
                    for z in range(hmc.total_site[2]):
                        for t in range(hmc.total_site[3]):
                            elems = field.get_elems([x,y,z,t])
                            update_phi_sq_dist(elems,hmc.vev,norm_factor)
                            update_phi_i_dist(elems[1],hmc.vev,norm_factor)
                            update_phi_i_dist(elems[2],hmc.vev,norm_factor)
                            update_phi_i_dist(elems[3],hmc.vev,norm_factor)
                            update_theta_dist(elems,norm_factor)
        if traj%50 == 0:
            hmc.save_field()
            save_observables()
    
    # Saves the final field configuration so that the next run can be 
    # started where this one left off
    hmc.save_field()
    save_observables()

# Stores the trajectory number for debugging purposes
trajs = []
# Stores the average phi^2 for each trajectory
psq_list=[]
psq_pred_list=[]
# Stores the average values of each phi_i for each trajectory
phi_list=[]
phi_pred_list=[]
# Stores the timeslice sums of each phi_i for each trajectory
timeslices=[]
timeslices_pred=[]
# Stores timeslices sums calculated using only high modes
hm_timeslices=[]
hm_timeslices_pred=[]
# Stores the timeslice sums of the time component of each of the three
# axial currents for each trajectory
ax_cur_timeslices=[]
ax_cur_timeslices_pred=[]
# Stores the timeslice sums of the polar-coordinate fields
polar_timeslices=[]
# Save the acceptance rates
accept_rates=[]
fields=[]
forces=[]
fields_pred=[]
momentums=[]
phi_sq_dist=[0.0]*64
phi_i_dist=[0.0]*64
theta_dist=[0.0]*64

# The lattice dimensions
total_site = [4,4,4,8]

# The multiplicity of the scalar field
mult = 4

# The number of trajectories to calculate
n_traj = 1000
# The number of steps to take in a single trajectory
steps = 20
# The factor by which to scale down the force when setting a lower limit
# on Fourier acceleration masses
mass_force_coef = 1000.0

# Use action for a Euclidean scalar field. The Lagrangian will be:
# (1/2)*[sum i]|dphi_i|^2 + (1/2)*m_sq*[sum i]|phi_i|^2
#     + (1/24)*lmbd*([sum i]|phi_i|^2)^2
m_sq = -8.
lmbd = 32.0
alpha = 0.1

recalculate_masses = False
fresh_start = False

version = "1-6"
date = datetime.datetime.now().date()

for i in range(1,len(sys.argv),2):
    try:
        if(sys.argv[i]=="-d"):
            a = sys.argv[i+1].split("x")
            total_site = [int(a[j]) for j in range(4)]
        elif(sys.argv[i]=="-n"):
            mult = int(sys.argv[i+1])
        elif(sys.argv[i]=="-t"):
            n_traj = int(sys.argv[i+1])
        elif(sys.argv[i]=="-m"):
            m_sq = float(sys.argv[i+1])
        elif(sys.argv[i]=="-l"):
            lmbd = float(sys.argv[i+1])
        elif(sys.argv[i]=="-a"):
            alpha = float(sys.argv[i+1])
        elif(sys.argv[i]=="-s"):
            steps = int(sys.argv[i+1])
        elif(sys.argv[i]=="-f"):
            mass_force_coef = float(sys.argv[i+1])
        elif(sys.argv[i]=="-r"):
            recalculate_masses = True
        elif(sys.argv[i]=="-R"):
            fresh_start = True
    except:
        raise Exception("Invalid arguments: use -d for lattice dimensions, -n for multiplicity, -t for number of trajectories, -m for mass squared, -l for lambda, -a for alpha, -s for the number of steps in a trajectory, -f for the factor by which to scale down the force when setting a lower limit for Fourier acceleration masses, -r to force recalculating the masses, and -R to force recalculating the masses and the initial field. e.g. python hmc-pions.py -l 8x8x8x16 -n 4 -t 50 -m -1.0 -l 1.0 -a 0.1 -f 100.0")

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 2, 2, 2],
        [2, 2, 2, 2],
        [2, 2, 2, 4]]

q.begin(sys.argv, size_node_list)

q.show_machine()

q.qremove_all_info("results")

main()

q.timer_display()

q.displayln_info(f"CHECK: Simulation finished successfully.")

q.end()
