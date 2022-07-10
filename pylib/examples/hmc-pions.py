#!/usr/bin/env python3

import sys
import numpy as np
import pickle
import datetime

import qlat as q

class Field_fft:
    def __init__(self, geo, mult=1):
        self.field = q.Field("double", geo, mult)
        self.field_ft = q.Field("Complex", geo, mult)
        self.updated = True
        self.updated_ft = True
        self.fft = q.mk_fft(True, is_normalizing=True)
        self.ifft = q.mk_fft(False, is_normalizing=True)
        self.V = geo.total_volume()
    
    def geo(self):
        return self.field.geo()
    
    def get_field(self):
        if(not self.updated):
            self.field.set_double_from_complex(self.ifft*self.field_ft)
            self.updated = True
        return self.field
    
    def get_field_ft(self):
        if(not self.updated_ft):
            self.field_ft.set_complex_from_double(self.field)
            self.field_ft = self.fft*self.field_ft
            self.updated_ft = True
        return self.field_ft
    
    def set_unit(self):
        q.set_unit(self.field)
        self.updated_ft = False
        self.updated = True
    
    def set_rand_momentum(self, action, masses, rng_state):
        action.hmc_set_rand_momentum(self.field_ft, masses, rng_state)
        # Performing the inverse Fourier transform this way projects
        # to real momenta
        self.field.set_double_from_complex(self.ifft*self.field_ft)
        self.field *= 2**0.5
        self.updated_ft = False
        self.updated = True
    
    def set_field(self, f):
        self.field @= f
        self.updated_ft = False
        self.updated = True
    
    def hmc_evolve(self, action, momentum, masses, dt):
        if(not self.updated_ft):
            self.field_ft.set_complex_from_double(self.field)
            self.field_ft = self.fft*self.field_ft
            self.updated_ft = True
        action.hmc_field_evolve(self.field_ft, momentum.get_field_ft(), masses, dt)
        self.updated = False
    
    def hmc_set_force(self, action, field):
        action.hmc_set_force(self.field, field)
        self.updated = True
        self.updated_ft = False
    
    def hmc_predict_field(self, action, momentum, masses, vev):
        if(not self.updated_ft):
            self.field_ft.set_complex_from_double(self.field)
            self.field_ft = self.fft*self.field_ft
        action.hmc_predict_field(self.field_ft, momentum.get_field_ft(), masses, vev)
        self.updated_ft = True
        self.updated = False
    
    def vacuum_subtract(self, vev):
        if(not self.updated_ft):
            self.field_ft.set_complex_from_double(self.field)
            self.field_ft = self.fft*self.field_ft
            self.updated_ft = True
        self.field_ft.set_elem([0,0,0,0],0,np.array([self.field_ft.get_elem([0,0,0,0],0)-vev*self.V**0.5], dtype='c16').tobytes())
        self.updated = False
    
    def load(self, path):
        self.field.load_double(path)
        self.updated_ft = False
        self.updated = True
    
    def add(self, f2):
        if(not self.updated):
            self.field.set_double_from_complex(self.ifft*self.field_ft)
            self.updated = True
        self.field += f2
        self.updated_ft = False
    
    def multiply(self, factor):
        if(not self.updated):
            self.field.set_double_from_complex(self.ifft*self.field_ft)
            self.updated = True
        self.field *= factor
        self.updated_ft = False
    
    def add_ft(self, f2):
        if(not self.updated_ft):
            self.field_ft.set_complex_from_double(self.field)
            self.field_ft = self.fft*self.field_ft
            self.updated_ft = True
        self.field_ft += f2
        self.updated = False
    
    def glb_sum(self):
        if(not self.updated):
            self.field.set_double_from_complex(self.ifft*self.field_ft)
            self.updated = True
        return self.field.glb_sum()

class HMC:
    def __init__(self, m_sq, lmbd, alpha, total_site, mult, steps):
        self.action = q.ScalarAction(m_sq, lmbd, alpha)
        self.m_sq = m_sq
        self.lmbd = lmbd
        self.alpha = alpha
        self.mult = mult
        self.steps = steps
        # Create the geometry for the field
        geo = q.Geometry(total_site, mult)
        # Create a random number generator that can be split between 
        # different portions of the lattice
        self.rs = q.RngState("test_hmc_pions-{}x{}x{}x{}".format(total_site[0], total_site[1], total_site[2], total_site[3]))
        
        # Save the spacial volume and the total volume of the lattice for 
        # future use
        self.Vx = total_site[0]*total_site[1]*total_site[2]
        self.V = self.Vx*total_site[3]
        self.total_site = total_site
        
        # Create the scalar field and set all field values to 1
        self.field = Field_fft(geo,mult)
        self.field.set_unit();
        #self.field.load(f"output_data/hmc-pions-sigma-pi-corrs_{total_site[0]}x{total_site[3]}_msq_{m_sq}_lmbd_{lmbd}_alph_{alpha}.field")
        
        self.f0 = Field_fft(self.field.geo(), self.mult)
        
        # Create a field to store field configurations predicted based on 
        # the initial momenta (with the assumption of harmonic evolution)
        self.field_predicted = Field_fft(geo,mult)
        self.field_predicted.set_unit()
        
        # Create a field to store the masses used for Fourier acceleration
        self.masses = q.Field("double",geo,mult)
        q.set_unit(self.masses);
        # Create axuillary fields to store mass estimates
        self.masses_est = q.Field("double",geo,mult)
        self.masses_avg = q.Field("double",geo,mult)
        # Create a field to keep track of which modes may have unreliable
        # mass estimates
        self.mask = q.Field("double",geo,mult)
        # Create a field to store the divisors when taking averages
        self.divisors = q.Field("double",geo,mult)
        self.divisors.set_unit()
        # The number of trajectories to calculate before taking measurements
        self.start_measurements = 0
        self.init_length = 10
        self.block_length = 45
        self.num_blocks = 2
        self.block_length2 = 50
        self.num_blocks2 = 2
        # A variable to store the estimated vacuum expectation value of sigma
        self.vev = 0

        # True if the system is approximately thermalized (as estimated from
        # the stability of the vacuum expectation value of sigma)
        self.thermalized = False
        self.traj = 1
        self.estimate_masses = True
        self.perform_metro = False
    
    @q.timer_verbose
    def run_traj(self):
        if(self.traj<self.init_length):
            self.run_hmc(self.masses, self.rs.split("hmc-{}".format(self.traj)))
            self.vevs = [self.vev]
            self.masses_avg @= self.masses
        elif(self.traj<self.init_length+self.num_blocks*self.block_length):
            self.perform_metro = True
            if((self.traj-self.init_length) % self.block_length == 0):
                self.divisors.invert_double()
                self.masses_avg.multiply_double(self.divisors)
                self.masses @= self.masses_avg
                q.set_zero(self.masses_avg)
                q.set_zero(self.divisors)
                self.vev=np.mean(self.vevs)
                self.vevs=[]
            if((self.traj-self.init_length) % self.block_length < self.init_length):
                self.estimate_masses = False
                self.run_hmc(self.masses, self.rs.split("hmc-{}".format(self.traj)))
            else:
                self.run_hmc_w_mass_est()
        elif(self.traj<self.init_length+self.num_blocks*self.block_length+self.num_blocks2*self.block_length2):
            if((self.traj - self.init_length - self.num_blocks*self.block_length) % self.block_length2 == 0):
                self.divisors.invert_double()
                self.masses_avg.multiply_double(self.divisors)
                self.masses @= self.masses_avg
                q.set_zero(self.masses_avg)
                q.set_zero(self.divisors)
                self.vev=np.mean(self.vevs)
                self.vevs=[]
            if((self.traj-self.init_length-self.num_blocks*self.block_length) % self.block_length2 < self.init_length):
                self.estimate_masses = False
                self.run_hmc(self.masses, self.rs.split("hmc-{}".format(self.traj)))
            else:
                self.run_hmc_w_mass_est2()
        else:
            self.estimate_masses = False
            if(self.traj==self.init_length+self.num_blocks*self.block_length+self.num_blocks2*self.block_length2):
                self.divisors.invert_double()
                self.masses_avg.multiply_double(self.divisors)
                self.masses @= self.masses_avg
                self.vev=np.mean(self.vevs)
                self.masses.save_double(f"output_data/masses_{self.total_site[0]}x{self.total_site[3]}_msq_{self.m_sq}_lmbd_{self.lmbd}_alph_{self.alpha}_{datetime.datetime.now().date()}.field")
            self.run_hmc(self.masses, self.rs.split("hmc-{}".format(self.traj)))
        self.traj += 1
    
    def run_hmc_w_mass_est(self):
        self.estimate_masses = True
        # Keep using the same estimated masses for every trajectory 
        # in this block
        self.masses_est @= self.masses
        self.run_hmc(self.masses_est, self.rs.split("hmc-est-mass{}".format(self.traj)))
        # Take the average of all the estimated masses in this block
        self.masses_est.multiply_double(self.mask)
        self.masses_avg+=self.masses_est
        self.divisors+=self.mask
        self.vevs.append(self.field.glb_sum()[0]/self.V)
        self.display_masses(self.masses_est)
   
    def run_hmc_w_mass_est2(self):
        # The mass estimation won't happen within self.run_hmc
        self.estimate_masses = False
        momentum0 = Field_fft(self.field.geo(), self.mult)
        momentum0.set_rand_momentum(self.action, self.masses, self.rs.split(f"set_rand_momentum{self.traj}"))
        
        self.masses_est.set_double_from_complex(momentum0.get_field_ft())
        
        self.run_hmc(self.masses, self.rs.split("hmc-est-mass{}".format(self.traj)), momentum0)
        
        vac_sub = q.Field("double", self.field.geo(), self.mult)
        vac_sub.set_double_from_complex(self.field.get_field_ft())
        q.displayln_info("=============================================================================+")
        q.displayln_info(vac_sub.get_elem([0,0,0,0],0))
        vac_sub.set_elem([0,0,0,0],0,np.array(vac_sub.get_elem([0,0,0,0],0)-self.vev*self.V**0.5, dtype="float").tobytes())
        q.displayln_info(vac_sub.get_elem([0,0,0,0],0))
        self.masses_est.multiply_double(vac_sub)
        self.masses_est*=2.0/np.pi
        self.masses_avg+=self.masses_est
        
        vac_sub.multiply_double(vac_sub)
        self.divisors+=vac_sub
        
        self.vevs.append(self.field.glb_sum()[0]/self.V)
        self.display_masses(self.masses_est)
    
    @q.timer_verbose
    def run_hmc(self, masses, rs, momentum0=False):
        # Create a copy of the scalar field
        self.f0.set_field(self.field.get_field())
        
        if(not momentum0):
            momentum = Field_fft(self.field.geo(), self.mult)
            
            # Create a random momentum field, distributed according to a 
            # special distribution that is Gaussian for any given point in
            # the lattice momentum space, but varries in width from point to 
            # point because of the momentum-dependent mass term used for Fourier
            # acceleration
            momentum.set_rand_momentum(self.action, masses, rs.split("set_rand_momentum"))
        else:
            momentum = momentum0
        
        momentum_ft = momentum.get_field_ft()
        momentums.append([[momentum_ft.get_elem([0,0,0,0],0),momentum_ft.get_elem([1,0,0,0],0),momentum_ft.get_elem([0,2,0,0],0),momentum_ft.get_elem([3,0,0,0],0)],
                          [momentum_ft.get_elem([0,0,0,0],1),momentum_ft.get_elem([1,0,0,0],1),momentum_ft.get_elem([0,2,0,0],1),momentum_ft.get_elem([3,0,0,0],1)]])
        
        # Predicts the field value at the end of the trajectory based on the
        # assumption that the evolution is a perfect harmonic oscillator
        self.field_predicted.hmc_predict_field(self.action, momentum, masses, self.vev)

        field_predicted_ft = self.field_predicted.get_field_ft()
        fields_pred.append([[field_predicted_ft.get_elem([0,0,0,0],0),field_predicted_ft.get_elem([1,0,0,0],0),field_predicted_ft.get_elem([0,2,0,0],0),field_predicted_ft.get_elem([3,0,0,0],0)],
                            [field_predicted_ft.get_elem([0,0,0,0],1),field_predicted_ft.get_elem([1,0,0,0],1),field_predicted_ft.get_elem([0,2,0,0],1),field_predicted_ft.get_elem([3,0,0,0],1)]])
        
        # Evolve the field over time md_time using the given momenta and 
        # the Hamiltonian appropriate for the given action
        delta_h = self.run_hmc_evolve(self.f0, momentum, masses, rs)
        
        # Decide whether to accept or reject the field update using the 
        # metropolis algorithm
        flag, accept_prob = self.metropolis_accept(delta_h, self.traj, rs.split("metropolis_accept"))

        accept_rates.append(accept_prob)
        field_ft = self.f0.get_field_ft()
        fields.append([[field_ft.get_elem([0,0,0,0],0),field_ft.get_elem([1,0,0,0],0),field_ft.get_elem([0,2,0,0],0),field_ft.get_elem([3,0,0,0],0)],
                       [field_ft.get_elem([0,0,0,0],1),field_ft.get_elem([1,0,0,0],1),field_ft.get_elem([0,2,0,0],1),field_ft.get_elem([3,0,0,0],1)]])
        
        # If the field update is accepted or we are within the first few 
        # trajectories, save the field update
        if flag or not self.perform_metro:
            q.displayln_info(f"run_hmc: update field (traj={self.traj})")
            self.field.set_field(self.f0.get_field())
    
    @q.timer_verbose
    def run_hmc_evolve(self, field, momentum, masses, rs):
        # Calculate the value of the molecular dynamics Hamiltonian for the 
        # initial field and momentum configuration
        energy = self.action.hmc_m_hamilton_node(momentum.get_field_ft(), masses) + self.action.action_node(field.get_field())
        
        # Evolve the field forward in molecular dynamics time using the 
        # given momenta and the Hamiltonian appropriate for the action
        dt = 1/self.steps
        masses_new, mask = self.hmc_evolve(field, momentum, masses)
        
        # Calculate the change in the value of the molecular dynamics 
        # Hamilton after the evolution 
        delta_h = self.action.hmc_m_hamilton_node(momentum.get_field_ft(), masses) + self.action.action_node(field.get_field()) - energy;
        
        if(self.estimate_masses):
            # Save the new estimated masses
            masses@=masses_new
            self.mask@=mask
        
        # Sum over delta_h for every parallel node (each node handles part 
        # of the lattice)
        delta_h = q.glb_sum(delta_h)
        
        return delta_h
    
    @q.timer_verbose
    def hmc_evolve(self, field, momentum, masses):
        # Evolve the field according to the given action using the force 
        # gradient algorithm
        lam = 0.5 * (1.0 - 1.0 / 3.0**0.5);
        theta = (2.0 - 3.0**0.5) / 48.0;
        dt = 1/self.steps
        ttheta = theta * dt * dt * dt;
        # The Fourier transformed field is updated, and then the field is 
        # updated based on the new Fourier transformed field
        field.hmc_evolve(self.action, momentum, masses, lam*dt)
        # Save a list of the vacuum expectation value of sigma at each point
        # on the trajectory
        if(self.estimate_masses):
            vevs=[]
            geo = field.geo()
            # Fields to store the average of the modulus of the field
            # and the forces over the course of the trajectory
            field_mod_av = q.Field("double",geo,geo.multiplicity())
            force_mod_av = q.Field("double",geo,geo.multiplicity())
            field_mod_av.set_zero()
            force_mod_av.set_zero()
            # Fields to store the modulus of the field and the forces over 
            # the course of the trajectory
            field_mod = q.Field("double",geo,geo.multiplicity())
            force_mod = q.Field("double",geo,geo.multiplicity())
            # A field to keep track of which modes may have unreliable
            # mass estimates
            mask = q.Field("double",geo,geo.multiplicity())
        for i in range(self.steps):
            self.sm_evolve(field, momentum, 4.0 * ttheta / dt, 0.5 * dt)
            field.hmc_evolve(self.action, momentum, masses, (1.0 - 2.0 * lam) * dt)
            force = self.sm_evolve(field, momentum, 4.0 * ttheta / dt, 0.5 * dt)
            if i < steps - 1:
                field.hmc_evolve(self.action, momentum, masses, 2.0 * lam * dt)
            else:
                field.hmc_evolve(self.action, momentum, masses, lam * dt)
            if(self.estimate_masses):
                # Add the modulus of the current field configuration to the 
                # average (for the field, we subtract the vaccum expectation
                # value of sigma from the zero mode)
                field.vacuum_subtract(self.vev)
                field_mod.set_abs_from_complex(field.get_field_ft())
                # Add back the vacuum expectation value of sigma to the field
                field.vacuum_subtract(-self.vev)
                field_mod*=1/steps
                field_mod_av+=field_mod
                force_mod.set_abs_from_complex(force.get_field_ft())
                force_mod*=1/steps
                force_mod_av+=force_mod
                # Save data to estimate the v.e.v. of sigma
                field_sum = field.glb_sum()
                vevs.append(field_sum[0]/field.V)
        if(self.estimate_masses):
            # Estimate the masses we should use in order to evolve each field 
            # mode by half of its period
            masses_new = q.Field("double",field.geo(),self.mult)
            masses_new.set_ratio_double(force_mod_av, field_mod_av)
            # If any of the values of field_mod_av are small, the estimated
            # mass may be unreliable. Using the inverse masses to set the 
            # scale, we keep track of where the possibly unreliable masses
            # are in the mask field (a value of 1 means reliable, and 0 means
            # unreliable)
            masses.invert_double()
            masses*=0.01
            masses.less_than_double(field_mod_av, mask)
            masses*=100.0
            masses.invert_double()
            # After multiplying the ratio of force_mod_av/field_mod_av by
            # (pi/2)**(-2), we have our estimated masses
            masses_new *= 4/np.pi**2
            self.vev = np.mean(vevs)
            return masses_new, mask
        else:
            return 0, 0
    
    @q.timer_verbose
    def sm_evolve(self, field_init, momentum, fg_dt, dt):
        # Evolve the momentum field according to the given action using the  
        # force gradient algorithm
        field = q.Field("double",field_init.geo(),self.mult)
        field = field_init.get_field().copy()
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
    
    def display_masses(self, masses):
        q.displayln_info("Estmiated masses:")
        q.displayln_info([masses.get_elem([0,0,0,0],0),masses.get_elem([1,0,0,0],0),masses.get_elem([5,0,0,0],0)])

def phi_squared(field,action):
    # Calculate the average value of phi^2
    phi_sq = action.sum_sq(field) # Returns sum of field^2/2
    phi_sq = q.glb_sum(phi_sq) # Sums over all nodes
    geo = field.geo()
    return phi_sq/geo.total_volume()/geo.multiplicity()

@q.timer_verbose
def main():
    hmc = HMC(m_sq,lmbd,alpha,total_site,mult,steps)
    
    # Create the geometry for the axial current field
    geo_cur = q.Geometry(total_site, 3)
    # This field will store the calculated axial currents
    axial_current = q.Field("double",geo_cur)
    # This field will store the calculated sigma and pion fields
    sigma_pions = q.Field("double",hmc.field.geo())
    
    for traj in range(1,n_traj+1):
        # Run the HMC algorithm to update the field configuration
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
        
        #
        hmc.action.get_sigma_pions(sigma_pions, hmc.field.get_field())
        tslices_sp = sigma_pions.glb_sum_tslice()
        
        q.displayln_info("Average sigma/pion:")
        field_sum2 = sigma_pions.glb_sum()
        phi2=[field_sum2[i]/hmc.field.V for i in range(mult)]
        q.displayln_info([sigma_pions.get_elem([4,0,0,0],0),sigma_pions.get_elem([4,0,0,0],1),sigma_pions.get_elem([4,0,0,0],2),sigma_pions.get_elem([4,0,0,0],3)])
        #q.displayln_info(phi2)
        
        hmc.action.get_sigma_pions(sigma_pions, hmc.field_predicted.get_field())
        tslices_sp_predicted = sigma_pions.glb_sum_tslice()
        
        #q.displayln_info("Analytic masses (free case):")
        #ms=[calc_mass([0,0,0,0]),calc_mass([1,0,0,0]),calc_mass([2,0,0,0]),calc_mass([3,0,0,0]),calc_mass([4,0,0,0]),calc_mass([5,0,0,0])]
        #q.displayln_info(ms)
        
        q.displayln_info("Estmiated masses:")
        ms=[hmc.masses.get_elem([0,0,0,0],0),hmc.masses.get_elem([1,0,0,0],0),hmc.masses.get_elem([2,0,0,0],0),hmc.masses.get_elem([3,0,0,0],0),hmc.masses.get_elem([4,0,0,0],0)]
        q.displayln_info(ms)
        ms=[hmc.masses.get_elem([0,0,0,0],1),hmc.masses.get_elem([1,0,0,0],1),hmc.masses.get_elem([2,0,0,0],1),hmc.masses.get_elem([3,0,0,0],1),hmc.masses.get_elem([4,0,0,0],1)]
        q.displayln_info(ms)
        
        tslices = hmc.field.get_field().glb_sum_tslice()
        tslices_predicted = hmc.field_predicted.get_field().glb_sum_tslice()
        
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
            sp_list.append(tslices_sp.to_numpy())
            sp_pred_list.append(tslices_sp_predicted.to_numpy())
    
    # Saves the final field configuration so that the next run can be 
    # started where this one left off
    hmc.field.get_field().save_double(f"output_data/hmc-pions-sigma-pi-corrs_{total_site[0]}x{total_site[3]}_msq_{m_sq}_lmbd_{lmbd}_alph_{alpha}.field")

# Stores the average phi^2 for each trajectory
psq_list=[]
psq_pred_list=[]
# Stores the average values of each phi_i for each trajectory
phi_list=[]
phi_pred_list=[]
# Stores the timeslice sums of each phi_i for each trajectory
timeslices=[]
timeslices_pred=[]
# Stores the timeslice sums of the time component of each of the three
# axial currents for each trajectory
ax_cur_timeslices=[]
ax_cur_timeslices_pred=[]
# Sigma pion timeslices
sp_list = []
sp_pred_list = []
# Save the acceptance rates
accept_rates=[]
fields=[]
fields_pred=[]
momentums=[]

# The lattice dimensions
total_site = [8,8,8,16]

# The multiplicity of the scalar field
mult = 4

# The number of trajectories to calculate
n_traj = 500
# The number of steps to take in a single trajectory
steps = 20

# Use action for a Euclidean scalar field. The Lagrangian will be:
# (1/2)*[sum i]|dphi_i|^2 + (1/2)*m_sq*[sum i]|phi_i|^2
#     + (1/24)*lmbd*([sum i]|phi_i|^2)^2
m_sq = -1.0
lmbd = 1.0
alpha = 0.1

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
    except:
        raise Exception("Invalid arguments: use -d for lattice dimensions, -n for multiplicity, -t for number of trajectories, -m for mass squared, -l for lambda, -a for alpha, and -s for the number of steps in a trajectory. e.g. python hmc-pions.py -l 8x8x8x16 -n 4 -t 50 -m -1.0 -l 1.0 -a 0.1")

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

with open(f"output_data/sigma_pion_corrs_{total_site[0]}x{total_site[3]}_msq_{m_sq}_lmbd_{lmbd}_alph_{alpha}_{datetime.datetime.now().date()}.bin", "wb") as output:
    pickle.dump([accept_rates,psq_list,phi_list,sp_list,timeslices,ax_cur_timeslices,psq_pred_list,phi_pred_list,sp_pred_list,timeslices_pred,ax_cur_timeslices_pred,fields,momentums,fields_pred],output)

q.timer_display()

q.end()
