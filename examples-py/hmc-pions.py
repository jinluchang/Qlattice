#!/usr/bin/env python3

import sys
import numpy as np
import pickle
import datetime
import glob
import fnmatch

import qlat as q

class Field_fft:
    def __init__(self, geo, mult=1):
        self.field = q.Field(q.ElemTypeRealD, geo, mult)
        self.field_ft = q.Field(q.ElemTypeComplexD, geo, mult)
        self.updated = True
        self.updated_ft = True
        self.fft = q.mk_fft(True, is_normalizing=True)
        self.ifft = q.mk_fft(False, is_normalizing=True)
        self.V = geo.total_volume
        self.mult = mult

    @property
    def geo(self):
        return self.field.geo

    @property
    def multiplicity(self):
        return self.field.multiplicity

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
        self.field_ft.set_elem_xg(q.Coordinate([0,0,0,0]),0,np.array([self.field_ft.get_elem_xg([0,0,0,0],0)-vev*self.V**0.5], dtype='c16'))
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
        total_site = self.geo.total_site
        self.get_field_ft()
        for m in range(self.mult):
            for x in [[0,0,0,0],[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[total_site[0]-1,0,0,0],[0,total_site[1]-1,0,0],[0,0,total_site[2]-1,0],[0,0,0,total_site[3]-1]]:
                self.field_ft.set_elem_xg(q.Coordinate(x),m,np.array([0.0j], dtype='c16'))
        self.updated = False

    def get_representatives_ft(self):
        field_ft = self.get_field_ft()
        return [[field_ft.get_elem_xg(q.Coordinate([0,0,0,0]),0)[:].item(),
                 field_ft.get_elem_xg(q.Coordinate([1,0,0,0]),0)[:].item(),
                 field_ft.get_elem_xg(q.Coordinate([0,2,0,0]),0)[:].item(),
                 field_ft.get_elem_xg(q.Coordinate([3,0,0,0]),0)[:].item()],
                [field_ft.get_elem_xg(q.Coordinate([0,0,0,0]),1)[:].item(),
                 field_ft.get_elem_xg(q.Coordinate([1,0,0,0]),1)[:].item(),
                 field_ft.get_elem_xg(q.Coordinate([0,2,0,0]),1)[:].item(),
                 field_ft.get_elem_xg(q.Coordinate([3,0,0,0]),1)[:].item()]]

class HMC:
    def __init__(self, m_sq, lmbd, alpha, total_site, mult, steps, mass_force_coef, recalculate_masses, fresh_start, block_sizes, version, date):
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
        self.fileidwc_no_version = f"{self.total_site[0]}x{self.total_site[3]}_msq_{self.m_sq}_lmbd_{self.lmbd}_alph_{self.alpha}_*"

        # The number of trajectories to calculate before taking measurements
        self.init_length = block_sizes[0]
        self.block_init_length = block_sizes[1]
        self.block_length = block_sizes[2]
        self.num_blocks = block_sizes[3]
        self.final_block_length = block_sizes[4]
        # A variable to store the estimated vacuum expectation value of sigma
        self.vev = 0
        self.vevs=[self.vev]
        self.vev_est = 0
        self.accept_prob = 0.0

        self.traj = 1
        self.estimate_masses = True
        self.safe_estimate_masses = True
        self.masses_loaded = False
        self.perform_metro = False

        self.action = q.ScalarAction(m_sq, lmbd, alpha)
        geo = q.Geometry(total_site)
        # Create a random number generator that can be split between
        # different portions of the lattice
        self.rs = q.RngState("test_hmc_pions-{}x{}x{}x{}".format(total_site[0], total_site[1], total_site[2], total_site[3]))

        # Create the scalar field and set all field values to 1
        self.field = Field_fft(geo,mult)
        self.init_momentum = Field_fft(geo,mult)

        # Create a field to store field configurations predicted based on
        # the initial momenta (with the assumption of harmonic evolution)
        self.field_predicted = Field_fft(geo,mult)
        self.field_predicted.set_unit()

        # Create a field to store the masses used for Fourier acceleration
        self.masses = q.Field(q.ElemTypeRealD,geo,mult)
        # Create a field to store the estimated optimal Fourier accleration
        # masses before lower bounds are applied
        self.masses_est = q.Field(q.ElemTypeRealD,geo,mult)

        if(fresh_start):
            self.masses.set_unit()
            self.field.set_unit()
        elif(recalculate_masses):
            self.masses.set_unit()
            self.load_field()
        else:
            self.load_masses()
            self.load_field()

        # Create an auxiliary field to store the field as it evolves
        self.f0 = Field_fft(self.field.geo, self.mult)

        # Fields to store everything we need to do linear regression to
        # determine what HMC masses to use for Fourier acceleration
        self.field_av = q.Field(q.ElemTypeRealD,geo,mult)
        self.force_av = q.Field(q.ElemTypeRealD,geo,mult)
        self.field_sq_av = q.Field(q.ElemTypeRealD,geo,mult)
        self.force_mod_av = q.Field(q.ElemTypeRealD,geo,mult)
        self.field_mod_av = q.Field(q.ElemTypeRealD,geo,mult)
        self.field_force_cor = q.Field(q.ElemTypeRealD,geo,mult)
        self.divisor = 0
        self.mask = q.Field(q.ElemTypeRealD,geo,mult)
        self.aux1 = q.Field(q.ElemTypeRealD,geo,mult)
        self.aux2 = q.Field(q.ElemTypeRealD,geo,mult)
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
        self.force_av *= self.field_av
        self.force_av*=1/self.divisor
        self.field_force_cor-=self.force_av
        self.field_av *= self.field_av
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
        field2 *= self.mask
        self.aux1.set_unit()
        q.field_double.less_than_double(self.mask, self.aux1, self.mask)
        field1 *= self.mask
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
        if(filename==""):
            filename, self.traj = self.find_latest_traj(f"output_data/fields/hmc_pions_traj_*_{self.fileidwc_no_version}.field")
        self.init_length+=self.traj-1
        if(not filename==""):
            self.field.load(filename)
            # If a field is loaded, avoid running any trajectories 
            # without a metropolis accept/reject step
            self.init_length = 0
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
                        datenos.append(float(date[0])+float(date[1])/12.0+float(date[2])/365.25)
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
        # Only update the trajectory number if the loaded field was produced
        # using the same version
        traj = max(trajnos)
        if(fnmatch.fnmatch(files[trajnos.index(traj)], f"output_data/fields/hmc_pions_traj_*_{self.fileidwc}.field")):
            traj_new = traj
        else:
            traj_new = 1
        return files[trajnos.index(traj)], traj_new

    def run_hmc_w_mass_est(self):
        self.estimate_masses = True
        # Keep using the same estimated masses for every trajectory
        # in this block
        self.run_hmc(self.rs.split("hmc-est-mass{}".format(self.traj)))
        self.vevs.append(self.field.glb_sum()[0, 0].item()/self.V)

    @q.timer_verbose
    def run_hmc(self, rs):
        # Create a copy of the scalar field
        self.f0.set_field(self.field.get_field())

        momentum = Field_fft(self.field.geo, self.mult)

        # Create a random momentum field, distributed according to a
        # special distribution that is Gaussian for any given point in
        # the lattice momentum space, but varries in width from point to
        # point because of the momentum-dependent mass term used for Fourier
        # acceleration
        momentum.set_rand_momentum(self.action, self.masses, rs.split("set_rand_momentum"))
        self.init_momentum.set_field(momentum.get_field())

        # Predicts the field value at the end of the trajectory based on the
        # assumption that the evolution is a perfect harmonic oscillator
        self.field_predicted.hmc_predict_field(self.action, momentum, self.masses, self.vev)
        
        # Evolve the field over time md_time using the given momenta and
        # the Hamiltonian appropriate for the given action
        delta_h = self.run_hmc_evolve(self.f0, momentum, rs)

        # Decide whether to accept or reject the field update using the
        # metropolis algorithm
        flag, self.accept_prob = self.metropolis_accept(delta_h, self.traj, rs.split("metropolis_accept"))

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
            if i < self.steps - 1:
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
                self.aux2 *= self.aux1
                self.field_force_cor+=self.aux2
                self.aux1 *= self.aux1
                self.field_sq_av+=self.aux1
                self.divisor+=1

    @q.timer_verbose
    def sm_evolve(self, field_init, momentum, fg_dt, dt):
        # Evolve the momentum field according to the given action using the
        # force gradient algorithm
        field = q.Field(q.ElemTypeRealD,field_init.geo,self.mult)
        field @= field_init.get_field()
        force = Field_fft(field_init.geo, self.mult)
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
        q.displayln_info([masses.get_elem_xg(q.Coordinate([0,0,0,0]),0)[:].item(),
                          masses.get_elem_xg(q.Coordinate([1,0,0,0]),0)[:].item(),
                          masses.get_elem_xg(q.Coordinate([4,0,0,0]),0)[:].item()])
        q.displayln_info([masses.get_elem_xg(q.Coordinate([0,0,0,0]),1)[:].item(),
                          masses.get_elem_xg(q.Coordinate([1,0,0,0]),1)[:].item(),
                          masses.get_elem_xg(q.Coordinate([4,0,0,0]),1)[:].item()])

class Measurements:
    def __init__(self, total_site, field_geo, multiplicity, save_file):
        self.save_file = save_file
        # Stores the trajectory number for debugging purposes
        self.traj_list = []
        # Save the acceptance rates
        self.accept_rates=[]
        # Stores the average phi^2 for each trajectory
        self.psq_list=[]
        self.psq_pred_list=[]
        # Stores the average values of each phi_i for each trajectory
        self.phi_list=[]
        self.phi_pred_list=[]
        # Stores the timeslice sums of each phi_i for each trajectory
        self.timeslices=[]
        self.timeslices_pred=[]
        # Stores timeslices sums calculated using only high modes
        self.hm_timeslices=[]
        self.hm_timeslices_pred=[]
        # Stores the timeslice sums of the time component of each of the three
        # axial currents for each trajectory
        self.ax_cur_timeslices=[]
        self.ax_cur_timeslices_pred=[]
        # Stores the timeslice sums of the polar-coordinate fields
        self.polar_timeslices=[]
        #
        self.kinematic_ms = [[1,0,0,0],
                            [0,1,0,0],
                            [0,0,1,0],
                            [1,1,0,0],
                            [1,0,1,0],
                            [0,1,1,0],
                            [1,1,1,0],
                            [-1,1,1,0],
                            [1,-1,1,0],
                            [1,1,-1,0],
                            [2,0,0,0],
                            [0,2,0,0],
                            [0,0,2,0],
                            [3,0,0,0],
                            [0,3,0,0],
                            [0,0,3,0]]
        self.timeslices_m = []
        #
        self.fields=[]
        self.forces=[]
        self.fields_pred=[]
        self.momentums=[]
        #
        self.psq_dist_center = 0
        self.phi_dist_center = 0
        self.phi_sq_dist=[0.0]*64
        self.phi_i_dist=[0.0]*64
        self.theta_dist=[0.0]*64
        #
        # Create fields to store only the high modes
        self.hm_field = Field_fft(field_geo,4)
        # Create a field to store the "spherical" version of the field
        self.polar_field = q.Field(q.ElemTypeRealD,field_geo, multiplicity)
        # Auxillary fields for use in calculations
        self.auxc = q.Field(q.ElemTypeComplexD,field_geo, multiplicity)
        self.auxd = q.Field(q.ElemTypeRealD,field_geo, multiplicity)
        # Create the geometry for the axial current field
        geo_cur = q.Geometry(total_site)
        # This field will store the calculated axial currents
        self.axial_current = q.Field(q.ElemTypeRealD, geo_cur, 3)
        # Create fields to project out momentum states
        self.mom_factors = []
        geo_m = q.Geometry(total_site)
        for m in self.kinematic_ms:
            self.mom_factors.append(q.mk_phase_field(geo_m, m))
    
    def measure(self, hmc):
        self.traj_list.append(hmc.traj)
        self.accept_rates.append(hmc.accept_prob)
        
        # Calculate the expectation values of phi^2
        self.psq_list.append(self.calc_psq(hmc.field.get_field(), hmc.action))
        self.psq_pred_list.append(self.calc_psq(hmc.field_predicted.get_field(), hmc.action))
        # Calculate the expectation value of phi
        field_sum = hmc.field.get_field().glb_sum()[0]
        self.phi_list.append([field_sum[i]/hmc.field.V for i in range(hmc.mult)])
        field_sum = hmc.field_predicted.get_field().glb_sum()[0]
        self.phi_pred_list.append([field_sum[i]/hmc.field.V for i in range(hmc.mult)])
        #
        self.momentums.append(hmc.init_momentum.get_representatives_ft())
        self.fields_pred.append(hmc.field_predicted.get_representatives_ft())
        self.fields.append(hmc.f0.get_representatives_ft())
        #
        self.timeslices.append(hmc.field.get_field().glb_sum_tslice().to_numpy())
        self.timeslices_pred.append(hmc.field_predicted.get_field().glb_sum_tslice().to_numpy())
        #
        self.hm_field.set_field_ft(hmc.field.get_field_ft())
        self.hm_field.remove_low_modes()
        self.hm_timeslices.append(self.hm_field.get_field().glb_sum_tslice().to_numpy())
        #
        self.hm_field.set_field_ft(hmc.field_predicted.get_field_ft())
        self.hm_field.remove_low_modes()
        self.hm_timeslices_pred.append(self.hm_field.get_field().glb_sum_tslice().to_numpy())
        #
        hmc.action.get_polar_field(self.polar_field, hmc.field.get_field())
        self.polar_timeslices.append(self.polar_field.glb_sum_tslice().to_numpy())
        #
        tslices_m = []
        for m in self.mom_factors:
            q.field_double.set_complex_from_double(self.auxc, hmc.field.get_field())
            self.auxc*=m
            tslices_m.append(self.auxc.glb_sum_tslice().to_numpy())
        self.timeslices_m.append(tslices_m)
        
        # Calculate the axial current of the current field configuration
        # and save it in axial_current
        hmc.action.axial_current_node(self.axial_current, hmc.field.get_field())
        self.ax_cur_timeslices.append(self.axial_current.glb_sum_tslice().to_numpy())
        hmc.action.axial_current_node(self.axial_current, hmc.field_predicted.get_field())
        self.ax_cur_timeslices_pred.append(self.axial_current.glb_sum_tslice().to_numpy())
        
        # Update lists that give histograms of phi, phi^2, and theta distributions
        if hmc.traj>hmc.init_length+hmc.num_blocks*hmc.block_length+hmc.final_block_length:
            field = hmc.field.get_field()
            if not self.psq_dist_center:
                self.psq_dist_center = self.psq_list[-1]
                q.field_double.set_complex_from_double(self.auxc, field)
                q.field_double.set_abs_from_complex(self.auxd, self.auxc)
                field_sum_abs = self.auxd.glb_sum()[0]
                phi_abs=[field_sum_abs[i]/hmc.field.V for i in range(hmc.mult)]
                self.phi_dist_center = (phi_abs[1]+phi_abs[2]+phi_abs[3])/3.0
            for x in range(hmc.total_site[0]):
                for y in range(hmc.total_site[1]):
                    for z in range(hmc.total_site[2]):
                        for t in range(hmc.total_site[3]):
                            elems = field.get_elems_xg(q.Coordinate([x,y,z,t]))[0]
                            self.update_phi_sq_dist(elems,hmc.V)
                            self.update_phi_i_dist(np.abs(elems[1]),hmc.V)
                            self.update_phi_i_dist(np.abs(elems[2]),hmc.V)
                            self.update_phi_i_dist(np.abs(elems[3]),hmc.V)
                            self.update_theta_dist(elems,hmc.V)
    
    def display_measurements(self):
        q.displayln_info("Average phi:")
        q.displayln_info(self.phi_list[-1])
        q.displayln_info("Average phi^2:")
        q.displayln_info(self.psq_list[-1])
        q.displayln_info("Average phi^2 predicted:")
        q.displayln_info(self.psq_pred_list[-1])
    
    def calc_psq(self,field,action):
        # Calculate the average value of phi^2
        phi_sq = action.sum_sq(field) # Returns sum of field^2/2
        phi_sq = q.glb_sum(phi_sq) # Sums over all nodes
        geo = field.geo
        return phi_sq/geo.total_volume
        
    def histogram_bin(self,val,midpoint,n):
        j = 1
        for i in range(n):
            if(val<j*midpoint/2**i):
                j = j*2 - 1
            else:
                j = j*2 + 1
        return int((j-1)/2.0)

    def update_phi_sq_dist(self,elems,norm_factor):
        phi_sq = 0.0
        for elem in elems:
            phi_sq+=elem**2
        self.phi_sq_dist[self.histogram_bin(phi_sq,self.psq_dist_center,6)]+=1.0/norm_factor

    def update_phi_i_dist(self,phi,norm_factor):
        self.phi_i_dist[self.histogram_bin(np.abs(phi),self.phi_dist_center,6)]+=1.0/norm_factor

    def update_theta_dist(self,elems,norm_factor):
        phi_sq = 0.0
        for elem in elems:
            phi_sq+=elem**2
        for elem in elems[1:]:
            self.theta_dist[self.histogram_bin(np.arccos(elems[0]/phi_sq**0.5),np.pi/2,6)]+=1.0/norm_factor
    
    def save(self):
        with open(self.save_file, "wb") as output:
            pickle.dump({"traj_list": self.traj_list,
                        "accept_rates": self.accept_rates,
                        "psq_list": self.psq_list,
                        "phi_list": self.phi_list,
                        "timeslices": self.timeslices,
                        "timeslices_m": self.timeslices_m,
                        "kinematic_ms": self.kinematic_ms,
                        "hm_timeslices": self.hm_timeslices,
                        "ax_cur_timeslices": self.ax_cur_timeslices,
                        "polar_timeslices": self.polar_timeslices,
                        "psq_dist_center": self.psq_dist_center,
                        "phi_dist_center": self.phi_dist_center,
                        "phi_sq_dist": self.phi_sq_dist,
                        "phi_i_dist": self.phi_i_dist,
                        "theta_dist": self.theta_dist,
                        "psq_pred_list": self.psq_pred_list,
                        "phi_pred_list": self.phi_pred_list,
                        "timeslices_pred": self.timeslices_pred,
                        "hm_timeslices_pred": self.hm_timeslices_pred,
                        "ax_cur_timeslices_pred": self.ax_cur_timeslices_pred,
                        "fields": self.fields,
                        "momentums": self.momentums,
                        "field_pred": self.fields_pred},output)
    
    def load(self):
        if len(glob.glob(self.save_file)):
            with open(self.save_file,"rb") as input:
                data = pickle.load(input)
                if(data["kinematic_ms"]!=self.kinematic_ms):
                    raise Exception("Combining data with different kinematic factors")
                self.traj_list.extend(data["traj_list"])
                self.accept_rates.extend(data["accept_rates"])
                self.psq_list.extend(data["psq_list"])
                self.phi_list.extend(data["phi_list"])
                self.timeslices.extend(data["timeslices"])
                self.timeslices_m.extend(data["timeslices_m"])
                self.hm_timeslices.extend(data["hm_timeslices"])
                self.ax_cur_timeslices.extend(data["ax_cur_timeslices"])
                self.polar_timeslices.extend(data["polar_timeslices"])
                self.psq_pred_list.extend(data["psq_pred_list"])
                self.phi_pred_list.extend(data["phi_pred_list"])
                self.timeslices_pred.extend(data["timeslices_pred"])
                self.hm_timeslices_pred.extend(data["hm_timeslices_pred"])
                self.ax_cur_timeslices_pred.extend(data["ax_cur_timeslices_pred"])
                self.fields.extend(data["fields"])
                self.momentums.extend(data["momentums"])
                self.fields_pred.extend(data["field_pred"])
                #
                self.psq_dist_center = data["psq_dist_center"]
                self.phi_dist_center = data["phi_dist_center"]
                self.phi_sq_dist = [self.phi_sq_dist[i] + data["phi_sq_dist"][i] for i in range(len(self.phi_sq_dist))]
                self.phi_i_dist = [self.phi_i_dist[i] + data["phi_i_dist"][i] for i in range(len(self.phi_i_dist))]
                self.theta_dist = [self.theta_dist[i] + data["theta_dist"][i] for i in range(len(self.theta_dist))]

@q.timer_verbose
def main():
    # The lattice dimensions
    total_site = q.Coordinate([4,4,4,8])
    # The multiplicity of the scalar field
    mult = 4
    # Use action for a Euclidean scalar field. The Lagrangian will be:
    # (1/2)*[sum i]|dphi_i|^2 + (1/2)*m_sq*[sum i]|phi_i|^2
    #     + (1/24)*lmbd*([sum i]|phi_i|^2)^2
    m_sq = -8.
    lmbd = 32.0
    alpha = 0.1
    # The number of trajectories to calculate
    n_traj = 1000
    #
    version = "3-1"
    date = datetime.datetime.now().date()
    # The number of steps to take in a single trajectory
    steps = 20
    # The factor by which to scale down the force when setting a lower limit
    # on Fourier acceleration masses
    mass_force_coef = 1000.0
    #
    init_length = 20
    block_init_length = 20
    block_length = 95
    num_blocks = 4
    final_block_length = 200
    recalculate_masses = False
    fresh_start = False

    for i in range(1,len(sys.argv)):
        try:
            if(sys.argv[i]=="-d"):
                a = sys.argv[i+1].split("x")
                total_site = q.Coordinate([int(a[j]) for j in range(4)])
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
            elif(sys.argv[i]=="-i"):
                init_length = int(sys.argv[i+1])
            elif(sys.argv[i]=="-I"):
                block_init_length = int(sys.argv[i+1])
            elif(sys.argv[i]=="-b"):
                block_length = int(sys.argv[i+1])
            elif(sys.argv[i]=="-N"):
                num_blocks = int(sys.argv[i+1])
            elif(sys.argv[i]=="-B"):
                final_block_length = int(sys.argv[i+1])
            elif(sys.argv[i]=="-S"):
                save_frequency = int(sys.argv[i+1])
        except:
            raise Exception("Invalid arguments: use -d for lattice dimensions, -n for multiplicity, -t for number of trajectories, -m for mass squared, -l for lambda, -a for alpha, -s for the number of steps in a trajectory, -f for the factor by which to scale down the force when setting a lower limit for Fourier acceleration masses, -r to force recalculating the masses, -R to force recalculating the masses and the initial field, -i for the number of trajectories to do at the beginning without a Metropolis step, -I for the number of trajectories to omit from the start of each HMC mass estimation block, -b for the number of trajectories in one HMC mass estimation block, -N for the number of HMC mass estimation blocks (excluding the final block), -B for the number of trajectories in the final mass estimation block, and -S for the number of trajectories between each save. e.g. python hmc-pions.py -l 8x8x8x16 -n 4 -t 50 -m -1.0 -l 1.0 -a 0.1 -f 100.0")
    
    hmc = HMC(m_sq,lmbd,alpha,total_site,mult,steps,mass_force_coef,recalculate_masses,fresh_start,[init_length,block_init_length,block_length,num_blocks,final_block_length],version,date)
    measurements = Measurements(total_site, hmc.field.geo, hmc.field.multiplicity, f"output_data/measurements_{hmc.fileid}.bin")
    
    # If observables have been saved from a previous calculation (on the
    # same day), then load that file first
    measurements.load()
    
    save_frequency = 50
    while hmc.traj <= n_traj:
        # Run the HMC algorithm to update the field configuration
        hmc.run_traj()
        hmc.display_masses("Masses:", hmc.masses)
        q.displayln_info("vev: ")
        q.displayln_info(hmc.vev)
        measurements.measure(hmc)
        
        if hmc.traj%save_frequency == 0:
            hmc.save_field()
            measurements.save()

    # Saves the final field configuration so that the next run can be
    # started where this one left off
    hmc.save_field()
    measurements.save()
    
    q.displayln_info(f"CHECK: The vacuum expectation value of phi_0 is {round(np.mean(measurements.phi_list[int(n_traj/2):], axis=0)[0],2)}.")
    q.displayln_info(f"CHECK: The vacuum expectation value of phi^2 is {round(np.mean(measurements.psq_list[int(n_traj/2):]),2)}.")

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 2, 2, 2],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        [2, 2, 2, 8],
        [2, 4, 4, 4]]

q.begin_with_mpi(size_node_list)

q.show_machine()

q.qremove_all_info("results")

main()

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: Simulation finished successfully.")
