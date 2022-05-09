#!/usr/bin/env python3

import sys
import math as m
import numpy as np
import pickle
import datetime

import qlat as q

def phi_squared(field,action):
    # Calculate the average value of phi^2
    phi_sq = action.sum_sq(field) # Returns sum of field^2/2
    phi_sq = q.glb_sum(phi_sq) # Sums over all nodes
    geo = field.geo()
    return phi_sq/geo.total_volume()/geo.multiplicity()

def calc_mass(pg):
	return 4/np.pi**2*(m_sq + 8 - 2*(np.cos(2*np.pi*pg[0]/total_site[0]) +
									np.cos(2*np.pi*pg[1]/total_site[1]) + 
									np.cos(2*np.pi*pg[2]/total_site[2]) + 
									np.cos(2*np.pi*pg[3]/total_site[3])))

@q.timer_verbose
def sm_evolve(momentum_ft, field_init, action, fg_dt, dt, fft, ifft):
    # Evolve the momentum field according to the given action using the  
    # force gradient algorithm
    geo = field_init.geo()
    V = geo.total_volume()
    field = q.Field("double",geo)
    field @= field_init
    momentum = q.Field("double",geo)
    momentum.set_double_from_complex(ifft*momentum_ft)
    momentum*=1/V**0.5
    #
    force = q.Field("double",geo)
    #
    action.hmc_set_force(force, field)
    #
    force*=fg_dt
    field+=force
    force*=1/fg_dt
    #
    action.hmc_set_force(force, field)
    #
    force*=-dt
    momentum+=force
    force*=-1/dt
    # Save the Fourier transform of the new momentum
    momentum_ft.set_complex_from_double(momentum)
    momentum_ft *= 1/V**0.5
    momentum_ft @= fft*momentum_ft
    return force

@q.timer_verbose
def hmc_evolve(field, momentum_ft, field_ft, action, masses, steps, dt, V, fft, ifft):
    # Evolve the field according to the given action using the force 
    # gradient algorithm
    lam = 0.5 * (1.0 - 1.0 / m.sqrt(3.0));
    theta = (2.0 - m.sqrt(3.0)) / 48.0;
    ttheta = theta * dt * dt * dt;
    # The Fourier transformed field is updated, and then the field is 
    # updated based on the new Fourier transformed field
    action.hmc_field_evolve(field_ft, momentum_ft, masses, lam * dt)
    # Perform the inverse Fourier transform
    field.set_double_from_complex(ifft*field_ft)
    field*=1/V**0.5
    #
    force_ft = q.Field("Complex", field.geo())
    # Save a list of the vacuum expectation value of sigma at each point
    # on the trajectory
    vev = []
    for i in range(steps):
        sm_evolve(momentum_ft, field, action, 4.0 * ttheta / dt, 0.5 * dt, fft, ifft);
        action.hmc_field_evolve(field_ft, momentum_ft, masses, (1.0 - 2.0 * lam) * dt);
        field.set_double_from_complex(ifft*field_ft)
        field*=1/V**0.5
        force = sm_evolve(momentum_ft, field, action, 4.0 * ttheta / dt, 0.5 * dt, fft, ifft);
        if i < steps - 1:
            action.hmc_field_evolve(field_ft, momentum_ft, masses, 2.0 * lam * dt);
            field.set_double_from_complex(ifft*field_ft)
            field*=1/V**0.5
        else:
            action.hmc_field_evolve(field_ft, momentum_ft, masses, lam * dt);
            field.set_double_from_complex(ifft*field_ft)
            field*=1/V**0.5
        field_sum = field.glb_sum()
        vev.append(field_sum[0]/V)
    force_ft.set_complex_from_double(force)
    force_ft=fft*force_ft
    force_ft*=1/V**0.5
    # Estimate the masses we should use in order to evolve each field 
    # mode by half of its period
    geo = field_ft.geo()
    masses_new = q.Field("double",geo,geo.multiplicity())
    action.hmc_estimate_mass(masses_new, field_ft, force_ft, np.mean(vev))
    print(f"Vev: {np.mean(vev)}")
    return masses_new

@q.timer_verbose
def run_hmc_evolve(field, momentum_ft, field_ft, action, masses, rs, steps, md_time, V, fft, ifft):
    # Calculate the value of the molecular dynamics Hamiltonian for the 
    # initial field and momentum configuration
    energy = action.hmc_m_hamilton_node(momentum_ft, masses) + action.action_node(field)
    
    # Evolve the field forward in molecular dynamics time using the 
    # given momenta and the Hamiltonian appropriate for the action
    dt = float(md_time) / float(steps)
    masses_new = hmc_evolve(field, momentum_ft, field_ft, action, masses, steps, dt, V, fft, ifft)
    
    # Calculate the change in the value of the molecular dynamics 
    # Hamilton after the evolution 
    delta_h = action.hmc_m_hamilton_node(momentum_ft, masses) + action.action_node(field) - energy;
    
    # Save the new estimated masses
    masses@=masses_new
    
    # Sum over delta_h for every parallel node (each node handles part 
    # of the lattice)
    delta_h = q.glb_sum(delta_h)
    
    return delta_h

@q.timer_verbose
def metropolis_accept(delta_h, traj, rs):
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
            accept_prob = m.exp(-delta_h)
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

@q.timer_verbose
def run_hmc(field, geo, action, masses, traj, rs, fft, ifft):
    # Create a copy of the scalar field
    f0 = field.copy()
    
    # Creates axillary complex fields to store the Fourier 
    # transformations of the scalar field and the momentum field
    field_ft = q.Field("Complex",geo,mult)
    field_ft.set_complex_from_double(field)
    field_ft = fft*field_ft
    field_ft*=1/geo.total_volume()**0.5
    momentum_ft = q.Field("Complex",geo,mult)
    
    # Create a random momentum field, distributed according to a 
    # special distribution that is Gaussian for any given point in
    # the lattice momentum space, but varries in width from point to 
    # point because of the momentum-dependent mass term used for Fourier
    # acceleration
    #momentum = q.Field("double",geo,mult)
    #momentum.set_rand_g(rs.split("set_rand_momentum"), 0.0, 1.0)
    #momentum_ft.set_complex_from_double(momentum)
    #momentum_ft = fft*momentum_ft
    #momentum_ft*=1/geo.total_volume()**0.5
    action.hmc_set_rand_momentum(momentum_ft, masses, rs.split("set_rand_momentum"))
    momentum_ft = ifft*momentum_ft
    momentum_ft*=1/geo.total_volume()**0.5
    momentum = q.Field("double",geo,mult)
    momentum.set_double_from_complex(momentum_ft)
    momentum_ft.set_complex_from_double(momentum)
    momentum_ft = fft*momentum_ft
    momentum_ft*=1/geo.total_volume()**0.5
    
    # The number of steps to take in a single trajectory
    steps = 100
    # The length of a single trajectory in molecular dynamics time
    md_time = 1.0
    
    # Evolve the field over time md_time using the given momenta and 
    # the Hamiltonian appropriate for the given action
    delta_h = run_hmc_evolve(f0, momentum_ft, field_ft, action, masses, rs, steps, md_time, geo.total_volume(), fft, ifft)
    
    # Decide whether to accept or reject the field update using the 
    # metropolis algorithm
    flag, accept_prob = metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    accept_rates.append(accept_prob)
    
    # If the field update is accepted or we are within the first few 
    # trajectories, save the field update
    init_len = 20
    if flag or traj <= init_len:
        q.displayln_info("run_hmc: update field (traj={:d})".format(traj))
        field @= f0

@q.timer_verbose
def test_hmc(total_site, action, mult, n_traj):
    # Create the geometry for the field
    geo = q.Geometry(total_site, mult)
    # Create objects to perform Fourier and inverse Fourier transformations
    fft = q.mk_fft(True)
    ifft = q.mk_fft(False)
    # Save the spacial volume and the total volume of the lattice for 
    # future use
    Vx = total_site[0]*total_site[1]*total_site[2]
    V = Vx*total_site[3]
    
    # Create a random number generator that can be split between 
    # different portions of the lattice
    rs = q.RngState("test_hmc_pions-{}x{}x{}x{}".format(total_site[0], total_site[1], total_site[2], total_site[3]))
    
    # Create the scalar field and set all field values to 1
    field = q.Field("double",geo,mult)
    q.set_unit(field);
    #field.load_double(f"output_data/hmc-pions-sigma-pi-corrs_{total_site[0]}x{total_site[3]}_msq_{m_sq}_lmbd_{lmbd}_alph_{alpha}.field")
    
    # Create a field to store the masses used for Fourier acceleration
    masses = q.Field("double",geo,mult)
    q.set_unit(masses);
    
    # Create the geometry for the axial current field
    geo_cur = q.Geometry(total_site, mult-1)
    # This field will store the calculated axial currents
    axial_current = q.Field("double",geo_cur)
    
    # Stores the index of the current trajectory
    traj = 0
    # The number of trajectories to calculate before taking measurements
    start_measurements = 0;
    
    for i in range(n_traj):
        traj += 1
        
        # Run the HMC algorithm to update the field configuration
        run_hmc(field, geo, action, masses, traj, rs.split("hmc-{}".format(traj)), fft, ifft)
        
        # Calculate the expectation values of phi and phi^2
        q.displayln_info("Average phi^2:")
        psq = phi_squared(field, action)
        q.displayln_info(psq)
        
        q.displayln_info("Average phi:")
        field_sum = field.glb_sum()
        phi=[field_sum[i]/V for i in range(mult)]
        q.displayln_info(phi)
        
        #q.displayln_info("Analytic masses (free case):")
        #ms=[calc_mass([0,0,0,0]),calc_mass([1,0,0,0]),calc_mass([2,0,0,0]),calc_mass([3,0,0,0]),calc_mass([4,0,0,0]),calc_mass([5,0,0,0])]
        #q.displayln_info(ms)
        
        q.displayln_info("Estmiated masses:")
        ms=[masses.get_elem([0,0,0,0],0),masses.get_elem([1,0,0,0],0),masses.get_elem([2,0,0,0],0),masses.get_elem([3,0,0,0],0),masses.get_elem([4,0,0,0],0),masses.get_elem([5,0,0,0],0)]
        q.displayln_info(ms)
        ms=[masses.get_elem([0,0,0,0],1),masses.get_elem([1,0,0,0],1),masses.get_elem([2,0,0,0],1),masses.get_elem([3,0,0,0],1),masses.get_elem([4,0,0,0],1),masses.get_elem([5,0,0,0],1)]
        q.displayln_info(ms)
        
        tslices = field.glb_sum_tslice()
        
        # Calculate the axial current of the current field configuration
        # and save it in axial_current
        action.axial_current_node(axial_current, field)
        tslices_ax_cur = axial_current.glb_sum_tslice()
        
        if i>start_measurements:
            psq_list.append(psq)
            phi_list.append(phi)
            timeslices.append(tslices.to_numpy())
            ax_cur_timeslices.append(tslices_ax_cur.to_numpy())
    
    # Saves the final field configuration so that the next run can be 
    # started where this one left off
    field.save_double(f"output_data/hmc-pions-sigma-pi-corrs_{total_site[0]}x{total_site[3]}_msq_{m_sq}_lmbd_{lmbd}_alph_{alpha}.field")

@q.timer_verbose
def main():
    action = q.ScalarAction(m_sq, lmbd, alpha)
    
    test_hmc(total_site, action, mult, n_traj)

# Stores the average phi^2 for each trajectory
psq_list=[]
# Stores the average values of each phi_i for each trajectory
phi_list=[]
# Stores the timeslice sums of each phi_i for each trajectory
timeslices=[]
# Stores the timeslice sums of the time component of each of the three
# axial currents for each trajectory
ax_cur_timeslices=[]
# Save the acceptance rates
accept_rates=[]

# The lattice dimensions
total_site = [8,8,8,8]

# The multiplicity of the scalar field
mult = 4

# The number of trajectories to calculate
n_traj = 10

# Use action for a Euclidean scalar field. The Lagrangian will be:
# (1/2)*[sum i]|dphi_i|^2 + (1/2)*m_sq*[sum i]|phi_i|^2
#     + (1/24)*lmbd*([sum i]|phi_i|^2)^2
m_sq = -1.0
lmbd = 1.0
alpha = 0.1

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 2, 2, 2],
        [2, 2, 2, 2],
        [2, 2, 2, 4]]

q.begin(sys.argv, size_node_list)

# q.show_machine()

q.qremove_all_info("results")

main()

with open(f"output_data/sigma_pion_corrs_{total_site[0]}x{total_site[3]}_msq_{m_sq}_lmbd_{lmbd}_alph_{alpha}_{datetime.datetime.now().date()}.bin", "wb") as output:
    pickle.dump([accept_rates,psq_list,phi_list,timeslices,ax_cur_timeslices],output)

q.timer_display()

q.end()
