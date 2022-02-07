#!/usr/bin/env python3

import sys
import math as m
import numpy as np

import qlat as q

def phi_squared(field):
    # Calculate the average value of phi^2
    phi_sq = q.sm_hamilton_node(field) # Returns sum of field^2/2
    phi_sq = q.glb_sum(phi_sq) # Sums over all nodes
    geo = field.geo()
    return phi_sq*2/geo.total_volume()/geo.multiplicity()

@q.timer_verbose
def sm_evolve_fg(momentum, field_init, action, fg_dt, dt):
    # Evolve the momentum field according to the given action using the  
    # force gradient algorithm
    geo = field_init.geo()
    field = q.Field("double",geo)
    field @= field_init
    #
    force = q.Field("double",geo)
    #
    q.set_sm_force(force, field, action)
    #
    q.sf_evolve(field, force, fg_dt)
    #
    q.set_sm_force(force, field, action)
    #
    # q.display_gm_force_magnitudes(gm_force, 5)
    # q.displayln_info(q.get_gm_force_magnitudes(gm_force, 5))
    #
    force *= dt
    momentum -= force

@q.timer_verbose
def sm_evolve(momentum, field, action, dt):
    # Evolve the momentum field according to the given action using the  
    # leapfrog algorithm
    geo = field.geo()
    #
    force = q.Field("double",geo)
    #
    q.set_sm_force(force, field, action)
    #
    # q.display_gm_force_magnitudes(gm_force, 5)
    # q.displayln_info(q.get_gm_force_magnitudes(gm_force, 5))
    #
    force *= dt
    momentum -= force

@q.timer_verbose
def hmc_evolve(momentum, field, action, steps, dt):
    # Evolve the field according to the given action using the force 
    # gradient algorithm
    lam = 0.5 * (1.0 - 1.0 / m.sqrt(3.0));
    theta = (2.0 - m.sqrt(3.0)) / 48.0;
    ttheta = theta * dt * dt * dt;
    #
    q.sf_evolve(field, momentum, lam * dt)
    #
    for i in range(steps):
        sm_evolve_fg(momentum, field, action, 4.0 * ttheta / dt, 0.5 * dt);
        q.sf_evolve(field, momentum, (1.0 - 2.0 * lam) * dt);
        sm_evolve_fg(momentum, field, action, 4.0 * ttheta / dt, 0.5 * dt);
        if i < steps - 1:
            q.sf_evolve(field, momentum, 2.0 * lam * dt);
        else:
            q.sf_evolve(field, momentum, lam * dt);
    return momentum

@q.timer_verbose
def hmc_evolve_leapfrog(momentum, field, action, steps, dt):
    # Evolve the field according to the given action using the leapfrog 
    # algorithm
    for i in range(steps):
        sm_evolve(momentum, field, action, 0.5 * dt)
        q.sf_evolve(field, momentum, dt)
        sm_evolve(momentum, field, action, 0.5 * dt)
    return momentum

@q.timer_verbose
def run_hmc_evolve(momentum, field, action, rs, steps, md_time = 1.0):
    # Calculate the value of the molecular dynamics Hamiltonian for the 
    # initial field and momentum configuration
    energy = q.sm_hamilton_node(momentum) + q.sf_hamilton_node(field, action)
    
    # Evolve the field forward in molecular dynamics time using the 
    # given momenta and the Hamiltonian appropriate for the action
    dt = float(md_time) / float(steps)
    momentum = hmc_evolve(momentum, field, action, steps, dt)
    
    # Calculate the change in the value of the molecular dynamics 
    # Hamilton after the evolution 
    delta_h = q.sm_hamilton_node(momentum) + q.sf_hamilton_node(field, action) - energy;
    
    # Sum over delta_h for every parallel node (each node handles part 
    # of the lattice)
    delta_h = q.glb_sum(delta_h)
    
    return delta_h

@q.timer_verbose
def metropolis_accept(delta_h, traj, rs):
    flag_d = 0.0
    accept_prob = 0.0
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
    flag = flag_d > 0.5
    q.displayln_info("metropolis_accept: flag={:d} with accept_prob={:.1f}% delta_h={:.16f} traj={:d}".format(
        flag, accept_prob * 100.0, delta_h, traj))
    return flag, accept_prob

@q.timer_verbose
def run_hmc(field, action, traj, rs):
    # Get the field geometry
    geo = field.geo()
    # Create a copy of the scalar field
    f0 = field.copy()
    
    # Create a random momentum field, distributed according to a 
    # gaussian distribution centered at 0.0 with sigma=1.0
    momentum = q.Field("double",geo)
    momentum.set_rand_g(rs.split("set_rand_momentum"), 0.0, 1.0)
    
    # The number of steps to take in a single trajectory
    steps = 100
    # The length of a single trajectory in molecular dynamics time
    md_time = 1.0
    
    # Evolve the field over time md_time using the given momenta and 
    # the Hamiltonian appropriate for the given action
    delta_h = run_hmc_evolve(momentum, f0, action, rs, steps, md_time)
    
    # Test reversibility
    #if traj < 3:
    #    gm_r = q.GaugeMomentum(geo)
    #    gm_r @= gm
    #    gf0_r = q.GaugeField(geo)
    #    gf0_r @= gf0
    #    delta_h_rev = run_hmc_evolve(gm_r, gf0_r, ga, rs, steps, -md_time)
    #    gf0_r -= gf;
    #    q.displayln_info("run_hmc_evolve reversed delta_diff: {} / {}".format(delta_h + delta_h_rev, delta_h))
    #    q.displayln_info("run_hmc_evolve reversed gf_diff: {} / {}".format(q.qnorm(gf0_r), q.qnorm(gf0)))
    #
    
    # Decide whether to accept or reject the field update using the 
    # metropolis algorithm
    flag, accept_prob = metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    
    # If the field update is accepted or we are within the first few 
    # trajectories, save the field update
    init_len = 0
    if flag or traj <= init_len:
        q.displayln_info("run_hmc: update field (traj={:d})".format(traj))
        field @= f0

@q.timer_verbose
def test_hmc(total_site, action, mult, n_traj):
    # Create the geometry for the field
    geo = q.Geometry(total_site, mult)
    
    # Create a random number generator that can be split between 
    # different portions of the lattice
    rs = q.RngState("test_hmc_pions-{}x{}x{}x{}".format(total_site[0], total_site[1], total_site[2], total_site[3]))
    
    # Create the scalar field and set all field values to 1
    field = q.Field("double",geo,mult)
    q.set_unit(field);
    
    traj = 0
    for i in range(n_traj):
        traj += 1
        
        # Run the HMC algorithm to update the field configuration
        #momentum @= run_hmc(field, momentum, action, traj, rs.split("hmc-{}".format(traj)))
        run_hmc(field, action, traj, rs.split("hmc-{}".format(traj)))
        
        # Calculate the expectation values of phi and phi^2
        q.displayln_info("Average phi^2:")
        psq = phi_squared(field)
        q.displayln_info(psq)
        q.displayln_info("Average phi:")
        phi = sum(field.glb_sum())/geo.total_volume()/geo.multiplicity()
        q.displayln_info(phi)
        if i % 1 == 0:
            psq_list.append(psq)
            phi_list.append(phi)

@q.timer_verbose
def main():
	# The lattice dimensions
    total_site = [8, 8, 8, 8]
    
    # The multiplicity of the scalar field
    mult = 1
    
    # The number of trajectories to calculate
    n_traj = 20
    
    # Use action for a Euclidean scalar field. The Lagrangian will be: 
    # (1/2)*[sum fields]|dphi|^2 + (1/2)*m_sq*[sum fields]|phi|^2 
    #     + (1/24)*lmbd*([sum fields]|phi|^2)^2
    m_sq = -5.0
    lmbd = 24.0/16.0
    action = q.ScalarAction(m_sq, lmbd)
    
    test_hmc(total_site, action, mult, n_traj)

psq_list=[]
phi_list=[]

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

q.displayln_info("Expectation value of phi^2 on all trajectories:")
q.displayln_info(psq_list)
q.displayln_info("Expectation value of phi on all trajectories:")
q.displayln_info(phi_list)

#q.timer_display()

q.end()
