#!/usr/bin/env python3

import sys
import math as m
import numpy as np

import qlat as q

def phi_squared(field):
	geo = field.geo()
	return q.sm_hamilton_node(field)*2/geo.total_volume()/geo.multiplicity()

@q.timer_verbose
def sm_evolve_fg(momentum, field_init, action, fg_dt, dt):
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
def run_hmc_evolve(momentum, field, action, rs, steps, md_time = 1.0):
    energy = q.sm_hamilton_node(momentum) + q.sf_hamilton_node(field, action)
    #
    dt = md_time / steps
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
    #
    delta_h = q.sm_hamilton_node(momentum) + q.sf_hamilton_node(field, action) - energy;
    # Global sum
    delta_h = q.glb_sum(delta_h)
    #
    return delta_h
    return 1

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
    #
    is_reverse_test = traj < 3
    #
    geo = field.geo()
    f0 = field.copy()
    #
    momentum = q.Field("double",geo)
    momentum.set_rand_g(rs.split("set_rand_momentum"), 0.0, 1.0)
    #
    steps = 100
    md_time = 1.0
    #
    delta_h = run_hmc_evolve(momentum, f0, action, rs, steps, md_time)
    #
    #if is_reverse_test:
    #    gm_r = q.GaugeMomentum(geo)
    #    gm_r @= gm
    #    gf0_r = q.GaugeField(geo)
    #    gf0_r @= gf0
    #    delta_h_rev = run_hmc_evolve(gm_r, gf0_r, ga, rs, steps, -md_time)
    #    gf0_r -= gf;
    #    q.displayln_info("run_hmc_evolve reversed delta_diff: {} / {}".format(delta_h + delta_h_rev, delta_h))
    #    q.displayln_info("run_hmc_evolve reversed gf_diff: {} / {}".format(q.qnorm(gf0_r), q.qnorm(gf0)))
    #
    flag, accept_prob = metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    #
    if flag or traj <= 20:
        q.displayln_info("run_hmc: update field (traj={:d})".format(traj))
        field @= f0

@q.timer_verbose
def test_hmc(total_site, action):
    geo = q.Geometry(total_site, 1)
    rs = q.RngState("test_hmc_pions-{}x{}x{}x{}".format(total_site[0], total_site[1], total_site[2], total_site[3]))
    field = q.Field("double",geo,1)
    q.set_unit(field);
    traj = 0
    for i in range(100):
        traj += 1
        run_hmc(field, action, traj, rs.split("hmc-{}".format(traj)))
        print("Average phi^2:")
        psq = phi_squared(field)
        print(psq)
        a.append(psq)

@q.timer_verbose
def main():
	# Set the lattice dimensions
    total_site = [8, 8, 8, 8]
    
    # Use action for a Euclidean scalar field. The Lagrangian will be: 
    # (1/2)*[sum fields]|dphi|^2 + (1/2)*m_sq*[sum fields]|phi|^2 
    #     + (1/24)*lmbd*([sum fields]|phi|^2)^2
    m_sq = -5.0
    lmbd = 1.0
    action = q.ScalarAction(m_sq, lmb)
    
    test_hmc(total_site, action)

a=[]

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

q.timer_display()

print(a)

q.end()
