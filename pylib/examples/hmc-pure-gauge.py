#!/usr/bin/env python3

import sys
import qlat as q
from math import sqrt


def gm_evolve_fg(gm, gf_init, ga, fg_dt, dt):
    geo = gf_init.geo()
    gf = q.Field("ColorMatrix", geo, 4)
    gf @= gf_init
    #
    gm_force = q.Field("ColorMatrix", geo, 4)
    #
    q.set_gm_force(gm_force, gf, ga)
    #
    q.gf_evolve(gf, gm_force, fg_dt)
    #
    q.set_gm_force(gm_force, gf, ga)
    #
    gm_force *= dt
    gm += gm_force
    

def run_hmc_evolve(gm, gf, ga, rs, steps, md_time = 1.0):
    energy = q.gm_hamilton_node(gm) + q.gf_hamilton_node(gf, ga)
    #
    dt = md_time / steps
    lam = 0.5 * (1.0 - 1.0 / sqrt(3.0));
    theta = (2.0 - sqrt(3.0)) / 48.0;
    ttheta = theta * dt * dt * dt;
    #
    q.gf_evolve(gf, gm, lam * dt)
    #
    for i in range(steps):
        gm_evolve_fg(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt);
        q.gf_evolve(gf, gm, (1.0 - 2.0 * lam) * dt);
        gm_evolve_fg(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt);
        if i < steps - 1:
            q.gf_evolve(gf, gm, 2.0 * lam * dt);
        else:
            q.gf_evolve(gf, gm, lam * dt);
    #
    q.unitarize(gf)
    #
    delta_h = q.gm_hamilton_node(gm) + q.gf_hamilton_node(gf, ga) - energy;
    delta_h = q.glb_sum(delta_h)
    #
    return delta_h

def run_hmc(gf, ga, traj, rs):
    geo = gf.geo()
    gf0 = q.Field("ColorMatrix", geo, 4)
    gf0 @= gf
    #
    gm = q.Field("ColorMatrix", geo, 4)
    q.set_rand_gauge_momentum(gm, 1.0, rs.split("set_rand_gauge_momentum"))
    #
    steps = 6
    md_time = 1.0
    #
    delta_h = run_hmc_evolve(gm, gf0, ga, rs, steps, md_time)
    #
    # TODO metropolis_accept
    #
    gf @= gf0

def test_hmc(total_site, ga):
    geo = q.Geometry(total_site, 1)
    rs = q.RngState("test_hmc-{}x{}x{}x{}".format(total_site[0], total_site[1], total_site[2], total_site[3]))
    gf = q.Field("ColorMatrix", geo, 4)
    q.set_unit(gf);
    traj = 0
    for i in range(30):
        traj += 1
        run_hmc(gf, ga, traj, rs.split("hmc-{}".format(traj)))
        plaq_avg = q.gf_avg_plaq(gf)
        if q.get_id_node() == 0:
            print(traj, plaq_avg)

q.begin(sys.argv, [
    (1, 1, 1, 1),
    (1, 1, 1, 2),
    (1, 1, 1, 4),
    (1, 1, 1, 8),
    (2, 2, 2, 2),
    (2, 2, 2, 4)])

print("id_node: {:4} / {} ; coor_node: {:9} / {}".format(
    q.get_id_node(),
    q.get_num_node(),
    str(q.get_coor_node()),
    str(q.get_size_node())))

total_site = (4, 4, 4, 8)

ga = q.GaugeAction(2.13, -0.331)

test_hmc(total_site, ga)

ga = q.GaugeAction(5.5, 0.0)

test_hmc(total_site, ga)

q.timer_display()

q.end()
