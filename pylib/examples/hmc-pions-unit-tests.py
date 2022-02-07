#!/usr/bin/env python3

import sys
import math as m
import numpy as np

import qlat as q

def phi_squared(field):
	geo = field.geo()
	return q.sm_hamilton_node(field)*2/geo.total_volume()/geo.multiplicity()

def main():
    rs = q.RngState("hmc-pions-unit-tests")
    
	# Creates some fields to use for tests
    unit_geo = q.Geometry([3,3,1,1], 1)
    unit_field = q.Field("double",unit_geo,1)
    q.set_unit(unit_field)
    rand_geo = q.Geometry([24,24,24,24], 1)
    rand_field = q.Field("double",rand_geo,1)
    rand_field.set_rand_g(rs.split("set_rand_g field"), 0.0, 2.0)
    checkers_geo = q.Geometry([6,6,6,6], 1)
    checkers_field = q.Field("double",checkers_geo,1)
    checkers_field.set_checkers()
    
    q.displayln_info("Average phi^2 of the unit field (should be 1): ")
    q.displayln_info(phi_squared(unit_field))
    
    q.displayln_info("Average phi^2 of the checkers field (should be 1): ")
    q.displayln_info(phi_squared(checkers_field))
    
    q.displayln_info("Average phi^2 of the random field (should be 4): ")
    q.displayln_info(phi_squared(rand_field))
    
    q.displayln_info("Average phi of the unit field (should be 1): ")
    q.displayln_info(sum(unit_field.glb_sum())/unit_geo.total_volume()/unit_geo.multiplicity())
    
    q.displayln_info("Average phi of the checkers field (should be 0): ")
    q.displayln_info(sum(checkers_field.glb_sum())/checkers_geo.total_volume()/checkers_geo.multiplicity())
    
    q.displayln_info("Average phi of the random field (should be 0): ")
    q.displayln_info(sum(rand_field.glb_sum())/rand_geo.total_volume()/rand_geo.multiplicity())
    
    # Creates some actions to use for tests
    action0 = q.ScalarAction(0.0, 0.0)
    action1 = q.ScalarAction(1.0, 0.0)
    action2 = q.ScalarAction(1.0, 24.0)
    
    q.displayln_info("Average derivative squared term for unit field (should be 0): ")
    h = q.sf_hamilton_node(unit_field, action0)
    h = q.glb_sum(h)
    q.displayln_info(h/unit_geo.total_volume()/unit_geo.multiplicity())
    
    q.displayln_info("Average derivative squared term for checkers field (should be dim*dim/2=8: ")
    h = q.sf_hamilton_node(checkers_field, action0)
    h = q.glb_sum(h)
    q.displayln_info(h/checkers_geo.total_volume()/checkers_geo.multiplicity())
    
    q.displayln_info("Average derivative squared term for random field (should be dim*sigma^2=16): ")
    h = q.sf_hamilton_node(rand_field, action0)
    h = q.glb_sum(h)
    q.displayln_info(h/rand_geo.total_volume()/rand_geo.multiplicity())
    
    q.displayln_info("Average free Hamiltonian for unit field (should be 0+1/2): ")
    h = q.sf_hamilton_node(unit_field, action1)
    h = q.glb_sum(h)
    q.displayln_info(h/unit_geo.total_volume()/unit_geo.multiplicity())
    
    q.displayln_info("Average free Hamiltonian for checkers field (should be 8+1/2: ")
    h = q.sf_hamilton_node(checkers_field, action1)
    h = q.glb_sum(h)
    q.displayln_info(h/checkers_geo.total_volume()/checkers_geo.multiplicity())
    
    q.displayln_info("Average free Hamiltonian for random field (should be 16+sigma^2/2=18): ")
    h = q.sf_hamilton_node(rand_field, action1)
    h = q.glb_sum(h)
    q.displayln_info(h/rand_geo.total_volume()/rand_geo.multiplicity())
    
    q.displayln_info("Average Hamiltonian for unit field (should be 0+1/2+1): ")
    h = q.sf_hamilton_node(unit_field, action2)
    h = q.glb_sum(h)
    q.displayln_info(h/unit_geo.total_volume()/unit_geo.multiplicity())
    
    q.displayln_info("Average Hamiltonian for checkers field (should be 8+1/2+1: ")
    h = q.sf_hamilton_node(checkers_field, action2)
    h = q.glb_sum(h)
    q.displayln_info(h/checkers_geo.total_volume()/checkers_geo.multiplicity())
    
    q.displayln_info("Average Hamiltonian for random field (should be 16+2+3*sigma^4=66): ")
    h = q.sf_hamilton_node(rand_field, action2)
    h = q.glb_sum(h)
    q.displayln_info(h/rand_geo.total_volume()/rand_geo.multiplicity())

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 2, 2, 2],
        [2, 2, 2, 2],
        [2, 2, 2, 4]]

q.begin(sys.argv, size_node_list)

q.qremove_all_info("results")

main()

q.end()
