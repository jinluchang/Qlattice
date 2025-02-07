#!/usr/bin/env python3

import sys
import numpy as np
import pickle
import datetime
import glob
import fnmatch
#import matplotlib.pyplot as plt

import qlat as q

class HMC:
    def __init__(self, action, action_id, total_site, mult, steps, init_length, date, version, fresh_start):
        self.action = action
        self.Vx = total_site[0]*total_site[1]*total_site[2]
        self.V = self.Vx*total_site[3]
        self.total_site = total_site
        self.mult = mult
        self.steps = steps
        self.steps_initial = steps*10
        # The number of trajectories to calculate before doing accept/reject step
        self.init_length = init_length
        
        self.fileid = f"{self.total_site[0]}x{self.total_site[3]}_{action_id}_{date}_{version}"
        self.fileidwc = f"{self.total_site[0]}x{self.total_site[3]}_{action_id}_*_{version}"
        self.fileidwc_no_version = f"{self.total_site[0]}x{self.total_site[3]}_{action_id}_*"
        
        self.traj = 1
        self.accept_prob = 0.0
        
        geo = q.Geometry(total_site)
        # Create a random number generator that can be split between
        # different portions of the lattice
        self.rs = q.RngState(f"false_vacuum_decay-{self.fileidwc}")
        
        # Create the scalar field
        self.field = q.Field(q.ElemTypeRealD, geo, mult)
        self.momentum = q.Field(q.ElemTypeRealD, geo, mult)
        self.force = q.Field(q.ElemTypeRealD, geo, mult)
        # Create an auxiliary field to store the field as it evolves
        self.f0 = q.Field(q.ElemTypeRealD, geo, mult)
        self.aux = q.Field(q.ElemTypeRealD, geo, mult)
        
        if(fresh_start):
            self.init_field()
        else:
            self.load_field()
        
        self.perform_metro = False
    
    @q.timer_verbose
    def run_traj(self):
        if(self.traj<self.init_length):
            self.perform_metro = False
            self.run_hmc(self.rs.split(f"hmc-{self.traj}"))
        else:
            self.perform_metro = True
            self.run_hmc(self.rs.split(f"hmc-{self.traj}"))
        self.traj += 1
    
    def load_field(self):
        filename, self.traj = self.find_latest_traj(f"output_data/fields/fvd_traj_*_{self.fileidwc}.field")
        # if(filename==""):
        #    filename, self.traj = self.find_latest_traj(f"output_data/fields/fvd_traj_*_{self.fileidwc_no_version}.field")
        self.init_length+=self.traj-1
        if(not filename==""):
            self.field.load_double(filename)
            # If a field is loaded, avoid running any trajectories 
            # without a metropolis accept/reject step
            self.init_length = 0
        else:
            self.init_field()
    
    def init_field(self):
        self.field.set_unit()
        self.field *= -3.0
    
    def save_field(self):
        self.field.save_double(f"output_data/fields/fvd_traj_{self.traj}_{self.fileid}.field")
    
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
        if(fnmatch.fnmatch(files[trajnos.index(traj)], f"output_data/fields/fvd_traj_*_{self.fileidwc}.field")):
            traj_new = traj
        else:
            traj_new = 1
        return files[trajnos.index(traj)], traj_new
    
    @q.timer_verbose
    def run_hmc(self, rs):
        # Create a copy of the scalar field
        self.f0 @= self.field
        
        # Create a random momentum field, distributed according to a
        # special distribution that is Gaussian for any given point in
        # the lattice momentum space, but varries in width from point to
        # point because of the momentum-dependent mass term used for Fourier
        # acceleration
        self.action.hmc_set_rand_momentum(self.momentum, rs.split("set_rand_momentum"))
        #q.displayln_info(f"Initial momentum:")
        #q.displayln_info([self.momentum.get_elem_xg(q.Coordinate([0,0,0,0]),0)[:].item(),
        #        self.momentum.get_elem_xg(q.Coordinate([0,0,0,5]),0)[:].item(),
        #        self.momentum.get_elem_xg(q.Coordinate([0,0,0,10]),0)[:].item(),
        #        self.momentum.get_elem_xg(q.Coordinate([0,0,0,15]),0)[:].item(),
        #        self.momentum.get_elem_xg(q.Coordinate([0,0,0,20]),0)[:].item()])
        
        # Evolve the field over time md_time using the given momenta and
        # the Hamiltonian appropriate for the given action
        delta_h = self.run_hmc_evolve(self.f0, self.momentum, rs)

        # Decide whether to accept or reject the field update using the
        # metropolis algorithm
        flag, self.accept_prob = self.metropolis_accept(delta_h, self.traj, rs.split("metropolis_accept"))

        # If the field update is accepted or we are within the first few
        # trajectories, save the field update
        if flag or not self.perform_metro:
            q.displayln_info(f"run_hmc: update field (traj={self.traj})")
            self.field @= self.f0
    
    @q.timer_verbose
    def run_hmc_evolve(self, field, momentum, rs):
        # Calculate the value of the molecular dynamics Hamiltonian for the
        # initial field and momentum configuration
        energy = self.action.hmc_m_hamilton_node(momentum) + self.action.action_node(field)
        
        # Evolve the field forward in molecular dynamics time using the
        # given momenta and the Hamiltonian appropriate for the action
        self.hmc_evolve(field, momentum)
        
        # Calculate the change in the value of the molecular dynamics
        # Hamilton after the evolution
        delta_h = self.action.hmc_m_hamilton_node(momentum) + self.action.action_node(field) - energy;
        
        # Sum over delta_h for every parallel node (each node handles part
        # of the lattice)
        delta_h = q.glb_sum(delta_h)
        
        return delta_h
    
    @q.timer_verbose
    def hmc_evolve(self, field, momentum):
        # Evolve the field according to the given action using the force
        # gradient algorithm
        if(not self.perform_metro):
            steps=self.steps_initial
        else:
            steps=self.steps
        lam = 0.5 * (1.0 - 1.0 / 3.0**0.5);
        theta = (2.0 - 3.0**0.5) / 48.0;
        dt = 1/steps
        ttheta = theta * dt * dt * dt;
        # The field is updated, and then the momentum is updated based 
        # on the new field, and finally the field is updated again
        self.action.hmc_field_evolve(field, momentum, lam*dt)
        for i in range(steps):
            self.sm_evolve(field, momentum, 4.0 * ttheta / dt, 0.5 * dt)
            self.action.hmc_field_evolve(field, momentum, (1.0 - 2.0 * lam) * dt)
            self.sm_evolve(field, momentum, 4.0 * ttheta / dt, 0.5 * dt)
            if i < steps - 1:
                self.action.hmc_field_evolve(field, momentum, 2.0 * lam * dt)
            else:
                self.action.hmc_field_evolve(field, momentum, lam * dt)
    
    @q.timer_verbose
    def sm_evolve(self, field_init, momentum, fg_dt, dt):
        # Evolve the momentum field according to the given action using the
        # force gradient algorithm
        self.aux @= field_init
        self.action.hmc_set_force(self.force, self.aux)
        
        self.force *= fg_dt
        self.aux += self.force
        self.force *= 1/fg_dt
        
        self.action.hmc_set_force(self.force, self.aux)
        
        self.force *= -dt
        momentum += self.force
    
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

class Measurements:
    def __init__(self, total_site, actions, save_file):
        self.save_file = save_file
        # Stores the trajectory number for debugging purposes
        self.trajs = []
        # Save the acceptance rates
        self.accept_rates=[]
        # Stores the average phi^2 for each trajectory
        self.psq_list=[]
        # Stores the average values of each phi_i for each trajectory
        self.phi_list=[]
        self.timeslices=[]
        #
        self.delta_actions = {}
        self.actions = {}
        for k in actions:
            self.actions[k] = {}
            self.delta_actions[k] = {}
            for a in actions[k]:
                self.actions[k][a] = actions[k][a]
                self.delta_actions[k][a] = []
    
    def measure(self, hmc, timeslices=False):
        self.trajs.append(hmc.traj)
        self.accept_rates.append(hmc.accept_prob)
        
        # Calculate the expectation values of phi^2
        self.psq_list.append(self.calc_psq(hmc.field, hmc.action))
        # Calculate the expectation value of phi
        field_sum = hmc.field.glb_sum()[0]
        self.phi_list.append([field_sum[i]/hmc.field.geo.total_volume for i in range(hmc.mult)])
        if(timeslices):
            self.timeslices.append(hmc.field.glb_sum_tslice().to_numpy())
        #
        for k in self.actions:
            for a in self.actions[k]:
                dS = hmc.action.action_node(hmc.field) - self.actions[k][a].action_node(hmc.field)
                self.delta_actions[k][a].append(q.glb_sum(dS))
    
    def display_measurements(self):
        q.displayln_info("Average phi:")
        q.displayln_info(self.phi_list[-1])
        q.displayln_info("Average phi^2:")
        q.displayln_info(self.psq_list[-1])
    
    def plot_measurements(self):
        pass
    
    def get_representatives(self, field):
        return [field.get_elem_xg(q.Coordinate([0,0,0,i]),0)[:].item() for i in range(0,10,2)]
    
    def calc_psq(self, field, action):
        # Calculate the average value of phi^2
        phi_sq = action.sum_sq(field) # Returns sum of field^2/2
        phi_sq = q.glb_sum(phi_sq) # Sums over all nodes
        return phi_sq/field.geo.total_volume
    
    def save(self):
        with open(self.save_file, "wb") as output:
            pickle.dump({"trajs": self.trajs,
                "accept_rates": self.accept_rates,
                "psq_list": self.psq_list,
                "phi_list": self.phi_list,
                "timeslices": self.timeslices,
                "delta_actions": self.delta_actions}, output)
    
    def load(self):
        if len(glob.glob(self.save_file)):
            with open(self.save_file,"rb") as input:
                q.displayln_info(f"Loading old measurements from {self.save_file}")
                data = pickle.load(input)
                self.trajs.extend(data["trajs"])
                self.accept_rates.extend(data["accept_rates"])
                self.psq_list.extend(data["psq_list"])
                self.phi_list.extend(data["phi_list"])
                self.timeslices.extend(data["timeslices"])
                for k in data["delta_actions"]:
                    if k in self.delta_actions:
                        for a in data["delta_actions"][k]:
                            if a in self.delta_actions[k]:
                                self.delta_actions[k][a].extend(data["delta_actions"][k][a])
                            else:
                                self.delta_actions[k][a] = data["delta_actions"][k][a]
                    else:
                        self.delta_actions[k] = data["delta_actions"][k]

@q.timer_verbose
def main():
    # The lattice dimensions
    Nt = 100
    total_site = q.Coordinate([1,1,1,Nt])
    # The multiplicity of the field
    mult = 1
    alpha = 1.0
    beta = 9.0
    FV_offset = 0.3
    TV_offset = 0.0
    barrier_strength = 1000.0
    M = 1.0
    L = 0.0
    P = 0.0
    epsilon = 1e-10
    t_full1 = 5
    t_full2 = 5
    t_TV = 30
    t_FV_mid = 6
    dt = 0.2
    # The number of trajectories to calculate
    n_traj = 50000
    proj_sq = False
    #
    version = "11-0"
    date = datetime.datetime.now().date()
    # The number of steps to take in a single trajectory
    steps = 10
    #
    init_length = 200
    fresh_start = False
    # The number of trajectories to run before each save
    save_frequency = 5000
    
    for i in range(1,len(sys.argv)):
        try:
            if(sys.argv[i]=="-a"):
                alpha = float(sys.argv[i+1])
            elif(sys.argv[i]=="-b"):
                beta = float(sys.argv[i+1])
            elif(sys.argv[i]=="-o"):
                FV_offset = float(sys.argv[i+1])
            elif(sys.argv[i]=="-O"):
                TV_offset = float(sys.argv[i+1])
            elif(sys.argv[i]=="-B"):
                barrier_strength = float(sys.argv[i+1])
            elif(sys.argv[i]=="-f"):
                t_full1 = int(sys.argv[i+1])
                t_full2 = int(sys.argv[i+1])
            elif(sys.argv[i]=="-f1"):
                t_full1 = int(sys.argv[i+1])
            elif(sys.argv[i]=="-f2"):
                t_full2 = int(sys.argv[i+1])
            elif(sys.argv[i]=="-t"):
                t_TV = int(sys.argv[i+1])
            elif(sys.argv[i]=="-m"):
                t_FV_mid = int(sys.argv[i+1])
            elif(sys.argv[i]=="-M"):
                M = float(sys.argv[i+1])
            elif(sys.argv[i]=="-L"):
                L = float(sys.argv[i+1])
            elif(sys.argv[i]=="-P"):
                P = float(sys.argv[i+1])
            elif(sys.argv[i]=="-e"):
                epsilon = float(sys.argv[i+1])
            elif(sys.argv[i]=="-p"):
                proj_sq=True
            elif(sys.argv[i]=="-D"):
                a = sys.argv[i+1].split("x")
                total_site = q.Coordinate([int(a[j]) for j in range(4)])
                Nt = total_site[3]
            elif(sys.argv[i]=="-d"):
                dt = float(sys.argv[i+1])
            elif(sys.argv[i]=="-T"):
                n_traj = int(sys.argv[i+1])
            elif(sys.argv[i]=="-s"):
                steps = int(sys.argv[i+1])
            elif(sys.argv[i]=="-S"):
                save_frequency = int(sys.argv[i+1])
            elif(sys.argv[i]=="-R"):
                fresh_start = True
            elif(sys.argv[i]=="-i"):
                init_length = int(sys.argv[i+1])
        except:
            raise Exception("Invalid arguments: use \
                            -a for alpha, \
                            -b for beta, \
                            -o to offset the barrier location in V_FV, \
                            -O to offset the barrier location in V_TV, \
                            -B for the barrier strength used in H_FV and H_TV, \
                            -f for the time to evolve with H_full (both first and second time), \
                            -f1 for the time to evolve with H_full the first time, \
                            -f2 for the time to evolve with H_full the second time, \
                            -t for the time to evolve with H_TV, \
                            -m for the time to evolve with H_FV_mid, \
                            -M to set the left barrier for H_TV, \
                            -L to set the right barrier for H_TV, \
                            -P to set the projection strength, \
                            -e to set epsilon (for V_proj), \
                            -p use V_bar^2 instead of V_bar in V_proj, \
                            -D for lattice dimensions, \
                            -d for the lattice spacing, \
                            -T for number of trajectories, \
                            -s for the number of steps in a trajectory, \
                            -S for the number of trajectories to run before each save, \
                            -R to force restarting with blank initial field, \
                            -i for the number of trajectories to do at the beginning without a Metropolis step.")
    
    t_FV_out = int((Nt - t_full1 - t_full2 - t_TV - t_FV_mid)/2)
    t_FV_mid = Nt - t_full1 - t_full2 - t_TV - 2*t_FV_out
    
    action = q.QMAction(alpha, beta, FV_offset, TV_offset, barrier_strength, M, L, P, epsilon, t_full1, t_full2, t_FV_out, t_FV_mid, 0, dt, proj_sq)
    hmc = HMC(action,f"alpha_{alpha}_beta_{beta}_dt_{dt}_FVoff_{FV_offset}_TVoff_{TV_offset}_bar_{barrier_strength}_M_{M}_L_{L}_P_{P}_eps_{epsilon}_tfull1_{t_full1}_tfull2_{t_full2}_proj2_{proj_sq}_tTV_{t_TV}_tFV_{t_FV_out*2+t_FV_mid}_tFVout_{t_FV_out}_tFVmid_{t_FV_mid}",total_site,mult,steps,init_length,date,version,fresh_start)
    
    steps = np.array([0.001, 0.002, 0.003, 0.004, 0.006, 0.008, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06, 0.07, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]) #np.concatenate([0.001*np.arange(1,100), 0.1 + 0.01*np.arange(0,91)])
    measure_Ms = steps[steps>M] #[round(min(max(M,0.001)*2**i, 1.0),5) for i in range(1,10)]
    measure_Ls = steps[steps>L] #[round(min(max(L,0.001)*2**i, 1.0),5) for i in range(1,10)]
    measure_Ps = steps[steps>P] #[round(min(max(L,0.001)*2**i, 1.0),5) for i in range(1,10)]
    
    actions = {"M": {}, "L": {}, "P": {}}
    if(L==1.0):
        for Mi in measure_Ms:
            actions["M"][f"{Mi}"] = q.QMAction(alpha, beta, FV_offset, TV_offset, barrier_strength, Mi, L, P, epsilon, t_full1, t_full2, t_FV_out, t_FV_mid, 0, dt, proj_sq)
    if(M==1.0):
        for Li in measure_Ls:
            actions["L"][f"{Li}"] = q.QMAction(alpha, beta, FV_offset, TV_offset, barrier_strength, M, Li, P, epsilon, t_full1, t_full2, t_FV_out, t_FV_mid, 0, dt, proj_sq)
    if(M==1.0 and L==0.0):
        for Pi in measure_Ps:
            actions["P"][f"{Pi}"] = q.QMAction(alpha, beta, FV_offset, TV_offset, barrier_strength, M, L, Pi, epsilon, t_full1, t_full2, t_FV_out, t_FV_mid, 0, dt, proj_sq)
    
    measurements = Measurements(total_site, actions, f"output_data/measurements_{hmc.fileid}.bin")
    
    # If observables have been saved from a previous calculation (on the
    # same day), then load that file first
    measurements.load()
    
    while hmc.traj <= n_traj:
        # Run the HMC algorithm to update the field configuration
        hmc.run_traj()
        measurements.measure(hmc, timeslices=(hmc.traj%1 == 0))
        measurements.display_measurements()
        q.displayln_info(f"{hmc.fileid}")
        #if hmc.traj%10 == 0:
        #    measurements.plot_measurements()
        
        #x = np.arange(-5,5,0.1)
        #for t in [0,t_full1,t_full1+t_FV_out,t_full1+t_FV_out+t_FV_mid,t_full1+t_FV_out*2+t_FV_mid]: #,t_full1+t_full2+t_FV_out*2+t_FV_mid]:
        #    plt.plot([min((action.V(i,t) - action.V(0,t))*Nt/20.0, 10.0) + t for i in x],x)
        #plt.plot(range(Nt), np.mean(measurements.timeslices,axis=0))
        #plt.show()
        
        if hmc.traj%save_frequency == 0:
            hmc.save_field()
            measurements.save()
    
    # Saves the final field configuration so that the next run can be
    # started where this one left off
    hmc.save_field()
    measurements.save()
    
    q.displayln_info(f"Acceptance rate: {np.mean(measurements.accept_rates[-int(n_traj/2):])}")
    for da in measurements.delta_actions["M"]:
        q.displayln_info(f"e^(Delta S) for M={da}: {np.mean(np.exp(measurements.delta_actions['M'][da][-int(n_traj/2):]))}+-{np.std(np.exp(measurements.delta_actions['M'][da][-int(n_traj/2):]))/(n_traj/2)**0.5}")
    for da in measurements.delta_actions["L"]:
        q.displayln_info(f"e^(Delta S) for L={da}: {np.mean(np.exp(measurements.delta_actions['L'][da][-int(n_traj/2):]))}+-{np.std(np.exp(measurements.delta_actions['L'][da][-int(n_traj/2):]))/(n_traj/2)**0.5}")
    for da in measurements.delta_actions["P"]:
        q.displayln_info(f"e^(Delta S) for P={da}: {np.mean(np.exp(measurements.delta_actions['P'][da][-int(n_traj/2):]))}+-{np.std(np.exp(measurements.delta_actions['P'][da][-int(n_traj/2):]))/(n_traj/2)**0.5}")
    
    q.displayln_info(f"CHECK: The vacuum expectation value of phi_0 is {round(np.mean(measurements.phi_list[int(n_traj/2):], axis=0)[0],2)}.")
    q.displayln_info(f"CHECK: The vacuum expectation value of phi^2 is {round(np.mean(measurements.psq_list[int(n_traj/2):]),2)}.")

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [1, 1, 1, 16],
        [1, 1, 1, 32],
        [1, 1, 1, 64]]

q.begin_with_mpi(size_node_list)

q.show_machine()

q.qremove_all_info("results")

main()

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: Simulation finished successfully.")
