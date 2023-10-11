import matplotlib.pyplot as plt
import pandas as pd
from data_analysis_blocked import *

lambdas_ = []
msqs = []
lattice_sizes_ = []
alphas_ = []
trajs = []
cutoffs = []
block_sizes = []
pion_masses = []
pion_mass_errors = []
pion_mass_ests = []
sigma_masses = []
sigma_mass_errors = []
sigma_mass_ests = []
fpis = []
fpi_errors = []
fpi_ests = []
mpi_from_fpis = []
mpi_from_fpi_errors = []
fpi2s = []
fpi2_errors = []
fpi2_ests = []
mpi_from_fpi2s = []
mpi_from_fpi2_errors = []
sigma_vevs = []
sigma_vev_errors = []
sigma_mass_effs = []
sigma_mass_eff_errors = []

lambdas=[10000.0]
msq_facts=[0.3, 0.2, 0.15, 0.099, 0.095, 0.09, 0.07]
lattice_sizes=[[12, 24]]
alphas=[0.015]

for lmbd in lambdas:
    for msq_fact in msq_facts:
        for lattice_size  in lattice_sizes:
            for alpha in alphas:
                data = LatticeData(lattice_size[0], lattice_size[1], -round(msq_fact*lmbd, 3), lmbd, alpha, "3-1", 2000, 200)
                try:
                    data.load_all_data()
                    #data.cut_old_data()
                except Exception as e:
                    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                    print("ERROR:")
                    print(e)
                    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                    print("Continuing...")
                    continue
                try:
                    lambdas_.append(lmbd)
                    msqs.append(-round(msq_fact*lmbd, 2))
                    lattice_sizes_.append(f"{lattice_size[0]}x{lattice_size[0]}x{lattice_size[0]}x{lattice_size[1]}")
                    alphas_.append(alpha)
                    trajs.append(data.data_len)
                    cutoffs.append(data.cutoff)
                    block_sizes.append(data.block_size)
                except Exception as e:
                    lambdas_.append(0.0)
                    msqs.append(0.0)
                    lattice_sizes_.append("-")
                    alphas_.append(0.0)
                    trajs.append(0)
                    cutoffs.append(0)
                    block_sizes.append(0)
                    print(e)
                try:
                    data.plot_phi()
                    plt.savefig(f"plots/phi_av/phi_av_{lattice_size[0]}x{lattice_size[1]}_msq_{-round(msq_fact*lmbd, 2)}_lmbd_{lmbd}_alph_{alpha}.png",facecolor=(1, 1, 1))
                    plt.clf()
                except Exception as e:
                    print(e)
                #
                try:
                    data.calc_quick_ests()
                    pion_mass_ests.append(data.pion_mass_est)
                    sigma_mass_ests.append(data.sigma_mass_est)
                    fpi_ests.append(data.fpi_est)
                except:
                    pion_mass_ests.append(0.0)
                    sigma_mass_ests.append(0.0)
                    fpi_ests.append(0.0)
                #
                try:
                    data.calc_sigma_vev()
                    sigma_vevs.append(data.vev_sigma)
                    sigma_vev_errors.append(data.vev_sigma_err)
                except Exception as e:
                    sigma_vevs.append(0.0)
                    sigma_vev_errors.append(0.0)
                    print(e)
                try:
                    data.calc_sigma_mass_eff()
                    sigma_mass_effs.append(data.sigma_mass_eff)
                    sigma_mass_eff_errors.append(data.sigma_mass_eff_err)
                except Exception as e:
                    sigma_mass_effs.append(0.0)
                    sigma_mass_eff_errors.append(0.0)
                    print(e)
                try:
                    data.calc_pion_mass()
                    pion_masses.append(data.pion_mass)
                    pion_mass_errors.append(data.pion_mass_err)
                    data.plot_pion_fit()
                    plt.savefig(f"plots/mpi/pion_fit_{lattice_size[0]}x{lattice_size[1]}_msq_{-round(msq_fact*lmbd, 2)}_lmbd_{lmbd}_alph_{alpha}.png",facecolor=(1, 1, 1))
                    plt.clf()
                except Exception as e:
                    pion_masses.append(0.0)
                    pion_mass_errors.append(0.0)
                    print(e)
                try:
                    data.calc_sigma_mass()
                    sigma_masses.append(data.sigma_mass)
                    sigma_mass_errors.append(data.sigma_mass_err)
                    data.plot_sigma_fit()
                    plt.savefig(f"plots/msigma/sigma_fit_{lattice_size[0]}x{lattice_size[1]}_msq_{-round(msq_fact*lmbd, 2)}_lmbd_{lmbd}_alph_{alpha}.png",facecolor=(1, 1, 1))
                    plt.clf()
                except Exception as e:
                    sigma_masses.append(0.0)
                    sigma_mass_errors.append(0.0)
                    print(e)
                try:
                    data.calc_fpi()
                    fpis.append(data.fpi)
                    fpi_errors.append(data.fpi_err)
                    mpi_from_fpis.append(data.mpi_from_fpi)
                    mpi_from_fpi_errors.append(data.mpi_from_fpi_err)
                    data.plot_fpi_fit()
                    plt.savefig(f"plots/fpi/fpi_fit_{lattice_size[0]}x{lattice_size[1]}_msq_{-round(msq_fact*lmbd, 2)}_lmbd_{lmbd}_alph_{alpha}.png",facecolor=(1, 1, 1))
                    plt.clf()
                except Exception as e:
                    fpis.append(0.0)
                    fpi_errors.append(0.0)
                    mpi_from_fpis.append(0.0)
                    mpi_from_fpi_errors.append(0.0)
                    print(e)
                try:
                    data.calc_fpi2()
                    fpi2s.append(data.fpi2)
                    fpi2_errors.append(data.fpi2_err)
                    mpi_from_fpi2s.append(data.mpi_from_fpi2)
                    mpi_from_fpi2_errors.append(data.mpi_from_fpi2_err)
                except Exception as e:
                    fpi2s.append(0.0)
                    fpi2_errors.append(0.0)
                    mpi_from_fpi2s.append(0.0)
                    mpi_from_fpi2_errors.append(0.0)
                    print(e)
                plt.close()

df = pd.DataFrame({
    "lattice size": lattice_sizes_,
    "m^2": msqs,
    "lambda": lambdas_,
    "alpha": alphas_,
    "pion mass": pion_masses,
    "pion mass error": pion_mass_errors,
    "estimated pion mass": pion_mass_ests,
    "sigma mass": sigma_masses,
    "sigma mass error": sigma_mass_errors,
    "estimated sigma mass": sigma_mass_ests,
    "fpi": fpis,
    "fpi error": fpi_errors,
    "fpi2": fpi2s,
    "fpi2 error": fpi2_errors,
    "estimated fpi": fpi_ests,
    "mpi from fpi fit": mpi_from_fpis,
    "mpi from fpi error": mpi_from_fpi_errors,
    "mpi from fpi 2 fit": mpi_from_fpi2s,
    "mpi from fpi 2 error": mpi_from_fpi2_errors,
    "sigma vev": sigma_vevs,
    "sigma vev error": sigma_vev_errors,
    "sigma effective mass": sigma_mass_effs,
    "sigma effective mass error": sigma_mass_eff_errors,
    "# trajectories": trajs,
    "cutoff": cutoffs,
    "block size": block_sizes
})

df.to_csv("results/simulation_results.csv")
