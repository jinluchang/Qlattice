import gpt as g
import numpy as np
import sys, os

g.message("Topological charge measurement with GPT")
g.message("by Christoph Lehner")
g.message("2024/05/17")

config = g.default.get_single("--config",None)

def smear(Uin, c1, epsilon):
    Usm = g.copy(Uin)
    if c1 == 0.0:
        action = g.qcd.gauge.action.wilson(2.0 * Uin[0].otype.shape[0])
    else:
        action = g.qcd.gauge.action.improved_with_rectangle(2.0 * Uin[0].otype.shape[0], c1)
    g.algorithms.integrator.euler(
        Usm, lambda: [g(-u) for u in action.gradient(Usm, Usm)], 1.0
    )(epsilon)
    return Usm

smear_info_list = [
    [ 20, 0.05, 0.0 ],
    [ 20, 0.05, 0.0 ],
    [ 20, 0.05, 0.0 ],
    [ 50, 0.01, -1.4008 ],
    [ 50, 0.01, -1.4008 ],
    [ 50, 0.01, -1.4008 ],
    [ 50, 0.01, -1.4008 ]
]

U = g.load(config)

vol3d = float(np.prod(U[0].grid.gdimensions[0:3]))

w = g.corr_io.writer(f"{config}.gluonic_luchang")
output = g.gpt_io.writer(f"{config}.gluonic_fields_luchang", mpi=[2,2,2,4])

# first plaquette
g.message("Plaquette")
P = g.slice(g.qcd.gauge.rectangle(U, 1, 1, field=True) / vol3d, 3)
w.write("P", P)

fQ = open(f"{config}.Q_luchang","wt")
fE = open(f"{config}.E_luchang","wt")

tau = 0.0
U_wf = U
c = {}
for nsmear, epsilon, c1 in smear_info_list:

    for i in range(nsmear):
        tau += epsilon
        U_wf = smear(U_wf, c1, epsilon)
        g.message("%g" % tau, g.qcd.gauge.plaquette(U_wf))

    g.message("Field Strength")
    E_field = g.qcd.gauge.energy_density(U_wf, field=True)
    E = g.slice(E_field / vol3d, 3)
    w.write("E(%g)" % tau, E)

    E = sum(E).real / len(E)
    t2E = tau**2 * E
    g.message("t2E = ", t2E)

    g.message("Topology")
    Q_field = g.qcd.gauge.topological_charge_5LI(U_wf, cache=c, field=True)
    Q = g.slice(Q_field / vol3d, 3)
    w.write("Q(%g)" % tau, Q)

    Q = sum(Q).real / len(Q)
    fQ.write("%g %.15g\n" % (tau, Q))
    fQ.flush()

    output.write({
        "t" : tau,
        "Q" : Q_field,
        "E" : E_field
    })
    output.flush()

    fE.write("%g %.15g\n" % (tau, E))
    fE.flush()

exit()
