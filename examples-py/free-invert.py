#!/usr/bin/env python3

import qlat as q

q.begin_with_mpi()

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
q.json_results_append(f"free-invert: geo.show()={geo.show()}")
rs = q.RngState("seed")

gf = q.GaugeField(geo)
gf.set_unit()
gf.show_info()

fa = q.FermionAction(mass=0.05, ls=16, m5=1.0)

qinv_free = q.InverterDwfFreeField(
    mass=fa.mass(), m5=fa.m5(), qtimer=q.Timer("py:InverterDwfFreeField")
)

qinv_dwf = q.InverterDomainWall(gf=gf, fa=fa, qtimer=q.Timer("py:InverterDomainWall"))

src_p = q.mk_point_src(geo, q.Coordinate([0, 0, 0, 0]))

src_r = q.Prop(geo)
src_r.set_rand(rs.split("src_r"))

for src in [src_p, src_r]:
    sol = qinv_free * src

    sol1 = qinv_dwf * src

    sol_diff = sol1.copy()
    sol_diff -= sol

    q.json_results_append("free-invert: sol qnorm", sol.qnorm(), 1e-10)
    q.json_results_append("free-invert: sol1 qnorm", sol1.qnorm(), 1e-10)
    q.json_results_append("free-invert: sol_diff qnorm", sol_diff.qnorm(), 1e-10)
    assert sol_diff.qnorm() < 1e-7

    sol2 = q.free_invert(
        src, mass=fa.mass(), m5=fa.m5(), momtwist=q.CoordinateD([0, 0, 0, 0])
    )
    sol2_diff = sol2.copy()
    sol2_diff -= sol
    assert sol2_diff.qnorm() < 1e-9

fa_z = q.FermionAction(
    mass=0.05,
    ls=12,
    m5=1.0,
    omega=[
        1.0903256131299373,
        0.9570283702230611,
        0.7048886040934104,
        0.48979921782791747,
        0.328608311201356,
        0.21664245377015995,
        0.14121112711957107,
        0.0907785101745156,
        complex(0.05608303440064219, -0.007537158177840385),
        complex(0.05608303440064219, 0.007537158177840385),
        complex(0.0365221637144842, -0.03343945161367745),
        complex(0.0365221637144842, 0.03343945161367745),
    ],
)

qinv_dwf_z = q.InverterDomainWall(
    gf=gf, fa=fa_z, qtimer=q.Timer("py:InverterDomainWall(ZMobius)")
)

for src in [src_p, src_r]:
    sol_free = q.free_invert(
        src, mass=fa_z.mass(), m5=fa_z.m5(), momtwist=q.CoordinateD([0, 0, 0, 0])
    )

    sol_z = qinv_dwf_z * src

    sol_z_diff = sol_z.copy()
    sol_z_diff -= sol_free

    q.json_results_append("free-invert: sol_free qnorm", sol_free.qnorm(), 1e-10)
    q.json_results_append("free-invert: sol_z qnorm", sol_z.qnorm(), 1e-10)
    q.json_results_append("free-invert: sol_z_diff qnorm", sol_z_diff.qnorm(), 1e-8)

q.check_log_json(__file__)
q.timer_display()
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
