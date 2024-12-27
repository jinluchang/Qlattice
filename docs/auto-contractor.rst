auto_contractor
===============

.. automodule:: auto_contractor
   :members:

Evaluation
----------

.. autosummary::
   :toctree: generated

   contract_simplify_compile
   filter_diagram_type
   contract_simplify
   compile_expr
   display_cexpr
   cexpr_code_gen_py

   CExpr

   cache_compiled_cexpr

   CCExpr

   eval_cexpr
   get_expr_names
   get_diagram_type_dict

   mk_sym
   mk_fac

   aff.rel_mod
   aff.rel_mod_sym
   aff.rel_mod_sym
   aff.c_rel_mod_sqr

Operators
---------

.. autosummary::
   :toctree: generated

   mk_scalar
   mk_scalar5
   mk_vec_mu
   mk_vec5_mu
   mk_meson
   mk_pi_0
   mk_pi_p
   mk_pi_m
   mk_a0_0
   mk_a0_p
   mk_a0_m
   mk_k_p
   mk_k_m
   mk_k_0
   mk_k_0_bar
   mk_sigma
   mk_kappa_p
   mk_kappa_m
   mk_kappa_0
   mk_kappa_0_bar
   mk_k_0_star_mu
   mk_k_0_star_bar_mu

   mk_pipi_i22
   mk_pipi_i21
   mk_pipi_i11
   mk_pipi_i20
   mk_pipi_i10
   mk_pipi_i0
   mk_kk_i11
   mk_kk_i10
   mk_kk_i0
   mk_k0k0bar
   mk_k0pi0
   mk_k0barpi0
   mk_kppim
   mk_kmpip
   mk_kpi_0_i1half
   mk_kpi_p_i1half
   mk_kpi_m_i3halves
   mk_kpi_0_i3halves
   mk_kpi_p1_i3halves
   mk_kpi_p2_i3halves

   mk_m

   mk_j5pi_mu
   mk_j5k_mu
   mk_j5km_mu
   mk_jpi_mu
   mk_jk_mu

   mk_j_mu
   mk_jl_mu
   mk_j0_mu
   mk_j10_mu
   mk_j11_mu
   mk_j1n1_mu

   mk_4qOp_VV
   mk_4qOp_VA
   mk_4qOp_AV
   mk_4qOp_AA
   mk_4qOp_SS
   mk_4qOp_SP
   mk_4qOp_PS
   mk_4qOp_PP
   mk_4qOp_LL
   mk_4qOp_LR
   mk_4qOp_RL
   mk_4qOp_RR
   mk_4qOp_LL_cmix
   mk_4qOp_LR_cmix
   mk_Qsub
   mk_Q1
   mk_Q2
   mk_Q3
   mk_Q4
   mk_Q5
   mk_Q6
   mk_Q7
   mk_Q8
   mk_Q9
   mk_Q10
   mk_Q0_b81
   mk_Q1_b81
   mk_Q2_b81
   mk_Q3_b81
   mk_Q4_b81
   mk_Q5_b81
   mk_Q6_b81
   mk_Q7_b81
   mk_Q8_b81

3-flavor operators in (8,1) representation
mk_Q1_b81
mk_Q2_b81
mk_Q3_b81
mk_Q4_b81

subtraction operators
mk_Q0_b81 ( = mk_Qsub )

charm-contained operators in (8,1) representation
mk_Q5_b81
mk_Q6_b81
mk_Q7_b81
mk_Q8_b81

:math:`Q_a^{e/o} = A_a^{e/o} Q_0^{e/o} + M_{a,i} Q_i^{e/o} ( i = 1, ... ,4; a = 5, ... ,8 )`

Tutorials
---------

Single propagator: ``examples-py/auto-contract-01.py``

.. literalinclude:: ../examples-py/auto-contract-01.py

Two propagators and a operator: ``examples-py/auto-contract-02.py``

.. literalinclude:: ../examples-py/auto-contract-02.py

Pion correlator and a operator: ``examples-py/auto-contract-03.py``

.. literalinclude:: ../examples-py/auto-contract-03.py

Examples
--------

Some examples: ``examples-py/auto-cexprs.py``

.. literalinclude:: ../examples-py/auto-cexprs.py

More examples: ``examples-py/auto-cexprs-more.py``

.. literalinclude:: ../examples-py/auto-cexprs-more.py

K->pipi example: ``examples-py/auto-cexprs-kpipi.py``

.. literalinclude:: ../examples-py/auto-cexprs-kpipi.py

A more complete example: ``examples-py-gpt/gpt-qlat-auto-simple.py``

.. literalinclude:: ../examples-py-gpt/gpt-qlat-auto-simple.py
