# all C/C++ level import

from qlat.everything cimport *

cdef extern from "qlat-ext/gauge-fixing.h" namespace "qlat":

    cdef cppclass GaugeFixingInfo:
        Bool converged
        Int n_iter
        Int max_iter
        RealD eps
        RealD alpha
        RealD functional
        RealD divergence
        std_vector[RealD] functional_list
        std_vector[RealD] divergence_list

    RealD gf_gauge_functional(const GaugeField& gf, const std_vector[Int]& dirs) except +

    RealD gf_gauge_divergence(const GaugeField& gf, const std_vector[Int]& dirs) except +

    RealD gf_gauge_unitarity_error(const GaugeField& gf) except +

    GaugeFixingInfo gf_fix_gauge(GaugeField& gf, const std_vector[Int]& dirs, const Int max_iter, const RealD eps, const RealD alpha, const Bool verbose) except +

    GaugeFixingInfo gf_fix_landau_gauge(GaugeField& gf, const Int max_iter, const RealD eps, const RealD alpha, const Bool verbose) except +

    GaugeFixingInfo gf_fix_coulomb_gauge(GaugeField& gf, const Int max_iter, const RealD eps, const RealD alpha, const Bool verbose) except +
