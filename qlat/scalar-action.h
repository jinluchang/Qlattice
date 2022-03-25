#pragma once

#include <qlat/config.h>

namespace qlat
{  //

struct ScalarAction {
  bool initialized;
  double m_sq;
  double lmbd;
  double alpha;
  //
  qacc void init()
  {
    initialized = false;
    m_sq = 1.0;
    lmbd = 1.0;
    alpha = 0.0;
  }
  //
  qacc ScalarAction() { init(); }
  qacc ScalarAction(const double m_sq_, const double lmbd_, const double alpha_)
  {
    init();
    initialized = true;
    m_sq = m_sq_;
    lmbd = lmbd_;
    alpha = alpha_;
  }
  
  qacc double action_point(const Field<double>& sf, const int multiplicity, Coordinate xl)
  {
    // TIMER("ScalarAction.action_point");
    double dphi_sq=0;
    double phi_sq=0;
    double phi_0=0;
    for (int m = 0; m < multiplicity; ++m) {
      double phi = sf.get_elem(xl,m);
      phi_sq += phi*phi;
      if (m==0) phi_0 += phi;
      for (int dir = 0; dir < 4; ++dir) {
        xl[dir]+=1;
        double phi_n = sf.get_elem(xl,m);
        dphi_sq += (phi_n-phi)*(phi_n-phi);
        xl[dir]-=1;
      }
    }
    return dphi_sq/2.0 + m_sq*phi_sq/2.0 + lmbd*phi_sq*phi_sq/24.0 + alpha*phi_0;
  }
    
  inline double action_node_no_comm(const Field<double>& sf)
  {
    TIMER("ScalarAction.action_node_no_comm");
    const Geometry geo = sf.geo();
    const Geometry geo_r = geo_reform(geo);
    FieldM<double, 1> fd;
    fd.init(geo_r);
    qacc_for(index, geo_r.local_volume(), {
      fd.get_elem(index) = action_point(sf,geo.multiplicity,geo_r.coordinate_from_index(index));
    });
    double sum = 0.0;
    for (long index = 0; index < geo_r.local_volume(); ++index) {
      sum += fd.get_elem(index);
    }
    return sum;
  }
  
  inline double action_node(const Field<double>& sf)
  {
	// Return the total Euclidean action (on the current node) associated 
	// with the given scalar field.
    TIMER("ScalarAction.action_node");
    // Get the values of the field sites that lie just outside of the 
    // current node. This is needed to compute the derivative term.
    const Coordinate expand_left(0, 0, 0, 0);
    const Coordinate expand_right(1, 1, 1, 1);
    const Geometry geo_ext = geo_resize(sf.geo(), expand_left, expand_right);
    Field<double> sf_ext;
    sf_ext.init(geo_ext);
    sf_ext = sf;
    refresh_expanded(sf_ext);
    // Calculate the energy of the scalar field
    return action_node_no_comm(sf_ext);
  }

  inline double sum_sq(const Field<double>& sf)
  {
	// Return the sum over lattice sites (x) and multiplicity (m) of 
	// phi_m(x)^2
    TIMER("ScalarAction.sum_sq");
    const Geometry geo = sf.geo();
    const Geometry geo_r = geo_reform(geo);
    FieldM<double, 1> fd;
    fd.init(geo_r);
    qacc_for(index, geo_r.local_volume(), {
      Coordinate xl = geo_r.coordinate_from_index(index);
      double s=0;
      for (int m = 0; m < geo.multiplicity; ++m) {
        double d = sf.get_elem(xl,m);
        s += d*d;
      }
      fd.get_elem(index) = s;
    });
    double sum = 0.0;
    for (long index = 0; index < geo_r.local_volume(); ++index) {
      sum += fd.get_elem(index);
    }
    return sum;
  }

  inline double hmc_m_hamilton_node(const Field<double>& sm)
  {
    // Return the part of an HMC Hamiltonian due to the given momentum 
    // field (on the current node).
    TIMER("ScalarAction.hmc_m_hamilton_node");
    return sum_sq(sm)/2.0;
  }
  
  inline void hmc_set_force_no_comm(Field<double>& sm_force, const Field<double>& sf)
  {
    TIMER("ScalarAction.hmc_set_sm_force_no_comm");
    const Geometry geo = sf.geo();
    sm_force.init(geo);
    qacc_for(index, geo.local_volume(), {
      Coordinate xl = geo.coordinate_from_index(index);
      Vector<double> sm_force_v = sm_force.get_elems(xl);
      int M = sm_force_v.size();
      qassert(M == geo.multiplicity);
      double sum_mult_sq = 0.0;
      for (int m = 0; m < M; ++m) {
        sum_mult_sq += pow(sf.get_elem(xl, m), 2);
      }
      for (int m = 0; m < M; ++m) {
        sm_force_v[m] =
            (2 * 4 + m_sq + lmbd / 6 * sum_mult_sq) * sf.get_elem(xl, m);
        if (m == 0) sm_force_v[m] += alpha;
        for (int dir = 0; dir < 4; ++dir) {
          xl[dir] += 1;
          sm_force_v[m] -= sf.get_elem(xl, m);
          xl[dir] -= 2;
          sm_force_v[m] -= sf.get_elem(xl, m);
          xl[dir] += 1;
        }
      }
    });
  }
  
  inline void hmc_set_force(Field<double>&  sm_force, const Field<double>&  sf)
  {
	// Calculate and set the HMC force field based on the given field 
	// configuration.
    TIMER("ScalarAction.hmc_set_sm_force");
    Coordinate expand_left(1, 1, 1, 1);
    Coordinate expand_right(1, 1, 1, 1);
    const Geometry geo_ext = geo_resize(sf.geo(), expand_left, expand_right);
    Field<double> sf_ext;
    sf_ext.init(geo_ext);
    sf_ext = sf;
    refresh_expanded(sf_ext);
    hmc_set_force_no_comm(sm_force, sf_ext);
  }
  
  inline void hmc_sf_evolve(Field<double>& sf, const Field<double>& sm,
                            const double step_size)
  {
    TIMER("hmc_sf_evolve");
    const Geometry& geo = sf.geo();
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      Vector<double> sf_v = sf.get_elems(xl);
      const Vector<double> sm_v = sm.get_elems_const(xl);
      qassert(sf_v.size() == sm_v.size());
      for (int m = 0; m < sf_v.size(); ++m) {
        sf_v[m] = sf_v[m] + sm_v[m]*step_size;
      }
    });
  }
  
  inline void axial_current_node_no_comm(Field<double>&  axial_current, const Field<double>& sf)
  {
    TIMER("ScalarAction.axial_current_node_no_comm");
    const Geometry geo = sf.geo();
    const Geometry geo_r = geo_reform(geo, geo.multiplicity-1);
    axial_current.init(geo_r);
    qacc_for(index, geo.local_volume(), {
      Coordinate xl = geo.coordinate_from_index(index);
      Vector<double> ac_v = axial_current.get_elems(xl);
      int M = ac_v.size();
      qassert(M == geo.multiplicity-1);
      double p0;
      double pi;
      for (int m = 0; m < M; ++m) {
        p0 = sf.get_elem(xl, 0);
        pi = sf.get_elem(xl, m+1);
        xl[3]+=1;
        ac_v[m] = p0*sf.get_elem(xl, m+1) - sf.get_elem(xl, 0)*pi;
        xl[3]-=1;
      }
    });
  }
  
  inline void axial_current_node(Field<double>&  axial_current, const Field<double>& sf)
  {
	// Sets the axial_current field based on the provided field 
	// configuration sf. axial_current.get_elem(x,i) will give the time
	// component of the ith axial current vector at position x+a/2.
    TIMER("ScalarAction.axial_current_node");
    const Coordinate expand_left(0, 0, 0, 0);
    const Coordinate expand_right(1, 1, 1, 1);
    const Geometry geo_ext = geo_resize(sf.geo(), expand_left, expand_right);
    Field<double> sf_ext;
    sf_ext.init(geo_ext);
    sf_ext = sf;
    refresh_expanded(sf_ext);
    axial_current_node_no_comm(axial_current, sf_ext);
  }
  
};

}  // namespace qlat
