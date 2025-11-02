#pragma once

#include <qlat/field-fft.h>

namespace qlat
{  //

struct ScalarAction {
  bool initialized;
  RealD m_sq;
  RealD lmbd;
  RealD alpha;
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
  qacc ScalarAction(const RealD m_sq_, const RealD lmbd_, const RealD alpha_)
  {
    init();
    initialized = true;
    m_sq = m_sq_;
    lmbd = lmbd_;
    alpha = alpha_;
  }

  qacc RealD action_point(const Field<RealD>& sf, const Int multiplicity, Coordinate xl)
  {
    // Returns the contribution to the total action from a single lattice
    // point (including the relavent neighbor interactions)
    // TIMER("ScalarAction.action_point");
    // Stores [sum i, mu] [phi_i(x+mu)-phi_i(x)]^2
    RealD dphi_sq=0;
    // Stores [sum i] phi_i(x)^2
    RealD phi_sq=0;
    // Stores phi_0(x)
    RealD phi_0=0;
    for (Int m = 0; m < multiplicity; ++m) {
      RealD phi = sf.get_elem(xl,m);
      phi_sq += phi*phi;
      if (m==0) phi_0 = phi;
      for (Int mu = 0; mu < 4; ++mu) {
        xl[mu]+=1;
        RealD phi_mu = sf.get_elem(xl,m);
        dphi_sq += (phi_mu-phi)*(phi_mu-phi);
        xl[mu]-=1;
      }
    }
    return dphi_sq/2.0 + m_sq*phi_sq/2.0 + lmbd*phi_sq*phi_sq/24.0 + alpha*phi_0;
  }

  inline RealD action_node_no_comm(const Field<RealD>& sf)
  {
	// Returns the total action of the portion of the lattice on the
	// current node (assuming the necessary communication has already
	// been done)
    TIMER("ScalarAction.action_node_no_comm");
    const Geometry geo = sf.geo();
    // Creates a geometry that is the same as the field geometry, except
    // with multiplicity 1
    const Geometry geo_r = geo_resize(geo);
    // Creates a field to save the contribution to the total action from
    // each point
    FieldM<RealD, 1> fd;
    fd.init(geo_r);
    // Loops over every lattice point in the current node
    const Int multiplicity = sf.multiplicity;
    ScalarAction& sa = *this;
    qthread_for(index, geo_r.local_volume(), {
      const Geometry& geo = sf.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      fd.get_elem(index) = sa.action_point(sf, multiplicity, xl);
    });
    // Sums over the contributions to the total action from each point
    // (this cannot be done in the previous loops because the previous
    // loop runs in parallel)
    RealD sum = 0.0;
    for (Long index = 0; index < geo_r.local_volume(); ++index) {
      sum += fd.get_elem(index);
    }
    return sum;
  }

  inline RealD action_node(const Field<RealD>& sf)
  {
	// Return the total Euclidean action (on the current node) associated
	// with the given scalar field.
    TIMER("ScalarAction.action_node");
    // Get the values of the field sites that lie just outside of the
    // current node. This is needed to compute the derivative term.
    const Coordinate expand_left(0, 0, 0, 0);
    const Coordinate expand_right(1, 1, 1, 1);
    const Geometry geo_ext = geo_resize(sf.geo(), expand_left, expand_right);
    Field<RealD> sf_ext;
    sf_ext.init(geo_ext, sf.multiplicity);
    sf_ext = sf;
    refresh_expanded(sf_ext);
    // Calculate the energy of the scalar field
    return action_node_no_comm(sf_ext);
  }

  qacc RealD hmc_mass_p(const Coordinate& L, const Coordinate& pg)
  {
    // Returns the momentum-dependent mass factor for HMC Fourier
    // acceleration
    // TIMER("ScalarAction.hmc_mass_p");
    return 4/(PI*PI)*(m_sq + 8 - 2*(std::cos(2*PI*pg[0]/L[0]) +
									std::cos(2*PI*pg[1]/L[1]) +
									std::cos(2*PI*pg[2]/L[2]) +
									std::cos(2*PI*pg[3]/L[3])));
  }

  inline void hmc_estimate_mass(Field<RealD>& masses, const Field<ComplexD>& field_ft, const Field<ComplexD>& force_ft, const RealD phi0)
  {
    TIMER("ScalarAction.hmc_estimate_mass");
    const Geometry geo = field_ft.geo();
    const Int multiplicity = field_ft.multiplicity;
    (void) multiplicity;
    //masses.init(geo);
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = field_ft.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const Long gindex = geo.g_index_from_g_coordinate(xg);
      Vector<RealD> masses_v = masses.get_elems(xl);
      Int M = masses_v.size();
      qassert(M == multiplicity);
      for (Int m = 0; m < M; ++m) {
        ComplexD fld = field_ft.get_elem(xl, m);
        ComplexD frc = force_ft.get_elem(xl, m);
        if(gindex==0 && m==0){
          masses_v[m] = 4/(PI*PI)*std::pow((frc.real()*frc.real()+frc.imag()*frc.imag())/((fld.real()-phi0)*(fld.real()-phi0)+fld.imag()*fld.imag()), 0.5);
        } else {
          masses_v[m] = 4/(PI*PI)*std::pow((frc.real()*frc.real()+frc.imag()*frc.imag())/(fld.real()*fld.real()+fld.imag()*fld.imag()), 0.5);
        }
      }
    });
  }

  inline void to_mass_factor(Field<RealD>& sin_domega)
  {
    TIMER("ScalarAction.to_mass_factor");
    const Geometry geo = sin_domega.geo();
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = sin_domega.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      Vector<RealD> v = sin_domega.get_elems(xl);
      for (Int m = 0; m < v.size(); ++m) {
        v[m] = 1 + 2*std::asin(v[m])/PI;
        v[m] = v[m]*v[m];
      }
    });
  }

  inline RealD hmc_m_hamilton_node(const Field<ComplexD>& sm_complex, const Field<RealD>& masses)
  {
    // Return the part of an HMC Hamiltonian due to the given momentum
    // field (on the current node).
    TIMER("ScalarAction.hmc_m_hamilton_node");
    // Creates a complex copy of the real field so that we can compute
    // the Fourier transform
    // static Field<ComplexD> sm_complex;
    // set_complex_from_double(sm_complex, sm);
    // Computes the Fourier transform of the field
    // fft_complex_field(sm_complex);
    // Saves the field geometry
    const Geometry geo = sm_complex.geo();
    const Int multiplicity = sm_complex.multiplicity;
    // Creates a geometry that is the same as the field geometry, except
    // with multiplicity 1
    const Geometry geo_r = geo_resize(geo);
    // Creates a field to save the contribution to the sum
    // from each point
    FieldM<RealD, 1> fd;
    fd.init(geo_r);
    // Long V = geo.total_volume();
    qthread_for(index, geo_r.local_volume(), {
      const Geometry& geo = sm_complex.geo();
      Coordinate xl = geo.coordinate_from_index(index);
      RealD s=0;
      for (Int m = 0; m < multiplicity; ++m) {
        ComplexD c = sm_complex.get_elem(xl,m);
        s += (c.real()*c.real()+c.imag()*c.imag())/2/masses.get_elem(xl,m); // /hmc_mass_p(L,geo.coordinate_g_from_l(xl));
      }
      fd.get_elem(index) = s;
    });
    // Sums over the contributions to the sum from each point
    // (this cannot be done in the previous loops because the previous
    // loop runs in parallel)
    RealD sum = 0.0;
    for (Long index = 0; index < geo.local_volume(); ++index) {
      sum += fd.get_elem(index);
    }
    return sum;
  }

  inline void hmc_set_force_no_comm(Field<RealD>& sm_force, const Field<RealD>& sf)
  {
    TIMER("ScalarAction.hmc_set_sm_force_no_comm");
    const Geometry geo = sf.geo();
    const Int multiplicity = sf.multiplicity;
    (void) multiplicity;
    //sm_force.init(geo);
    ScalarAction& sa = *this;
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = sf.geo();
      Coordinate xl = geo.coordinate_from_index(index);
      Vector<RealD> sm_force_v = sm_force.get_elems(xl);
      Int M = sm_force_v.size();
      qassert(M == multiplicity);
      RealD sum_mult_sq = 0.0;
      for (Int m = 0; m < M; ++m) {
        sum_mult_sq += pow(sf.get_elem(xl, m), 2);
      }
      for (Int m = 0; m < M; ++m) {
        sm_force_v[m] =
            (2 * 4 + sa.m_sq + sa.lmbd / 6 * sum_mult_sq) * sf.get_elem(xl, m);
        if (m == 0) sm_force_v[m] += sa.alpha;
        for (Int dir = 0; dir < 4; ++dir) {
          xl[dir] += 1;
          sm_force_v[m] -= sf.get_elem(xl, m);
          xl[dir] -= 2;
          sm_force_v[m] -= sf.get_elem(xl, m);
          xl[dir] += 1;
        }
      }
    });
  }

  inline void hmc_set_force(Field<RealD>&  sm_force, const Field<RealD>&  sf)
  {
	// Calculate and set the HMC force field based on the given field
	// configuration.
    TIMER("ScalarAction.hmc_set_sm_force");
    Coordinate expand_left(1, 1, 1, 1);
    Coordinate expand_right(1, 1, 1, 1);
    const Geometry geo_ext = geo_resize(sf.geo(), expand_left, expand_right);
    Field<RealD> sf_ext;
    sf_ext.init(geo_ext, sf.multiplicity);
    sf_ext = sf;
    refresh_expanded(sf_ext);
    hmc_set_force_no_comm(sm_force, sf_ext);
  }

  inline void hmc_field_evolve(Field<ComplexD>& sf_complex, const Field<ComplexD>& sm_complex,
                               const Field<RealD>& masses, const RealD step_size)
  {
    TIMER("hmc_field_evolve");
    //Field<ComplexD> sf_complex;
    //Field<ComplexD> sm_complex;
    //set_complex_from_double(sf_complex, sf);
    //set_complex_from_double(sm_complex, sm);
    // Computes the Fourier transform of the fields
    const Geometry& geo = sf_complex.geo();
    // Long rt_V = std::pow(geo.total_volume(),0.5);
    //fft_complex_field(sf_complex,true);
    //sf_complex*=1/rt_V;
    //fft_complex_field(sm_complex,true);
    //sm_complex*=1/rt_V;
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = sf_complex.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      Vector<ComplexD> sf_v = sf_complex.get_elems(xl);
      const Vector<ComplexD> sm_v = sm_complex.get_elems_const(xl);
      qassert(sf_v.size() == sm_v.size());
      for (Int m = 0; m < sf_v.size(); ++m) {
        sf_v[m] = sf_v[m] + sm_v[m]*step_size/masses.get_elem(xl,m);
      }
    });
    //fft_complex_field(sf_complex, false);
    //sf_complex*=1/rt_V;
    //set_double_from_complex(sf, sf_complex);
  }

  inline void axial_current_node_no_comm(Field<RealD>&  axial_current, const Field<RealD>& sf)
  {
    TIMER("ScalarAction.axial_current_node_no_comm");
    const Geometry geo = sf.geo();
    const Int multiplicity = sf.multiplicity;
    (void) multiplicity;
    // const Geometry geo_r = geo_resize(geo);
    // axial_current.init(geo_r, sf.multiplicity-1);
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = sf.geo();
      Coordinate xl = geo.coordinate_from_index(index);
      Vector<RealD> ac_v = axial_current.get_elems(xl);
      Int M = ac_v.size();
      qassert(M == multiplicity-1);
      RealD p0;
      RealD pi;
      for (Int m = 0; m < M; ++m) {
        xl[3]-=1;
        p0 = sf.get_elem(xl, 0);
        pi = sf.get_elem(xl, m+1);
        xl[3]+=1;
        ac_v[m] = p0*sf.get_elem(xl, m+1) - sf.get_elem(xl, 0)*pi;
      }
    });
  }

  inline RealD sum_sq(const Field<RealD>& f)
  {
    // Returns the sum of f(x)^2 over lattice sites (on the current
    // node) and multiplicity
    TIMER("field_sum_sq");
    const Geometry geo = f.geo();
    const Int multiplicity = f.multiplicity;
    // Creates a geometry that is the same as the field geometry, except
    // with multiplicity 1
    const Geometry geo_r = geo_resize(geo);
    // Creates a field to save the contribution to the sum of squares
    // from each point
    FieldM<RealD, 1> fd;
    fd.init(geo_r);
    qthread_for(index, geo_r.local_volume(), {
      // const Geometry& geo = f.geo();
      Coordinate xl = geo_r.coordinate_from_index(index);
      RealD s=0;
      for (Int m = 0; m < multiplicity; ++m) {
        RealD d = f.get_elem(xl,m);
        s += d*d;
      }
      fd.get_elem(index) = s;
    });
    // Sums over the contributions to the sum of squares from each point
    // (this cannot be done in the previous loops because the previous
     // loop runs in parallel)
    RealD sum = 0;
    for (Long index = 0; index < geo_r.local_volume(); ++index) {
      sum += fd.get_elem(index);
    }
    return sum;
  }

  inline void axial_current_node(Field<RealD>& axial_current, const Field<RealD>& sf)
  {
	// Sets the axial_current field based on the provided field
	// configuration sf. axial_current.get_elem(x,i) will give the time
	// component of the ith axial current vector at position x-a/2.
    TIMER("ScalarAction.axial_current_node");
    const Coordinate expand_left(0, 0, 0, 1);
    const Coordinate expand_right(0, 0, 0, 0);
    const Geometry geo_ext = geo_resize(sf.geo(), expand_left, expand_right);
    Field<RealD> sf_ext;
    sf_ext.init(geo_ext, sf.multiplicity);
    sf_ext = sf;
    refresh_expanded(sf_ext);
    axial_current_node_no_comm(axial_current, sf_ext);
  }

  inline void hmc_set_rand_momentum(Field<ComplexD>& sm_complex, const Field<RealD>& masses, const RngState& rs)
  {
    TIMER("set_rand_momentum");
    // Note that momentum fields produced with this function need to be
    // projected to real momenta
    const Geometry& geo = sm_complex.geo();
    qthread_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const Long gindex = geo.g_index_from_g_coordinate(xg);
      RngState rsi = rs.newtype(gindex);
      Vector<ComplexD> v = sm_complex.get_elems(xl);
      for (Int m = 0; m < v.size(); ++m) {
        RealD sigma = std::pow(masses.get_elem(xl,m), 0.5);
        v[m] = ComplexD(g_rand_gen(rsi, 0, sigma), g_rand_gen(rsi, 0, sigma));
      }
    });
  }

  inline void hmc_predict_field(Field<ComplexD>& field_ft, const Field<ComplexD>& momentum_ft, const Field<RealD>& masses, const RealD vev_sigma)
  {
    TIMER("hmc_predict_field");
    const Geometry& geo = momentum_ft.geo();
    const Int multiplicity = momentum_ft.multiplicity;
    (void) multiplicity;
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = momentum_ft.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const Long gindex = geo.g_index_from_g_coordinate(xg);
      Vector<ComplexD> v = field_ft.get_elems(xl);
      Int M = v.size();
      qassert(M == multiplicity);
      qassert(M == masses.multiplicity);
      for (Int m = 0; m < M; ++m) {
        ComplexD mom = momentum_ft.get_elem(xl, m);
        if (gindex == 0 && m == 0) {
          v[m] = ComplexD(mom.real() / masses.get_elem(xl, m) / PI * 2 +
                             vev_sigma * std::pow(geo.total_volume(), 0.5),
                         mom.imag() / masses.get_elem(xl, m) / PI * 2);
        } else {
          v[m] = ComplexD(mom.real() / masses.get_elem(xl, m) / PI * 2,
                         mom.imag() / masses.get_elem(xl, m) / PI * 2);
        }
      }
    });
  }

  inline void get_polar_field(Field<RealD>& polar_fields, const Field<RealD>& sf)
  {
	// Sets the fields corresponding to polar coordinates based on the
    // provided cartesian field configuration sf.
    TIMER("ScalarAction.get_polar_field");
    const Geometry geo = sf.geo();
    const Int multiplicity = sf.multiplicity;
    Qassert(multiplicity == 4);
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = sf.geo();
      Coordinate xl = geo.coordinate_from_index(index);
      //const Vector<RealD> sf_v = sf.get_elems_const(xl);
      qassert(sf.get_elems_const(xl).size() == 4);
      Vector<RealD> pf_v = polar_fields.get_elems(xl);
      qassert(pf_v.size() == 4);
      RealD w = sf.get_elem(xl,0);
      RealD x = sf.get_elem(xl,1);
      RealD y = sf.get_elem(xl,2);
      RealD z = sf.get_elem(xl,3);
      // The formula to be inverted is
      // x_0 = r * cos(phi/r)
      // x_i = r * sin(phi/r) * phi_i/phi
      pf_v[0] = std::pow(w*w+x*x+y*y+z*z,0.5);
      RealD phi = pf_v[0]*std::acos(w/pf_v[0]);
      pf_v[1] = phi*x/(pf_v[0]*std::sin(phi/pf_v[0]));
      pf_v[2] = phi*y/(pf_v[0]*std::sin(phi/pf_v[0]));
      pf_v[3] = phi*z/(pf_v[0]*std::sin(phi/pf_v[0]));
    });
  }

};

}  // namespace qlat
