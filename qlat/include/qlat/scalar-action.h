#pragma once

#include <qlat/setup.h>
#include <qlat/field-fft.h>

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
    // Returns the contribution to the total action from a single lattice
    // point (including the relavent neighbor interactions)
    // TIMER("ScalarAction.action_point");
    // Stores [sum i, mu] [phi_i(x+mu)-phi_i(x)]^2
    double dphi_sq=0;
    // Stores [sum i] phi_i(x)^2
    double phi_sq=0;
    // Stores phi_0(x)
    double phi_0=0;
    for (int m = 0; m < multiplicity; ++m) {
      double phi = sf.get_elem(xl,m);
      phi_sq += phi*phi;
      if (m==0) phi_0 = phi;
      for (int mu = 0; mu < 4; ++mu) {
        xl[mu]+=1;
        double phi_mu = sf.get_elem(xl,m);
        dphi_sq += (phi_mu-phi)*(phi_mu-phi);
        xl[mu]-=1;
      }
    }
    return dphi_sq/2.0 + m_sq*phi_sq/2.0 + lmbd*phi_sq*phi_sq/24.0 + alpha*phi_0;
  }
    
  inline double action_node_no_comm(const Field<double>& sf)
  {
	// Returns the total action of the portion of the lattice on the 
	// current node (assuming the necessary communication has already 
	// been done)
    TIMER("ScalarAction.action_node_no_comm");
    const Geometry geo = sf.geo();
    // Creates a geometry that is the same as the field geometry, except
    // with multiplicity 1
    const Geometry geo_r = geo_reform(geo);
    // Creates a field to save the contribution to the total action from
    // each point
    FieldM<double, 1> fd;
    fd.init(geo_r);
    // Loops over every lattice point in the current node
    qacc_for(index, geo_r.local_volume(), {
      fd.get_elem(index) = action_point(sf,geo.multiplicity,geo_r.coordinate_from_index(index));
    });
    // Sums over the contributions to the total action from each point
    // (this cannot be done in the previous loops because the previous
    // loop runs in parallel)
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
  
  qacc double hmc_mass_p(const Coordinate& L, const Coordinate& pg)
  {
    // Returns the momentum-dependent mass factor for HMC Fourier 
    // acceleration
    // TIMER("ScalarAction.hmc_mass_p");
    return 4/(PI*PI)*(m_sq + 8 - 2*(std::cos(2*PI*pg[0]/L[0]) +
									std::cos(2*PI*pg[1]/L[1]) + 
									std::cos(2*PI*pg[2]/L[2]) + 
									std::cos(2*PI*pg[3]/L[3])));
  }
  
  inline void hmc_estimate_mass(Field<double>& masses, const Field<Complex>& field_ft, const Field<Complex>& force_ft, const double phi0)
  {
    TIMER("ScalarAction.hmc_estimate_mass");
    const Geometry geo = field_ft.geo();
    //masses.init(geo);
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const long gindex = geo.g_index_from_g_coordinate(xg);
      Vector<double> masses_v = masses.get_elems(xl);
      int M = masses_v.size();
      qassert(M == geo.multiplicity);
      for (int m = 0; m < M; ++m) {
		Complex fld = field_ft.get_elem(xl,m);
		Complex frc = force_ft.get_elem(xl,m);
        if(gindex==0 && m==0){
          masses_v[m] = 4/(PI*PI)*std::pow((frc.real()*frc.real()+frc.imag()*frc.imag())/((fld.real()-phi0)*(fld.real()-phi0)+fld.imag()*fld.imag()), 0.5);
	    }
	    else {
          masses_v[m] = 4/(PI*PI)*std::pow((frc.real()*frc.real()+frc.imag()*frc.imag())/(fld.real()*fld.real()+fld.imag()*fld.imag()), 0.5);
		}
      }
    });
  }
  
  inline void to_mass_factor(Field<double>& sin_domega)
  {
    TIMER("ScalarAction.to_mass_factor");
    const Geometry geo = sin_domega.geo();
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      Vector<double> v = sin_domega.get_elems(xl);
      for (int m = 0; m < v.size(); ++m) {
        v[m] = 1 + 2*std::asin(v[m])/PI;
        v[m] = v[m]*v[m];
      }
    });
  }

  inline double hmc_m_hamilton_node(const Field<Complex>& sm_complex, const Field<double>& masses)
  {
    // Return the part of an HMC Hamiltonian due to the given momentum 
    // field (on the current node).
    TIMER("ScalarAction.hmc_m_hamilton_node");
    // Creates a complex copy of the real field so that we can compute 
    // the Fourier transform
    // static Field<Complex> sm_complex;
    // set_complex_from_double(sm_complex, sm);
    // Computes the Fourier transform of the field
    // fft_complex_field(sm_complex);
    // Saves the field geometry
    const Geometry geo = sm_complex.geo();
    // Creates a geometry that is the same as the field geometry, except
    // with multiplicity 1
    const Geometry geo_r = geo_reform(geo);
    // Creates a field to save the contribution to the sum 
    // from each point
    FieldM<double, 1> fd;
    fd.init(geo_r);
    // long V = geo.total_volume();
    qacc_for(index, geo_r.local_volume(), {
      Coordinate xl = geo_r.coordinate_from_index(index);
      double s=0;
      for (int m = 0; m < geo.multiplicity; ++m) {
        Complex c = sm_complex.get_elem(xl,m);
        s += (c.real()*c.real()+c.imag()*c.imag())/2/masses.get_elem(xl,m); // /hmc_mass_p(L,geo.coordinate_g_from_l(xl));
      }
      fd.get_elem(index) = s;
    });
    // Sums over the contributions to the sum from each point
    // (this cannot be done in the previous loops because the previous
    // loop runs in parallel)
    double sum = 0.0;
    for (long index = 0; index < geo_r.local_volume(); ++index) {
      sum += fd.get_elem(index);
    }
    return sum;
  }
  
  inline void hmc_set_force_no_comm(Field<double>& sm_force, const Field<double>& sf)
  {
    TIMER("ScalarAction.hmc_set_sm_force_no_comm");
    const Geometry geo = sf.geo();
    //sm_force.init(geo);
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
  
  inline void hmc_field_evolve(Field<Complex>& sf_complex, const Field<Complex>& sm_complex,
                               const Field<double>& masses, const double step_size)
  {
    TIMER("hmc_field_evolve");
    //Field<Complex> sf_complex;
    //Field<Complex> sm_complex;
    //set_complex_from_double(sf_complex, sf);
    //set_complex_from_double(sm_complex, sm);
    // Computes the Fourier transform of the fields
    const Geometry& geo = sf_complex.geo();
    // long rt_V = std::pow(geo.total_volume(),0.5);
    //fft_complex_field(sf_complex,true);
    //sf_complex*=1/rt_V;
    //fft_complex_field(sm_complex,true);
    //sm_complex*=1/rt_V;
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      Vector<Complex> sf_v = sf_complex.get_elems(xl);
      const Vector<Complex> sm_v = sm_complex.get_elems_const(xl);
      qassert(sf_v.size() == sm_v.size());
      for (int m = 0; m < sf_v.size(); ++m) {
        sf_v[m] = sf_v[m] + sm_v[m]*step_size/masses.get_elem(xl,m);
      }
    });
    //fft_complex_field(sf_complex, false);
    //sf_complex*=1/rt_V;
    //set_double_from_complex(sf, sf_complex);
  }
  
  inline void axial_current_node_no_comm(Field<double>&  axial_current, const Field<double>& sf)
  {
    TIMER("ScalarAction.axial_current_node_no_comm");
    const Geometry geo = sf.geo();
    // const Geometry geo_r = geo_reform(geo, geo.multiplicity-1);
    //axial_current.init(geo_r);
    qacc_for(index, geo.local_volume(), {
      Coordinate xl = geo.coordinate_from_index(index);
      Vector<double> ac_v = axial_current.get_elems(xl);
      int M = ac_v.size();
      qassert(M == geo.multiplicity-1);
      double p0;
      double pi;
      for (int m = 0; m < M; ++m) {
        xl[3]-=1;
        p0 = sf.get_elem(xl, 0);
        pi = sf.get_elem(xl, m+1);
        xl[3]+=1;
        ac_v[m] = p0*sf.get_elem(xl, m+1) - sf.get_elem(xl, 0)*pi;
      }
    });
  }
  
  inline double sum_sq(const Field<double>& f)
  {
    // Returns the sum of f(x)^2 over lattice sites (on the current 
    // node) and multiplicity
    TIMER("field_sum_sq");
    const Geometry geo = f.geo();
    // Creates a geometry that is the same as the field geometry, except
    // with multiplicity 1
    const Geometry geo_r = geo_reform(geo);
    // Creates a field to save the contribution to the sum of squares 
    // from each point
    FieldM<double, 1> fd;
    fd.init(geo_r);
    qacc_for(index, geo_r.local_volume(), {
      Coordinate xl = geo_r.coordinate_from_index(index);
      double s=0;
      for (int m = 0; m < geo.multiplicity; ++m) {
        double d = f.get_elem(xl,m);
        s += d*d;
      }
      fd.get_elem(index) = s;
    });
    // Sums over the contributions to the sum of squares from each point
    // (this cannot be done in the previous loops because the previous
     // loop runs in parallel)
    double sum = 0;
    for (long index = 0; index < geo_r.local_volume(); ++index) {
      sum += fd.get_elem(index);
    }
    return sum;
  }
  
  inline void axial_current_node(Field<double>& axial_current, const Field<double>& sf)
  {
	// Sets the axial_current field based on the provided field 
	// configuration sf. axial_current.get_elem(x,i) will give the time
	// component of the ith axial current vector at position x-a/2.
    TIMER("ScalarAction.axial_current_node");
    const Coordinate expand_left(0, 0, 0, 1);
    const Coordinate expand_right(0, 0, 0, 0);
    const Geometry geo_ext = geo_resize(sf.geo(), expand_left, expand_right);
    Field<double> sf_ext;
    sf_ext.init(geo_ext);
    sf_ext = sf;
    refresh_expanded(sf_ext);
    axial_current_node_no_comm(axial_current, sf_ext);
  }
  
  inline void hmc_set_rand_momentum(Field<Complex>& sm_complex, const Field<double>& masses, const RngState& rs)
  {
    TIMER("set_rand_momentum");
    // Note that momentum fields produced with this function need to be
    // projected to real momenta
    const Geometry& geo = sm_complex.geo();
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const long gindex = geo.g_index_from_g_coordinate(xg);
      RngState rsi = rs.newtype(gindex);
      Vector<Complex> v = sm_complex.get_elems(xl);
      for (int m = 0; m < v.size(); ++m) {
        double sigma = std::pow(masses.get_elem(xl,m), 0.5);
        v[m] = Complex(g_rand_gen(rsi, 0, sigma), g_rand_gen(rsi, 0, sigma));
      }
    });
  }
  
  inline void hmc_predict_field(Field<Complex>& field_ft, const Field<Complex>& momentum_ft, const Field<double>& masses, const double vev_sigma)
  {
    TIMER("hmc_predict_field");
    const Geometry& geo = momentum_ft.geo();
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const long gindex = geo.g_index_from_g_coordinate(xg);
      Vector<Complex> v = field_ft.get_elems(xl);
      int M = v.size();
      qassert(M == geo.multiplicity);
      qassert(M == masses.geo().multiplicity);
      for (int m = 0; m < M; ++m) {
		Complex mom = momentum_ft.get_elem(xl,m);
        if(gindex==0 && m==0) {
          v[m] = Complex(mom.real()/masses.get_elem(xl,m)/PI*2 + 
                         vev_sigma*std::pow(geo.total_volume(),0.5), 
                         mom.imag()/masses.get_elem(xl,m)/PI*2);
		}
        else {
          v[m] = Complex(mom.real()/masses.get_elem(xl,m)/PI*2,
                         mom.imag()/masses.get_elem(xl,m)/PI*2);
	    }
      }
    });
  }
  
  inline void get_polar_field(Field<double>& polar_fields, const Field<double>& sf)
  {
	// Sets the fields corresponding to polar coordinates based on the 
    // provided cartesian field configuration sf.
    TIMER("ScalarAction.get_polar_field");
    const Geometry geo = sf.geo();
    qassert(geo.multiplicity == 4);
    qacc_for(index, geo.local_volume(), {
      Coordinate xl = geo.coordinate_from_index(index);
      Vector<double> pf_v = polar_fields.get_elems(xl);
      qassert(pf_v.size()==4);
      double w = sf.get_elem(xl,0);
      double x = sf.get_elem(xl,1);
      double y = sf.get_elem(xl,2);
      double z = sf.get_elem(xl,3);
      // The formula to be inverted is
      // x_0 = r * cos(phi/r)
      // x_i = r * sin(phi/r) * phi_i/phi
      pf_v[0] = std::pow(w*w+x*x+y*y+z*z,0.5);
      double phi = pf_v[0]*std::acos(w/pf_v[0]);
      pf_v[1] = phi*x/(pf_v[0]*std::sin(phi/pf_v[0]));
      pf_v[2] = phi*y/(pf_v[0]*std::sin(phi/pf_v[0]));
      pf_v[3] = phi*z/(pf_v[0]*std::sin(phi/pf_v[0]));
    });
  }
   
};

}  // namespace qlat
