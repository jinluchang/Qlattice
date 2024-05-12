#pragma once

namespace qlat
{  //

struct QMAction {
  bool initialized;
  double alpha;
  double beta;
  double center_TV;
  double center_bar;
  double barrier_strength;
  double M;
  double L;
  Long t_full1;
  Long t_full2;
  Long t_FV;
  double dt;
  //
  qacc void init()
  {
    initialized = false;
    alpha = 1.0;
    beta = 1.0;
    center_TV = 2.0;
    center_bar = 1.0;
    barrier_strength = 1.0;
    M = 1.0;
    L = 0.0;
    t_full1 = 10;
    t_full2 = 10;
    t_FV = 10;
    dt = 1.0;
  }
  //
  qacc QMAction() { init(); }
  qacc QMAction(const double alpha_, const double beta_,
                const double barrier_strength_, const double M_,
                const double L_, const Long t_full1_, const Long t_full2_, 
                const Long t_FV_, const double dt_)
  {
    init();
    initialized = true;
    alpha = alpha_;
    beta = beta_;
    center_TV = (2.0-2.0*std::pow(1-alpha, 0.5))/alpha; // (3.0+std::pow(9.0-8.0*alpha, 0.5))/2.0/alpha;
    center_bar = (3.0-std::pow(9.0-8.0*alpha, 0.5))/2.0/alpha;
    barrier_strength = barrier_strength_;
    M = M_;
    L = L_;
    t_full1 = t_full1_;
    t_full2 = t_full2_;
    t_FV = t_FV_;
    dt = dt_;
  }

  inline double V(const double x, const Long t)
  {
    // Returns the potential evaluated at point x
    if(t<t_full1)
        return V_full(x);
    else if(t<t_full1+t_FV)
      return V_FV(x);
    else if(t<t_full1+t_FV+t_full2)
      return V_full(x);
    else {
     //displayln(ssprintf("t (in V): %ld", t));
     //displayln(ssprintf("x (in V): %24.17E", x));
     //displayln(ssprintf("V_TV(x) (in V): %24.17E", V_TV(x)));
     return V_TV(x);
    }
  }

  inline double dV(const double x, const Long t)
  {
    // Returns the potential evaluated at point x
    if(t<t_full1)
        return dV_full(x);
    else if(t<t_full1+t_FV)
      return dV_FV(x);
    else if(t<t_full1+t_FV+t_full2)
      return dV_full(x);
    else
     return dV_TV(x);
  }

  inline double V_phi4(const double x)
  {
    return beta*(x*x/2.0-x*x*x/2.0+alpha*x*x*x*x/8.0);
  }

  inline double dV_phi4(const double x)
  {
    // Returns the derivative of the potential with respect to x
    return beta*(x-x*x*3.0/2.0+alpha*x*x*x/2.0);
  }

  inline double V_full(const double x)
  {
    // Returns the potential evaluated at point x
    if(x>center_TV) {
      return 0;
    }
    return V_phi4(x); // - V_phi4(center_TV);
  }

  inline double dV_full(const double x)
  {
    double rtn = dV_phi4(x);
    if(x>center_TV) {
      return 0.0;
    }
    return rtn;
  }

  inline double V_FV(const double x)
  {
    double rtn = V_full(x);
    if(x>center_bar)
      rtn += barrier_strength*(x-center_bar)*(x-center_bar);
    return rtn;
  }

  inline double dV_FV(const double x)
  {
    double rtn = dV_full(x);
    if(x>center_bar)
      rtn += 2.0*barrier_strength*(x-center_bar);
    return rtn;
  }

  inline double V_TV(const double x)
  {
    double rtn = V_full(x);
    if(x<center_bar)
      rtn += M*barrier_strength*(x-center_bar)*(x-center_bar);
    else
      rtn += L*barrier_strength*(x-center_bar)*(x-center_bar);
    return rtn;
  }

  inline double dV_TV(const double x)
  {
    double rtn = dV_full(x);
    if(x<center_bar)
      rtn += M*2.0*barrier_strength*(x-center_bar);
    else
      rtn += L*2.0*barrier_strength*(x-center_bar);
    return rtn;
  }

  qacc double action_point(QMAction& qma, const Field<double>& f, const Geometry& geo, Coordinate xl)
  {
    // Returns the contribution to the total action from a single lattice
    // point (including the relavent neighbor interactions)
    // TIMER("QMAction.action_point");
    double psi = f.get_elem(xl);
    xl[3]+=1;
    double psi_eps = f.get_elem(xl);
    xl[3]-=1;
    return (beta/2.0/qma.dt/qma.dt)*(psi_eps-psi)*(psi_eps-psi) + qma.V(psi, geo.coordinate_g_from_l(xl)[3]);
  }

  inline double action_node_no_comm(const Field<double>& f)
  {
	// Returns the total action of the portion of the lattice on the
	// current node (assuming the necessary communication has already
	// been done)
    TIMER("QMAction.action_node_no_comm");
    const Geometry geo = f.geo();
    // Creates a geometry that is the same as the field geometry, except
    // with multiplicity 1
    const Geometry geo_r = geo_resize(geo);
    // Creates a field to save the contribution to the total action from
    // each point
    FieldM<double, 1> fd;
    fd.init(geo_r);
    // Loops over every lattice point in the current node
    QMAction& qma = *this;
    qacc_for(index, geo_r.local_volume(), {
      const Geometry& geo = f.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      fd.get_elem(index) = qma.action_point(qma, f, geo, xl);
    });
    // Sums over the contributions to the total action from each point
    // (this cannot be done in the previous loops because the previous
    // loop runs in parallel)
    double sum = 0.0;
    for (Long index = 0; index < geo_r.local_volume(); ++index) {
      sum += fd.get_elem(index);
    }
    return sum*dt;
  }

  inline double action_node(const Field<double>& f)
  {
	// Return the total Euclidean action (on the current node) associated
	// with the given scalar field.
    TIMER("QMAction.action_node");
    // Get the values of the field sites that lie just outside of the
    // current node in the time direction. This is needed to compute the
    // derivative term.
    const Coordinate expand_left(0, 0, 0, 0);
    const Coordinate expand_right(0, 0, 0, 1);
    const Geometry geo_ext = geo_resize(f.geo(), expand_left, expand_right);
    Field<double> f_ext;
    f_ext.init(geo_ext, f.multiplicity);
    f_ext = f;
    refresh_expanded(f_ext);
    // Calculate the energy of the scalar field
    return action_node_no_comm(f_ext);
  }

  inline double hmc_m_hamilton_node(const Field<double>& m)
  {
    // Return the part of an HMC Hamiltonian due to the given momentum
    // field (on the current node).
    TIMER("QMAction.hmc_m_hamilton_node");
    double sum = sum_sq(m);
    return sum/2.0;
  }

  inline double sum_sq(const Field<double>& f)
  {
    // Returns the sum of f(x)^2 over lattice sites (on the current
    // node) and multiplicity
    TIMER("QMAction.sum_sq");
    const Geometry geo = f.geo();
    // Creates a geometry that is the same as the field geometry, except
    // with multiplicity 1
    const Geometry geo_r = geo_resize(geo);
    // Creates a field to save the contribution to the sum of squares
    // from each point
    FieldM<double, 1> fd;
    fd.init(geo_r);
    qacc_for(index, geo_r.local_volume(), {
      const Geometry& geo_r = fd.geo();
      Coordinate xl = geo_r.coordinate_from_index(index);
      double s=0;
      for (int m = 0; m < f.multiplicity; ++m) {
        double d = f.get_elem(xl,m);
        s += d*d;
      }
      fd.get_elem(index) = s;
    });
    // Sums over the contributions to the sum of squares from each point
    // (this cannot be done in the previous loops because the previous
     // loop runs in parallel)
    double sum = 0;
    for (Long index = 0; index < geo_r.local_volume(); ++index) {
      sum += fd.get_elem(index);
    }
    return sum;
  }

  inline void hmc_set_force_no_comm(Field<double>& force, const Field<double>& f)
  {
    TIMER("QMAction.hmc_set_sm_force_no_comm");
    const Geometry geo = f.geo();
    QMAction& qma = *this;
    qacc_for(index, geo.local_volume(), {
      const Geometry& geo = f.geo();
      Coordinate xl = geo.coordinate_from_index(index);
      Vector<double> force_v = force.get_elems(xl);
      qassert(force_v.size() == 1);
      double psi = f.get_elem(xl);
      force_v[0] = 2.0 * qma.beta / qma.dt / qma.dt * psi;
      force_v[0] += qma.dV(psi,geo.coordinate_g_from_l(xl)[3]);
      xl[3] += 1;
      force_v[0] -= qma.beta / qma.dt / qma.dt * f.get_elem(xl);
      xl[3] -= 2;
      force_v[0] -= qma.beta / qma.dt / qma.dt * f.get_elem(xl);
      xl[3] += 1;
      force_v[0] *= qma.dt;
    });
  }

  inline void hmc_set_force(Field<double>&  force, const Field<double>&  f)
  {
	// Calculate and set the HMC force field based on the given field
	// configuration.
    TIMER("QMAction.hmc_set_sm_force");
    Coordinate expand_left(0, 0, 0, 1);
    Coordinate expand_right(0, 0, 0, 1);
    const Geometry geo_ext = geo_resize(f.geo(), expand_left, expand_right);
    Field<double> f_ext;
    f_ext.init(geo_ext, f.multiplicity);
    f_ext = f;
    refresh_expanded(f_ext);
    hmc_set_force_no_comm(force, f_ext);
  }

  inline void hmc_field_evolve(Field<double>& f, const Field<double>& m,
                               const double step_size)
  {
    TIMER("QMAction.hmc_field_evolve");
    const Geometry& geo = f.geo();
    qacc_for(index, geo.local_volume(), {
      const Geometry& geo = f.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      Vector<double> f_v = f.get_elems(xl);
      qassert(f_v.size() == 1);
      f_v[0] = f_v[0] + m.get_elem(xl)*step_size;
    });
  }

  inline void hmc_set_rand_momentum(Field<double>& m, const RngState& rs)
  {
    TIMER("QMAction.set_rand_momentum");
    const Geometry& geo = m.geo();
    qthread_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const Long gindex = geo.g_index_from_g_coordinate(xg);
      //displayln_info(ssprintf("gindex: %ld", gindex));
      RngState rsi = rs.newtype(gindex);
      m.get_elem(xl) = g_rand_gen(rsi, 0, 1.0);
      //displayln_info(ssprintf("random num: %24.17E", m.get_elem(xl)));
    });
  }
};

}  // namespace qlat
