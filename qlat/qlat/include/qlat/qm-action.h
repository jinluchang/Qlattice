#pragma once

namespace qlat
{

struct QMAction {
  bool initialized;
  RealD alpha;
  RealD beta;
  RealD start_TV;
  RealD FV_offset;
  RealD TV_offset;
  RealD center_bar;
  RealD barrier_strength;
  RealD L;
  RealD M;
  RealD epsilon;
  RealD temp;
  Long t_FV_out;
  Long t_FV_mid;
  RealD dt;
  bool measure_offset_L;
  bool measure_offset_M;
  
  enum VType { GET_L, GET_M, GET_DTTV };
  VType vtype;
  
  qacc void init()
  {
    initialized = false;
    alpha = 1.0;
    beta = 1.0;
    start_TV = 2.0;
    FV_offset = 0.0;
    TV_offset = 0.0;
    center_bar = 1.0;
    barrier_strength = 1.0;
    L = 1.0;
    M = 0.0;
    epsilon = 0.0;
    t_FV_out = 10;
    t_FV_mid = 5;
    dt = 1.0;
    measure_offset_L = false;
    measure_offset_M = false;
  }
  //
  qacc QMAction(const RealD alpha_, const RealD beta_, const RealD FV_offset_,
                const RealD TV_offset_, const RealD barrier_strength_, const RealD L_,
                const RealD M_, const RealD epsilon_, const Long t_FV_out_, 
                const Long t_FV_mid_, const RealD dt_,
                const bool measure_offset_L_, const bool measure_offset_M_)
  {
    init();
    initialized = true;
    alpha = alpha_;
    beta = beta_;
    start_TV = (3.0+std::pow(9.0-8.0*alpha, 0.5))/2.0/alpha; // (2.0-2.0*std::pow(1-alpha, 0.5))/alpha + start_TV_;
    FV_offset = FV_offset_;
    TV_offset = TV_offset_;
    center_bar = (3.0-std::pow(9.0-8.0*alpha, 0.5))/2.0/alpha;
    barrier_strength = barrier_strength_;
    L = L_;
    M = M_;
    epsilon = epsilon_;
    t_FV_out = t_FV_out_;
    t_FV_mid = t_FV_mid_;
    dt = dt_;
    measure_offset_L = measure_offset_L_;
    measure_offset_M = measure_offset_M_;
    
    if(L == 0.0 && M == 0.0)
      vtype = GET_DTTV; // Calculate observables for Q(t_TV+1) / Q(t_TV)
    else if(L == 1.0)
        vtype = GET_M; // Calculate Q(M_{i+1}) / Q(M_i)
    else
      vtype = GET_L; // Calculate Q(L_{i+1}) / Q(L_i)
  }
  
  inline RealD V(const RealD x, const Long t)
  {
    switch(vtype) {
      case GET_L: return V_t_L(x,t);
      case GET_M: return V_t_M(x,t);
      default: return V_t_dtTV(x,t);
    }
  }
  
  inline RealD dV(const RealD x, const Long t)
  {
    switch(vtype) {
      case GET_L: return dV_t_L(x,t);
      case GET_M: return dV_t_M(x,t);
      default: return dV_t_dtTV(x,t);
    }
  }
  
  inline RealD V_t_L(const RealD x, const Long t)
  {
    // Returns the potential evaluated at point x
    // Start with H_proj
    if(t==0)
      return V_L(x, V_full(x) + V_proj(x));
    // Until t_FV_out has past, use H_FV_out
    else if(t<=t_FV_out)
      return V_FV_out(x);
    // Until t_FV_mid has past, use H_FV_mid
    else if(t<=t_FV_out+t_FV_mid)
      return V_FV_mid(x);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=2*t_FV_out+t_FV_mid)
      return V_FV_out(x);
    // Right after t_FV_out has past, use H_proj
    else if(t==2*t_FV_out+t_FV_mid+1)
      return V_L(x, V_full(x) + V_proj(x));
    // For the rest of the time, use H_TV
    else {
      return V_L(x, V_TV(x));
    }
  }
  
  inline RealD dV_t_L(const RealD x, const Long t)
  {
    // Returns the potential evaluated at point x
    // Start with H_proj
    if(t==0)
      return dV_L(x, V_full(x) + V_proj(x), dV_full(x) + dV_proj(x));
    // Until t_FV_out has past, use H_FV_out
    else if(t<=t_FV_out)
      return dV_FV_out(x);
    // Until t_FV_mid has past, use H_FV_mid
    else if(t<=t_FV_out+t_FV_mid)
      return dV_FV_mid(x);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=2*t_FV_out+t_FV_mid)
      return dV_FV_out(x);
    // Right after t_FV_out has past, use H_proj
    else if(t==2*t_FV_out+t_FV_mid+1)
      return dV_L(x, V_full(x) + V_proj(x), dV_full(x) + dV_proj(x));
    // For the rest of the time, use H_TV
    else {
      return dV_L(x, V_TV(x), dV_TV(x));
    }
  }
  
  inline RealD V_t_M(const RealD x, const Long t)
  {
    // Returns the potential evaluated at point x
    // Start with H_proj
    if(t==0)
      return V_M(x, V_full(x) + V_proj(x));
    // Until t_FV_out has past, use H_FV_out
    else if(t<=t_FV_out)
      return V_FV_out(x);
    // Until t_FV_mid has past, use H_FV_mid
    else if(t<=t_FV_out+t_FV_mid)
      return V_FV_mid(x);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=2*t_FV_out+t_FV_mid)
      return V_FV_out(x);
    // Right after t_FV_out has past, use H_proj
    else if(t==2*t_FV_out+t_FV_mid+1)
      return V_M(x, V_full(x) + V_proj(x));
    // For the rest of the time, use H_TV
    else {
      return V_M(x, V_TV(x));
    }
  }
  
  inline RealD dV_t_M(const RealD x, const Long t)
  {
    // Returns the potential evaluated at point x
    // Start with H_proj
    if(t==0)
      return dV_M(x, V_full(x) + V_proj(x), dV_full(x) + dV_proj(x));
    // Until t_FV_out has past, use H_FV_out
    else if(t<=t_FV_out)
      return dV_FV_out(x);
    // Until t_FV_mid has past, use H_FV_mid
    else if(t<=t_FV_out+t_FV_mid)
      return dV_FV_mid(x);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=2*t_FV_out+t_FV_mid)
      return dV_FV_out(x);
    // Right after t_FV_out has past, use H_proj
    else if(t==2*t_FV_out+t_FV_mid+1)
      return dV_M(x, V_full(x) + V_proj(x), dV_full(x) + dV_proj(x));
    // For the rest of the time, use H_TV
    else {
      return dV_M(x, V_TV(x), dV_TV(x));
    }
  }
  
  inline RealD V_t_dtTV(const RealD x, const Long t)
  {
    // Returns the potential evaluated at point x
    // Start with H_proj
    if(t==0)
      return V_full(x) + V_proj(x);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=t_FV_out)
      return V_FV_out(x);
    // Until t_FV_mid has past, use H_FV_mid
    else if(t<=t_FV_out+t_FV_mid)
      return V_FV_mid(x);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=2*t_FV_out+t_FV_mid)
      return V_FV_out(x);
    // Right after t_FV_out has past, use H_proj or V_max
    else if(t==2*t_FV_out+t_FV_mid+1) {
      if(measure_offset_L) return V_max(V_full(x) + V_proj(x), V_FV_out(x), 1.0);
      if(measure_offset_M) return V_max(V_FV_out(x), V_full(x) + V_proj(x), 1.0);
      else return V_full(x) + V_proj(x);
    }
    // One timeslice after t_FV_out has past, use H_TV or V_max
    else if(t==2*t_FV_out+t_FV_mid+2) {
      if(measure_offset_L) return V_max(V_TV(x), V_full(x) + V_proj(x), 1.0);
      if(measure_offset_M) return V_max(V_full(x) + V_proj(x), V_TV(x), 1.0);
      else return V_TV(x);
    }
    // For the rest of the time, use H_TV
    else {
      return V_TV(x);
    }
  }
  
  inline RealD dV_t_dtTV(const RealD x, const Long t)
  {
    // Returns the potential evaluated at point x
    // Start with H_proj
    if(t==0)
      return dV_full(x) + dV_proj(x);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=t_FV_out)
      return dV_FV_out(x);
    // Until t_FV_mid has past, use H_FV_mid
    else if(t<=t_FV_out+t_FV_mid)
      return dV_FV_mid(x);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=2*t_FV_out+t_FV_mid)
      return dV_FV_out(x);
    // Right after t_FV_out has past, use H_proj or V_max
    else if(t==2*t_FV_out+t_FV_mid+1) {
      if(measure_offset_L) return dV_max(V_full(x) + V_proj(x), dV_full(x) + dV_proj(x), V_FV_out(x), dV_FV_out(x), 1.0);
      if(measure_offset_M) return dV_max(V_FV_out(x), dV_FV_out(x), V_full(x) + V_proj(x), dV_full(x) + dV_proj(x), 1.0);
      else return dV_full(x) + dV_proj(x);
    }
    // One timeslice after t_FV_out has past, use H_TV or V_max
    else if(t==2*t_FV_out+t_FV_mid+2) {
      if(measure_offset_L) return dV_max(V_TV(x), dV_TV(x), V_full(x) + V_proj(x), dV_full(x) + dV_proj(x), 1.0);
      if(measure_offset_M) return dV_max(V_full(x) + V_proj(x), dV_full(x) + dV_proj(x), V_TV(x), dV_TV(x), 1.0);
      else return dV_TV(x);
    }
    // For the rest of the time, use H_TV
    else {
      return dV_TV(x);
    }
  }
  
  inline RealD V_zeroed(const RealD x, const Long t) {
    return V(x,t) - log(dt) / dt;
  }

  inline RealD V_phi4(const RealD x)
  {
    return beta*(x*x/2.0-x*x*x/2.0+alpha*x*x*x*x/8.0);
  }

  inline RealD dV_phi4(const RealD x)
  {
    // Returns the derivative of the potential with respect to x
    return beta*(x-x*x*3.0/2.0+alpha*x*x*x/2.0);
  }

  inline RealD V_full(const RealD x)
  {
    // Returns the potential evaluated at point x
    if(x>start_TV) {
      return 0;
    }
    return V_phi4(x) - V_phi4(start_TV);
  }

  inline RealD dV_full(const RealD x)
  {
    RealD rtn = dV_phi4(x);
    if(x>start_TV) {
      return 0.0;
    }
    return rtn;
  }

  inline RealD V_FV_out(const RealD x)
  {
    if(x>center_bar+FV_offset)
        return V_full(center_bar+FV_offset) + barrier_strength*(x-center_bar-FV_offset)*(x-center_bar-FV_offset);
    else
        return V_full(x);
  }

  inline RealD dV_FV_out(const RealD x)
  {
    if(x>center_bar+FV_offset)
        return 2.0*barrier_strength*(x-center_bar-FV_offset);
    else
        return dV_full(x);
  }

  inline RealD V_FV_mid(const RealD x)
  {
    if(x>center_bar)
      return V_full(center_bar) + barrier_strength*(x-center_bar)*(x-center_bar);
    return V_full(x);
  }

  inline RealD dV_FV_mid(const RealD x)
  {
    if(x>center_bar)
      return 2.0*barrier_strength*(x-center_bar);
    return dV_full(x);
  }

  inline RealD V_TV(const RealD x)
  {
    if(x < center_bar+TV_offset)
      return V_full(center_bar+TV_offset) + barrier_strength*(x-center_bar-TV_offset)*(x-center_bar-TV_offset);
    else
      return V_full(x);
  }

  inline RealD dV_TV(const RealD x)
  {
    if(x < center_bar+TV_offset)
      return 2.0*barrier_strength*(x-center_bar-TV_offset);
    else
      return dV_full(x);
  }
  
  inline RealD V_proj(const RealD x)
  {
    RealD rtn = -log((1-exp(-(V_FV_out(x) - V_full(x) + epsilon)*dt)) / dt) / dt;
    // When x is low enough that epsilon is relevant, remove V_full (which
    // will be added later) to avoid ergodicity issues
    if(x<center_bar+FV_offset) rtn += V_full(center_bar+FV_offset) - V_full(x) + barrier_strength*std::pow(center_bar+FV_offset-x, 0.5);
    return rtn;
    
  }
  
  inline RealD dV_proj(const RealD x)
  {
    RealD Vbar = V_FV_out(x) - V_full(x);
    RealD rtn = -((dV_FV_out(x) - dV_full(x))*exp(-(Vbar + epsilon)*dt))/(1-exp(-(Vbar + epsilon)*dt));
    if(x<center_bar+FV_offset) rtn += - dV_full(x) - 0.5*barrier_strength/std::pow(center_bar+FV_offset-x, 0.5);
    return rtn;
  }
  
  inline RealD V_max(const RealD V_D, const RealD V_N, const RealD P)
  {
    if(V_N < V_D) return V_D;
    else return (1-P)*V_D + P*V_N;
  }
  
  inline RealD dV_max(const RealD V_D, const RealD dV_D, const RealD V_N, const RealD dV_N, const RealD P)
  {
    if(V_N < V_D) return dV_D;
    else return (1-P)*dV_D + P*dV_N;
  }
  
  inline RealD V_FV_floored(const RealD x, const RealD P)
  {
    RealD v_fv_mid = V_FV_mid(x);
    RealD v_fv_min = V_full(0);
    RealD floor = v_fv_min + P*(V_FV_mid(center_bar+FV_offset)-v_fv_min);
    if(v_fv_mid < floor) return floor;
    else return v_fv_mid;
  }
  
  inline RealD dV_FV_floored(const RealD x, const RealD P)
  {
    RealD v_fv_min = V_full(0);
    RealD floor = v_fv_min + P*(V_FV_mid(center_bar+FV_offset)-v_fv_min);
    if(V_FV_mid(x) < floor) return 0;
    else return dV_FV_mid(x);
  }
  
  inline RealD V_L(const RealD x, const RealD V_N)
  {
    if(L<0.5) {
      return V_FV_floored(x, 2*L);
    }
    else return V_max(V_FV_floored(x, 1), V_N, 2*(L-0.5));
  }
  
  inline RealD dV_L(const RealD x, const RealD V_N, const RealD dV_N)
  {
    if(L<0.5) {
      return dV_FV_floored(x, 2*L);
    }
    else return dV_max(V_FV_floored(x, 1), dV_FV_floored(x, 1), V_N, dV_N, 2*(L-0.5));
  }
  
  inline RealD V_M(const RealD x, const RealD V_N)
  {
    return V_max(V_N, V_FV_floored(x, 1), M);
  }
  
  inline RealD dV_M(const RealD x, const RealD V_N, const RealD dV_N)
  {
    return dV_max(V_N, dV_N, V_FV_floored(x, 1), dV_FV_floored(x, 1), M);
  }

  inline RealD action_point(QMAction& qma, const Field<RealD>& f, const Geometry& geo, Coordinate xl)
  {
    // Returns the contribution to the total action from a single lattice
    // point (including the relavent neighbor interactions)
    // TIMER("QMAction.action_point");
    RealD psi = f.get_elem(xl);
    xl[3]+=1;
    RealD psi_eps = f.get_elem(xl);
    xl[3]-=1;
    return (beta/2.0/qma.dt/qma.dt)*(psi_eps-psi)*(psi_eps-psi) + qma.V_zeroed(psi, geo.coordinate_g_from_l(xl)[3]);
  }

  inline RealD action_node_no_comm(const Field<RealD>& f)
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
    FieldM<RealD, 1> fd;
    fd.init(geo_r);
    // Loops over every lattice point in the current node
    QMAction& qma = *this;
    qthread_for(index, geo_r.local_volume(), {
      const Geometry& geo = f.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      fd.get_elem(index) = qma.action_point(qma, f, geo, xl);
    });
    // Sums over the contributions to the total action from each point
    // (this cannot be done in the previous loops because the previous
    // loop runs in parallel)
    RealD sum = 0.0;
    for (Long index = 0; index < geo_r.local_volume(); ++index) {
      sum += fd.get_elem(index);
    }
    return sum*dt;
  }

  inline RealD action_node(const Field<RealD>& f)
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
    Field<RealD> f_ext;
    f_ext.init(geo_ext, f.multiplicity);
    f_ext = f;
    refresh_expanded(f_ext);
    // Calculate the energy of the scalar field
    return action_node_no_comm(f_ext);
  }

  inline RealD hmc_m_hamilton_node(const Field<RealD>& m)
  {
    // Return the part of an HMC Hamiltonian due to the given momentum
    // field (on the current node).
    TIMER("QMAction.hmc_m_hamilton_node");
    RealD sum = sum_sq(m);
    return sum/2.0;
  }

  inline RealD sum_sq(const Field<RealD>& f)
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
    FieldM<RealD, 1> fd;
    fd.init(geo_r);
    qthread_for(index, geo_r.local_volume(), {
      const Geometry& geo_r = fd.geo();
      Coordinate xl = geo_r.coordinate_from_index(index);
      RealD s=0;
      for (Int m = 0; m < f.multiplicity; ++m) {
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

  inline void hmc_set_force_no_comm(Field<RealD>& force, const Field<RealD>& f)
  {
    TIMER("QMAction.hmc_set_sm_force_no_comm");
    const Geometry geo = f.geo();
    QMAction& qma = *this;
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = f.geo();
      Coordinate xl = geo.coordinate_from_index(index);
      Vector<RealD> force_v = force.get_elems(xl);
      qassert(force_v.size() == 1);
      RealD psi = f.get_elem(xl);
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

  inline void hmc_set_force(Field<RealD>&  force, const Field<RealD>&  f)
  {
	// Calculate and set the HMC force field based on the given field
	// configuration.
    TIMER("QMAction.hmc_set_sm_force");
    Coordinate expand_left(0, 0, 0, 1);
    Coordinate expand_right(0, 0, 0, 1);
    const Geometry geo_ext = geo_resize(f.geo(), expand_left, expand_right);
    Field<RealD> f_ext;
    f_ext.init(geo_ext, f.multiplicity);
    f_ext = f;
    refresh_expanded(f_ext);
    hmc_set_force_no_comm(force, f_ext);
  }

  inline void hmc_field_evolve(Field<RealD>& f, const Field<RealD>& m,
                               const RealD step_size)
  {
    TIMER("QMAction.hmc_field_evolve");
    const Geometry& geo = f.geo();
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = f.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      Vector<RealD> f_v = f.get_elems(xl);
      qassert(f_v.size() == 1);
      f_v[0] = f_v[0] + m.get_elem(xl)*step_size;
    });
  }

  inline void hmc_set_rand_momentum(Field<RealD>& m, const RngState& rs)
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
