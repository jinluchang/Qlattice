#pragma once

namespace qlat
{

struct QMAction {
  bool initialized;
  RealD alpha;
  RealD beta;
  //RealD start_TV;
  RealD FV_offset;
  RealD TV_offset;
  RealD center_bar;
  RealD V_min_FV;
  RealD barrier_strength;
  RealD L;
  RealD M;
  RealD epsilon;
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
    //start_TV = 2.0;
    FV_offset = 0.0;
    TV_offset = 0.0;
    center_bar = 1.0;
    V_min_FV = 0.0;
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
  qacc QMAction(const RealD alpha_, const RealD beta_, const RealD V_FV_min_,
                const RealD FV_offset_, const RealD TV_offset_,
                const RealD barrier_strength_, const RealD L_,
                const RealD M_, const RealD epsilon_, const Long t_FV_out_, 
                const Long t_FV_mid_, const RealD dt_,
                const bool measure_offset_L_, const bool measure_offset_M_)
  {
    init();
    initialized = true;
    alpha = alpha_;
    beta = beta_;
    //start_TV = (3.0+std::pow(9.0-8.0*alpha, 0.5))/2.0/alpha;
    V_min_FV = V_FV_min_;//V_phi4(start_TV,0);
    FV_offset = FV_offset_;
    TV_offset = TV_offset_;
    // center_bar = (3.0-std::pow(9.0-8.0*alpha, 0.5))/2.0/alpha;
    center_bar = 3.0/(alpha*alpha);
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
  
  inline RealD V(const Vector<RealD>& x, const Long t)
  {
    switch(vtype) {
      case GET_L: return V_t_L(x,t);
      case GET_M: return V_t_M(x,t);
      default: return V_t_dtTV(x,t);
    }
  }
  
  inline RealD dV(const Vector<RealD>& x, const Long t, const int idx)
  {
    switch(vtype) {
      case GET_L: return dV_t_L(x,t,idx);
      case GET_M: return dV_t_M(x,t,idx);
      default: return dV_t_dtTV(x,t,idx);
    }
  }
  
  inline RealD V_t_L(const Vector<RealD>& x, const Long t)
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
  
  inline RealD dV_t_L(const Vector<RealD>& x, const Long t, const int idx)
  {
    // Returns the potential evaluated at point x
    // Start with H_proj
    if(t==0)
      return dV_L(x, V_full(x) + V_proj(x), dV_full(x,idx) + dV_proj(x,idx), idx);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=t_FV_out)
      return dV_FV_out(x, idx);
    // Until t_FV_mid has past, use H_FV_mid
    else if(t<=t_FV_out+t_FV_mid)
      return dV_FV_mid(x, idx);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=2*t_FV_out+t_FV_mid)
      return dV_FV_out(x, idx);
    // Right after t_FV_out has past, use H_proj
    else if(t==2*t_FV_out+t_FV_mid+1)
      return dV_L(x, V_full(x) + V_proj(x), dV_full(x,idx) + dV_proj(x,idx), idx);
    // For the rest of the time, use H_TV
    else {
      return dV_L(x, V_TV(x), dV_TV(x,idx), idx);
    }
  }
  
  inline RealD V_t_M(const Vector<RealD>& x, const Long t)
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
  
  inline RealD dV_t_M(const Vector<RealD>& x, const Long t, const int idx)
  {
    // Returns the potential evaluated at point x
    // Start with H_proj
    if(t==0)
      return dV_M(x, V_full(x) + V_proj(x), dV_full(x,idx) + dV_proj(x,idx), idx);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=t_FV_out)
      return dV_FV_out(x, idx);
    // Until t_FV_mid has past, use H_FV_mid
    else if(t<=t_FV_out+t_FV_mid)
      return dV_FV_mid(x, idx);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=2*t_FV_out+t_FV_mid)
      return dV_FV_out(x, idx);
    // Right after t_FV_out has past, use H_proj
    else if(t==2*t_FV_out+t_FV_mid+1)
      return dV_M(x, V_full(x) + V_proj(x), dV_full(x,idx) + dV_proj(x,idx), idx);
    // For the rest of the time, use H_TV
    else {
      return dV_M(x, V_TV(x), dV_TV(x,idx), idx);
    }
  }
  
  inline RealD V_t_dtTV(const Vector<RealD>& x, const Long t)
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
  
  inline RealD dV_t_dtTV(const Vector<RealD>& x, const Long t, const int idx)
  {
    // Returns the potential evaluated at point x
    // Start with H_proj
    if(t==0)
      return dV_full(x,idx) + dV_proj(x,idx);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=t_FV_out)
      return dV_FV_out(x,idx);
    // Until t_FV_mid has past, use H_FV_mid
    else if(t<=t_FV_out+t_FV_mid)
      return dV_FV_mid(x,idx);
    // Until t_FV_out has past, use H_FV_out
    else if(t<=2*t_FV_out+t_FV_mid)
      return dV_FV_out(x,idx);
    // Right after t_FV_out has past, use H_proj or V_max
    else if(t==2*t_FV_out+t_FV_mid+1) {
      if(measure_offset_L) return dV_max(V_full(x) + V_proj(x), dV_full(x,idx) + dV_proj(x,idx), V_FV_out(x), dV_FV_out(x,idx), 1.0);
      if(measure_offset_M) return dV_max(V_FV_out(x), dV_FV_out(x,idx), V_full(x) + V_proj(x), dV_full(x,idx) + dV_proj(x,idx), 1.0);
      else return dV_full(x,idx) + dV_proj(x,idx);
    }
    // One timeslice after t_FV_out has past, use H_TV or V_max
    else if(t==2*t_FV_out+t_FV_mid+2) {
      if(measure_offset_L) return dV_max(V_TV(x), dV_TV(x,idx), V_full(x) + V_proj(x), dV_full(x,idx) + dV_proj(x,idx), 1.0);
      if(measure_offset_M) return dV_max(V_full(x) + V_proj(x), dV_full(x,idx) + dV_proj(x,idx), V_TV(x), dV_TV(x,idx), 1.0);
      else return dV_TV(x,idx);
    }
    // For the rest of the time, use H_TV
    else return dV_TV(x, idx);
  }
  
  inline RealD V_zeroed(const Vector<RealD>& x, const Long t) {
    return V(x,t) - log(dt) / dt;
  }
  
  // Quantum Mechanics Potential========================================
  inline RealD order_param(const Vector<RealD> x)
  {
    if(x[0]<0) return x[1]*x[1];
    else return x[0]*x[0] + x[1]*x[1];
  }
  
  inline RealD d_order_param(const Vector<RealD> x, const int idx)
  {
    if (idx==1) return 2.0*x[1];
    else if (x[0]<0) return 0;
    else return 2.0*x[0];
  }

  inline RealD V_phi4(const RealD x, const RealD y)
  {
    // return beta*(x*x/2.0-x*x*x/2.0+alpha*x*x*x*x/8.0);
    // return beta*(x*x/2.0+y*y/2.0-x*x*x/2.0+alpha*x*x*x*x/8.0);
    return beta*(x*x+(1-alpha*x)*y*y);
  }

  inline RealD dV_phi4(const RealD x, const RealD y, const int idx)
  {
    // Adds the derivative of the potential with respect to x
    // if (idx==0) return beta*(x - 1.5*x*x + alpha*x*x*x/2.0);
    // else return beta*y;
    if (idx==0) return beta*(2.0*x-alpha*y*y);
    else return beta*2.0*(1-alpha*x)*y;
  }

  inline RealD V_full_xy(const RealD x, const RealD y)
  {
    const RealD rtn = V_phi4(x,y);
    if(rtn<V_min_FV) return 0;
    else return rtn - V_min_FV;
  }

  inline RealD dV_full_xy(const RealD x, const RealD y, const int idx)
  {
    if(V_phi4(x,y)<V_min_FV) return 0;
    else return dV_phi4(x,y,idx);
  }

  inline RealD V_full_op_fixed(const Vector<RealD>& x, RealD op)
  {
    const RealD norm = std::pow(order_param(x)/op,0.5);
    if(norm==0) return V_full_xy(op, 0);
    return V_full_xy(x[0]/norm, x[1]/norm); // V(div(x))
  }

  inline RealD dV_full_op_fixed(const Vector<RealD>& x, RealD op, const int idx)
  {
    const RealD norm = std::pow(order_param(x)/op,0.5); // div(x) = op^0.5*x_i / order_param(x)^0.5
    if(norm==0) return 0;
    return (1/norm - 0.5*x[idx]*d_order_param(x,idx)/norm/order_param(x)) * dV_full_xy(x[0]/norm, x[1]/norm, idx); // ddiv(x) * dV(div(x))
  }

  inline RealD V_full(const Vector<RealD>& x)
  {
    return V_full_xy(x[0], x[1]);
  }

  inline RealD dV_full(const Vector<RealD>& x, const int idx)
  {
    return dV_full_xy(x[0], x[1], idx);
  }
  
  // End Quantum Mechanics Potential====================================

  inline RealD V_FV_out(const Vector<RealD>& x)
  {
    if(order_param(x)>center_bar+FV_offset)
        return V_full_op_fixed(x,center_bar+FV_offset) + barrier_strength*(order_param(x)-center_bar-FV_offset)*(order_param(x)-center_bar-FV_offset);
    else
        return V_full(x);
  }

  inline RealD dV_FV_out(const Vector<RealD>& x, const int idx)
  {
    if(order_param(x)>center_bar+FV_offset) {
        return dV_full_op_fixed(x,center_bar+FV_offset,idx) + d_order_param(x,idx)*2.0*barrier_strength*(order_param(x)-center_bar-FV_offset);
    }
    else return dV_full(x,idx);
  }

  inline RealD V_FV_mid(const Vector<RealD>& x)
  {
    if(order_param(x)>center_bar)
      return V_full_op_fixed(x,center_bar) + barrier_strength*(order_param(x)-center_bar)*(order_param(x)-center_bar);
    else
      return V_full(x);
  }

  inline RealD dV_FV_mid(const Vector<RealD>& x, const int idx)
  {
    if(order_param(x)>center_bar) {
      return dV_full_op_fixed(x,center_bar,idx) + d_order_param(x,idx)*2.0*barrier_strength*(order_param(x)-center_bar);
    }
    else return dV_full(x,idx);
  }

  inline RealD V_FV_mid_op_fixed(const Vector<RealD>& x, RealD op)
  {
    if(op>center_bar)
      return V_full_op_fixed(x,center_bar) + barrier_strength*(op-center_bar)*(op-center_bar);
    else
      return V_full_op_fixed(x,op);
  }

  inline RealD dV_FV_mid_op_fixed(const Vector<RealD>& x, RealD op, const int idx)
  {
    if(op>center_bar)
      return dV_full_op_fixed(x,center_bar,idx);
    else
      return dV_full_op_fixed(x,op,idx);
  }

  inline RealD V_TV(const Vector<RealD>& x)
  {
    if(order_param(x) < center_bar+TV_offset)
      return V_full_op_fixed(x,center_bar+TV_offset) + barrier_strength*(order_param(x)-center_bar-TV_offset)*(order_param(x)-center_bar-TV_offset);
    else
      return V_full(x);
  }

  inline RealD dV_TV(const Vector<RealD>& x, const int idx)
  {
    if(order_param(x) < center_bar+TV_offset) {
      return dV_full_op_fixed(x,center_bar+TV_offset,idx) + d_order_param(x,idx)*2.0*barrier_strength*(order_param(x)-center_bar-TV_offset);
    }
    else return dV_full(x,idx);
  }
  
  inline RealD V_proj(const Vector<RealD>& x)
  {
    RealD rtn = -log((1-exp(-(V_FV_out(x) - V_full(x) + epsilon)*dt)) / dt) / dt;
    // When x is low enough that epsilon is relevant, remove V_full (which
    // will be added later) to avoid ergodicity issues
    if(order_param(x)<center_bar+FV_offset) {
      rtn += V_full_op_fixed(x, center_bar+FV_offset) - V_full(x) + barrier_strength*std::pow(center_bar+FV_offset-order_param(x), 0.5);
    }
    return rtn;
  }
  
  inline RealD dV_proj(const Vector<RealD>& x, const int idx)
  {
    RealD Vbar = V_FV_out(x) - V_full(x);
    RealD rtn = -((dV_FV_out(x,idx) - dV_full(x,idx))*exp(-(Vbar + epsilon)*dt))/(1-exp(-(Vbar + epsilon)*dt));
    if(order_param(x)<center_bar+FV_offset) {
      rtn += dV_full_op_fixed(x,center_bar+FV_offset,idx) - dV_full(x,idx) - d_order_param(x,idx)*0.5*barrier_strength/std::pow(center_bar+FV_offset-order_param(x), 0.5);
    }
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
  
  inline RealD V_FV_floored(const Vector<RealD>& x, const RealD P)
  {
    RealD v_fv_mid = V_FV_mid(x);
    RealD floor = V_min_FV + P*(V_FV_mid_op_fixed(x,center_bar+FV_offset)-V_min_FV);
    if(v_fv_mid < floor) return floor;
    else return v_fv_mid;
  }
  
  inline RealD dV_FV_floored(const Vector<RealD>& x, const RealD P, const int idx)
  {
    RealD floor = V_min_FV + P*(V_FV_mid_op_fixed(x,center_bar+FV_offset)-V_min_FV);
    if(V_FV_mid(x) < floor) return P*dV_FV_mid_op_fixed(x,center_bar+FV_offset,idx);
    else return dV_FV_mid(x, idx);
  }
  
  inline RealD V_L(const Vector<RealD>& x, const RealD V_N)
  {
    if(L<0.5) return V_FV_floored(x, 2*L);
    else return V_max(V_FV_floored(x, 1), V_N, 2*(L-0.5));
  }
  
  inline RealD dV_L(const Vector<RealD>& x, const RealD V_N, const RealD dV_N, const int idx)
  {
    if(L<0.5) return dV_FV_floored(x, 2*L, idx);
    else return dV_max(V_FV_floored(x, 1), dV_FV_floored(x, 1, idx), V_N, dV_N, 2*(L-0.5));
  }
  
  inline RealD V_M(const Vector<RealD>& x, const RealD V_N)
  {
    return V_max(V_N, V_FV_floored(x, 1), M);
  }
  
  inline RealD dV_M(const Vector<RealD>& x, const RealD V_N, const RealD dV_N, const int idx)
  {
    return dV_max(V_N, dV_N, V_FV_floored(x, 1), dV_FV_floored(x, 1, idx), M);
  }

  inline RealD action_point(QMAction& qma, Field<RealD>& f, const Int multiplicity, const Geometry& geo, Coordinate xl)
  {
    // Returns the contribution to the total action from a single lattice
    // point (including the relavent neighbor interactions)
    // TIMER("QMAction.action_point");
    Vector<RealD> psi = f.get_elems(xl);
    xl[3]+=1;
    Vector<RealD> psi_eps = f.get_elems(xl);
    xl[3]-=1;
    RealD dpsi_sq = 0;
    for (Int i = 0; i < multiplicity; ++i) {
      dpsi_sq += (psi_eps[i]-psi[i])*(psi_eps[i]-psi[i]);
    }
    return (beta/2.0/qma.dt/qma.dt)*dpsi_sq + qma.V_zeroed(psi, geo.coordinate_g_from_l(xl)[3]);
  }

  inline RealD action_node_no_comm(Field<RealD>& f)
  {
	// Returns the total action of the portion of the lattice on the
	// current node (assuming the necessary communication has already
	// been done)
    TIMER("QMAction.action_node_no_comm");
    const Geometry geo = f.geo();
    const Int multiplicity = f.multiplicity;
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
      fd.get_elem(index) = qma.action_point(qma, f, multiplicity, geo, xl);
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

  inline RealD action_node(Field<RealD>& f)
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

  inline RealD hmc_m_hamilton_node(const Field<RealD>& mf)
  {
    // Return the part of an HMC Hamiltonian due to the given momentum
    // field (on the current node).
    TIMER("QMAction.hmc_m_hamilton_node");
    RealD sum = sum_sq(mf);
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

  inline void hmc_set_force_no_comm(Field<RealD>& force, Field<RealD>& f)
  {
    TIMER("QMAction.hmc_set_sm_force_no_comm");
    const Geometry geo = f.geo();
    const Int multiplicity = f.multiplicity;
    QMAction& qma = *this;
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = f.geo();
      Coordinate xl = geo.coordinate_from_index(index);
      Vector<RealD> force_v = force.get_elems(xl);
      qassert(force_v.size() == multiplicity);
      Vector<RealD> psi = f.get_elems(xl);
      for (Int i = 0; i < multiplicity; ++i) {
          force_v[i] = 2.0 * qma.beta / qma.dt / qma.dt * psi[i];
          xl[3] += 1;
          force_v[i] -= qma.beta / qma.dt / qma.dt * f.get_elem(xl, i);
          xl[3] -= 2;
          force_v[i] -= qma.beta / qma.dt / qma.dt * f.get_elem(xl, i);
          xl[3] += 1;
          force_v[i] += qma.dV(psi, geo.coordinate_g_from_l(xl)[3], i);
          force_v[i] *= qma.dt;
      }
    });
  }

  inline void hmc_set_force(Field<RealD>&  force, Field<RealD>&  f)
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

  inline void hmc_field_evolve(Field<RealD>& f, const Field<RealD>& mf,
                               const RealD step_size)
  {
    TIMER("QMAction.hmc_field_evolve");
    const Geometry& geo = f.geo();
    const Int multiplicity = f.multiplicity;
    qthread_for(index, geo.local_volume(), {
      const Geometry& geo = f.geo();
      const Coordinate xl = geo.coordinate_from_index(index);
      Vector<RealD> f_v = f.get_elems(xl);
      for (Int i = 0; i < multiplicity; ++i) {
        f_v[i] = f_v[i] + mf.get_elem(xl,i)*step_size;
      }
    });
  }

  inline void hmc_set_rand_momentum(Field<RealD>& mf, const RngState& rs)
  {
    TIMER("QMAction.set_rand_momentum");
    const Geometry& geo = mf.geo();
    const Int multiplicity = mf.multiplicity;
    qthread_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const Long gindex = geo.g_index_from_g_coordinate(xg);
      //displayln_info(ssprintf("gindex: %ld", gindex));
      RngState rsi = rs.newtype(gindex);
      Vector<RealD> m_v = mf.get_elems(xl);
      for (Int i = 0; i < multiplicity; ++i) {
        m_v[i] = g_rand_gen(rsi, 0, 1.0);
      }
    });
  }
};

}  // namespace qlat
