#ifndef QUDA_PARA_H
#define QUDA_PARA_H

#pragma once

#include <qlat/qlat.h>

#include <quda.h>
#include <invert_quda.h>

#include <cstdlib>
#include "utils_float_type.h"

namespace qlat
{  //

static int mpi_rank_from_coords_x(const int* coords, void* fdata)
{
  int* dims = reinterpret_cast<int*>(fdata);
  //
  int rank;
  rank = coords[3];
  for (int i = 2; i >= 0; i--) {
    rank = dims[i] * rank + coords[i];
  }
  return rank;
}
//
static int mpi_rank_from_coords_t(const int* coords, void* fdata)
{
  int* dims = reinterpret_cast<int*>(fdata);
  //
  int rank;
  rank = coords[0];
  for (int i = 1; i <= 3; i++) {
    rank = dims[i] * rank + coords[i];
  }
  return rank;
}
//
inline void quda_begin(int mpi_layout[4], bool t = false)
{
  using namespace quda;
  // The following sets the MPI comm stuff.
  if (t) {
    initCommsGridQuda(4, mpi_layout, mpi_rank_from_coords_t,
                      reinterpret_cast<void*>(mpi_layout));
  } else {
    initCommsGridQuda(4, mpi_layout, mpi_rank_from_coords_x,
                      reinterpret_cast<void*>(mpi_layout));
  }
  // comm_set_gridsize(mpi_layout);
  initQuda(-1000);
  //initQuda(-1);
  printf(
      "initialized on quda rank #%03d (%03d,%03d,%03d,%03d), qlat rank #%03d "
      "(%03d,%03d,%03d,%03d).\n",
      comm_rank(), comm_coord(0), comm_coord(1), comm_coord(2), comm_coord(3),
      get_id_node(), get_coor_node()[0], get_coor_node()[1], get_coor_node()[2],
      get_coor_node()[3]);
  // Make sure there is no mismatch
  qassert(comm_rank() == get_id_node());
  for (int d = 0; d < 4; d++) {
    qassert(comm_coord(d) == get_coor_node()[d]);
  }
}

inline void quda_end()
{
  using namespace quda;
  endQuda();
}

template <class T>
void quda_convert_gauge(qlat::vector<T>& qgf, GaugeField& gf, int dir = 0)
{
  TIMER("quda_convert_gauge(qgf,gf)");
  const Geometry& geo = gf.geo();
  ColorMatrix* quda_pt = reinterpret_cast<ColorMatrix*>(qgf.data());
  qassert(geo.multiplicity == 4);
  long V = geo.local_volume();
  long Vh = V / 2;
  for (long qlat_idx = 0; qlat_idx < V; qlat_idx++) {
    Coordinate xl = geo.coordinate_from_index(qlat_idx);
    //const Vector<ColorMatrix> ms = gf.get_elems_const(xl);
    //Vector<ColorMatrix> ms = gf.get_elems(xl);
    int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    for (int mu = 0; mu < 4; mu++) {
      ColorMatrixT<T>& ms = gf.get_elem(qlat_idx*gf.geo().multiplicity+mu);
      long quda_idx = (qlat_idx / 2 + eo * Vh) * 4 + mu;
      if(dir == 0){quda_pt[quda_idx] = ms;}
      if(dir == 1){ms = quda_pt[quda_idx];}
    }
  }
}

QudaPrecision get_quda_precision(int byte)
{
  switch (byte) {
    case 8:
      return QUDA_DOUBLE_PRECISION;
      break;
    case 4:
      return QUDA_SINGLE_PRECISION; 
      break;
    case 2:
      return QUDA_HALF_PRECISION; 
      break;
    default:
      qassert(false);
      return QUDA_INVALID_PRECISION;
  }
}


//template <class T>
//void quda_convert_fermion_copy(FermionField5d& ff, const std::vector<T>& qff)
//{
//  TIMER("quda_convert_fermion(ff,qff)");
//  const Geometry& geo = ff.geo();
//  const WilsonVector* quda_pt =
//      reinterpret_cast<const WilsonVector*>(qff.data());
//  int Ls = geo.multiplicity;
//  qassert(Ls > 0);
//  long V = geo.local_volume();
//  long Vh = V / 2;
//
//  //
//  #pragma omp parallel for
//  for (long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
//    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
//    int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
//    Vector<WilsonVector> wvs = ff.get_elems(xl);
//    for (int s = 0; s < Ls; s++) {
//      int quda_idx = eo * Vh * Ls + s * Vh + qlat_idx_4d / 2;
//      wvs[s] = quda_pt[quda_idx];
//    }
//  }
//
//}

//template <class T>
//void quda_convert_fermion_copy(std::vector<T>& qff, const FermionField5d& ff)
//{
//  TIMER("quda_convert_fermion(qff,ff)");
//  const Geometry& geo = ff.geo();
//  WilsonVector* quda_pt = reinterpret_cast<WilsonVector*>(qff.data());
//  int Ls = geo.multiplicity;
//  qassert(Ls > 0);
//  long V = geo.local_volume();
//  long Vh = V / 2;
////
//#pragma omp parallel for
//  for (long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
//    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
//    int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
//    const Vector<WilsonVector> wvs = ff.get_elems_const(xl);
//    for (int s = 0; s < Ls; s++) {
//      int quda_idx = eo * Vh * Ls + s * Vh + qlat_idx_4d / 2;
//      quda_pt[quda_idx] = wvs[s];
//    }
//  }
//}

//template <typename Float> void constructCloverField(Float *res, double norm, double diag, long V)
//{
//
//  Float c = 2.0 * norm / RAND_MAX;
//
//  for (int i = 0; i < V; i++) {
//    for (int j = 0; j < 72; j++) { res[i * 72 + j] = c * rand() - norm; }
//
//    // impose clover symmetry on each chiral block
//    for (int ch = 0; ch < 2; ch++) {
//      res[i * 72 + 3 + 36 * ch] = -res[i * 72 + 0 + 36 * ch];
//      res[i * 72 + 4 + 36 * ch] = -res[i * 72 + 1 + 36 * ch];
//      res[i * 72 + 5 + 36 * ch] = -res[i * 72 + 2 + 36 * ch];
//      res[i * 72 + 30 + 36 * ch] = -res[i * 72 + 6 + 36 * ch];
//      res[i * 72 + 31 + 36 * ch] = -res[i * 72 + 7 + 36 * ch];
//      res[i * 72 + 32 + 36 * ch] = -res[i * 72 + 8 + 36 * ch];
//      res[i * 72 + 33 + 36 * ch] = -res[i * 72 + 9 + 36 * ch];
//      res[i * 72 + 34 + 36 * ch] = -res[i * 72 + 16 + 36 * ch];
//      res[i * 72 + 35 + 36 * ch] = -res[i * 72 + 17 + 36 * ch];
//    }
//
//    for (int j = 0; j < 6; j++) {
//      res[i * 72 + j] += diag;
//      res[i * 72 + j + 36] += diag;
//    }
//  }
//}

//void constructQudaCloverField(void *clover, double norm, double diag, QudaPrecision precision, long V)
//{
//  if (precision == QUDA_DOUBLE_PRECISION)
//    constructCloverField((double *)clover, norm, diag, V);
//  else
//    constructCloverField((float *)clover, norm, diag,  V);
//}

// Helper functions
//------------------------------------------------------
inline bool isPCSolution(QudaSolutionType solution_type)
{
  return (solution_type == QUDA_MATPC_SOLUTION || solution_type == QUDA_MATPC_DAG_SOLUTION
          || solution_type == QUDA_MATPCDAG_MATPC_SOLUTION);
}
//------------------------------------------------------

void constructWilsonTestSpinorParam(quda::ColorSpinorParam *cs_param, const QudaInvertParam *inv_param,
                                    const QudaGaugeParam *gauge_param)
{
  // Lattice vector spacetime/colour/spin/parity properties
  cs_param->nColor = 3;
  cs_param->nSpin = 4;
  if (inv_param->dslash_type == QUDA_DOMAIN_WALL_DSLASH || inv_param->dslash_type == QUDA_DOMAIN_WALL_4D_DSLASH
      || inv_param->dslash_type == QUDA_MOBIUS_DWF_DSLASH || inv_param->dslash_type == QUDA_MOBIUS_DWF_EOFA_DSLASH) {
    cs_param->nDim = 5;
    cs_param->x[4] = inv_param->Ls;
  } else if ((inv_param->dslash_type == QUDA_TWISTED_MASS_DSLASH || inv_param->dslash_type == QUDA_TWISTED_CLOVER_DSLASH)
             && (inv_param->twist_flavor == QUDA_TWIST_NONDEG_DOUBLET)) {
    cs_param->nDim = 5;
    cs_param->x[4] = 2;
  } else {
    cs_param->nDim = 4;
  }
  cs_param->pc_type = inv_param->dslash_type == QUDA_DOMAIN_WALL_DSLASH ? QUDA_5D_PC : QUDA_4D_PC;
  for (int d = 0; d < 4; d++) cs_param->x[d] = gauge_param->X[d];
  bool pc = isPCSolution(inv_param->solution_type);
  if (pc) cs_param->x[0] /= 2;
  cs_param->siteSubset = pc ? QUDA_PARITY_SITE_SUBSET : QUDA_FULL_SITE_SUBSET;

  // Lattice vector data properties
  cs_param->setPrecision(inv_param->cpu_prec);
  cs_param->pad = 0;
  cs_param->siteOrder = QUDA_EVEN_ODD_SITE_ORDER;
  cs_param->fieldOrder = QUDA_SPACE_SPIN_COLOR_FIELD_ORDER;
  cs_param->gammaBasis = inv_param->gamma_basis;
  cs_param->create = QUDA_ZERO_FIELD_CREATE;
  cs_param->location = QUDA_CPU_FIELD_LOCATION;
}

//// data reordering routines
//template <typename Out, typename In>
//void reorderQDPtoMILC(Out *milc_out, In *qdp_in, long V, int siteSize)
//{
//  qthread_for(i, V,{
//    for (int dir = 0; dir < 4; dir++) {
//      for (int j = 0; j < siteSize; j++) {
//        milc_out[(i * 4 + dir) * siteSize + j] = static_cast<Out>(qdp_in[dir][i * siteSize + j]);
//      }
//    }
//  });
//}

bool last_node_in_t()
{
  using namespace quda;
  // only apply T-boundary at edge nodes
  return commCoords(3) == commDim(3) - 1; 
}

int fullLatticeIndex(int dim[4], int index, int oddBit)
{

  int za = index / (dim[0] >> 1);
  int zb = za / dim[1];
  int x2 = za - zb * dim[1];
  int x4 = zb / dim[2];
  int x3 = zb - x4 * dim[2];

  return 2 * index + ((x2 + x3 + x4 + oddBit) & 1);
}

// given a "half index" i into either an even or odd half lattice (corresponding
// to oddBit = {0, 1}), returns the corresponding full lattice index.
// TODO check Z is the dimension
int fullLatticeIndex(int i, int oddBit, int Z[4])
{
  /*
    int boundaryCrossings = i/(Z[0]/2) + i/(Z[1]*Z[0]/2) + i/(Z[2]*Z[1]*Z[0]/2);
    return 2*i + (boundaryCrossings + oddBit) % 2;
  */

  int X1 = Z[0];
  int X2 = Z[1];
  int X3 = Z[2];
  // int X4 = Z[3];
  int X1h = X1 / 2;

  int sid = i;
  int za = sid / X1h;
  // int x1h = sid - za*X1h;
  int zb = za / X2;
  int x2 = za - zb * X2;
  int x4 = zb / X3;
  int x3 = zb - x4 * X3;
  int x1odd = (x2 + x3 + x4 + oddBit) & 1;
  // int x1 = 2*x1h + x1odd;
  int X = 2 * sid + x1odd;

  return X;
}


template <typename Ty>
void applyGaugeFieldScaling_long(Ty *gauge, long Vh, QudaGaugeParam *param, QudaDslashType dslash_type = QUDA_STAGGERED_DSLASH)
{
  int X1h = param->X[0] / 2;
  int X1 = param->X[0];
  int X2 = param->X[1];
  int X3 = param->X[2];
  int X4 = param->X[3];

  const long V = Vh*2;

  // rescale long links by the appropriate coefficient
  if (dslash_type == QUDA_ASQTAD_DSLASH) {
    #pragma omp parallel for
    for(long isp=0;isp < V; isp ++)
    for (int d = 0; d < 4 * 9; d++) {
      {
        gauge[isp*4*9 + d] /= (-24 * param->tadpole_coeff * param->tadpole_coeff);
      }
    }
  }

  // apply the staggered phases
  for (int d = 0; d < 3; d++) {

    // even
    #pragma omp parallel for
    for (int i = 0; i < Vh; i++) {

      int index = fullLatticeIndex(i, 0, param->X);
      int i4 = index / (X3 * X2 * X1);
      int i3 = (index - i4 * (X3 * X2 * X1)) / (X2 * X1);
      int i2 = (index - i4 * (X3 * X2 * X1) - i3 * (X2 * X1)) / X1;
      int i1 = index - i4 * (X3 * X2 * X1) - i3 * (X2 * X1) - i2 * X1;
      int sign = 1;

      if (d == 0) {
        if (i4 % 2 == 1) { sign = -1; }
      }

      if (d == 1) {
        if ((i4 + i1) % 2 == 1) { sign = -1; }
      }
      if (d == 2) {
        if ((i4 + i1 + i2) % 2 == 1) { sign = -1; }
      }

      long indexe = ((0 * Vh + i)*4 + d)*9 ;
      for (int j = 0; j < 9; j++) { gauge[indexe + j] *= sign; }
    }
    // odd
    #pragma omp parallel for
    for (int i = 0; i < Vh; i++) {
      int index = fullLatticeIndex(i, 1, param->X);
      int i4 = index / (X3 * X2 * X1);
      int i3 = (index - i4 * (X3 * X2 * X1)) / (X2 * X1);
      int i2 = (index - i4 * (X3 * X2 * X1) - i3 * (X2 * X1)) / X1;
      int i1 = index - i4 * (X3 * X2 * X1) - i3 * (X2 * X1) - i2 * X1;
      int sign = 1;

      if (d == 0) {
        if (i4 % 2 == 1) { sign = -1; }
      }

      if (d == 1) {
        if ((i4 + i1) % 2 == 1) { sign = -1; }
      }
      if (d == 2) {
        if ((i4 + i1 + i2) % 2 == 1) { sign = -1; }
      }

      long indexo = ((1 * Vh + i)*4 + d)*9 ;
      ////for (int j = 0; j < 18; j++) { gauge[d][(Vh + i) * gauge_site_size + j] *= sign; }
      for (int j = 0; j < 9; j++) { gauge[indexo + j] *= sign; }
    }
  }

  // Apply boundary conditions to temporal links
  if (param->t_boundary == QUDA_ANTI_PERIODIC_T && last_node_in_t()) {
    #pragma omp parallel for
    for (int j = 0; j < Vh; j++) {
      int sign = 1;
      if (dslash_type == QUDA_ASQTAD_DSLASH) {
        if (j >= (X4 - 3) * X1h * X2 * X3) { sign = -1; }
      } else {
        if (j >= (X4 - 1) * X1h * X2 * X3) { sign = -1; }
      }

      long indexe = ((0 * Vh + j)*4 + 3)*9 ;
      long indexo = ((1 * Vh + j)*4 + 3)*9 ;
      for (int i = 0; i < 9; i++) {
        gauge[indexe + i] *= sign;
        gauge[indexo + i] *= sign;
      }
    }
  }
}


template <class Ty, class T1>
void Ffield4d_to_quda_ff(T1*  quda_ff, qlat::FermionField4dT<Ty>& ff, int dir = 1)
{
  TIMER("Ffield4d_to_quda_ff(ff, quda_ff)");
  qassert(ff.initialized);
  const Geometry& geo = ff.geo();
  //const WilsonVector* quda_pt =
  //    reinterpret_cast<const WilsonVector*>(qff.data());
  long V = geo.local_volume();
  long Vh = V / 2;

  Ty* src = (Ty*) qlat::get_data(ff).data();

  //
  #pragma omp parallel for
  for (long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
    int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    for (int dc = 0; dc < 12; dc++)
    {
      long quda_idx = eo * Vh + qlat_idx_4d / 2;
      if(dir == 1){quda_ff[quda_idx*12 + dc] = src[qlat_idx_4d*12 + dc];}
      if(dir == 0){src[qlat_idx_4d*12 + dc] = quda_ff[quda_idx*12 + dc];}
    }
  }

}

template <class Ty, class T1>
void quda_ff_to_Ffield4d(qlat::FermionField4dT<Ty>& ff, T1* quda_ff)
{
  Ffield4d_to_quda_ff(quda_ff, ff, 0);
}

template <class Ty, class T1>
void qlat_cf_to_quda_cf(T1*  quda_cf, colorFT& qlat_cf, int dir = 1)
{
  TIMER("qlat_cf_to_quda_cf(qlat_cf, quda_cf)");
  qassert(qlat_cf.initialized);
  const Geometry& geo = qlat_cf.geo();
  long V = geo.local_volume();
  long Vh = V / 2;
  qassert(geo.multiplicity == 3);
  const long Dim = geo.multiplicity;

  Ty* src = (Ty*) qlat::get_data(qlat_cf).data();
  //
  #pragma omp parallel for
  for (long qlat_idx_4d = 0; qlat_idx_4d < V; qlat_idx_4d++) {
    const Coordinate xl = geo.coordinate_from_index(qlat_idx_4d);
    int eo = (xl[0] + xl[1] + xl[2] + xl[3]) % 2;
    for (int dc = 0; dc < Dim; dc++)
    {
      long quda_idx = eo * Vh + qlat_idx_4d / 2;
      if(dir == 1){quda_cf[quda_idx*Dim + dc] = src[qlat_idx_4d*Dim + dc];}
      if(dir == 0){src[qlat_idx_4d*Dim + dc] = quda_cf[quda_idx*Dim + dc];}
    }
  }

}

template <class Ty, class T1>
void quda_cf_to_qlat_cf(colorFT& qlat_cf, T1* quda_cf)
{
  qlat_cf_to_quda_cf(quda_cf, qlat_cf, 0);
}





}  // namespace qlat

#endif
