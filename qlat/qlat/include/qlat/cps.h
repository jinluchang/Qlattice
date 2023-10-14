#pragma once

// Use with CPS ( https://github.com/RBC-UKQCD/CPS_public ).
// Program needs arguments like: "-qmp-geom 2 2 2 4" to specify the MPI geometry.

#include <qlat/qlat.h>

#ifdef NO_CPS

namespace qlat
{  //

typedef InverterDomainWall InverterDomainWallCPS;

inline void cps_begin(int* argc, char** argv[], const Coordinate& total_site)
{
  begin(argc, argv);
}

inline void cps_end(const bool is_preserving_cache = false)
{
  end(is_preserving_cache);
}

}  // namespace qlat

#else

#define QLAT_CPS

#include <alg/alg_fix_gauge.h>
#include <alg/alg_rnd_gauge.h>
#include <util/gjp.h>

extern MPI_Comm QMP_COMM_WORLD;

namespace qlat
{  //

typedef InverterDomainWall InverterDomainWallCPS;

inline void set_do_arg(cps::DoArg& do_arg, const Coordinate& total_site)
{
  using namespace cps;
  do_arg.x_sites = total_site[0];
  do_arg.y_sites = total_site[1];
  do_arg.z_sites = total_site[2];
  do_arg.t_sites = total_site[3];
  do_arg.s_sites = 2;
  do_arg.dwf_height = 1.0;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_INPUT;
  do_arg.start_seed_value = 123121;
  do_arg.x_nodes = 0;
  do_arg.y_nodes = 0;
  do_arg.z_nodes = 0;
  do_arg.t_nodes = 0;
  do_arg.s_nodes = 0;
  do_arg.x_node_sites = 0;
  do_arg.y_node_sites = 0;
  do_arg.z_node_sites = 0;
  do_arg.t_node_sites = 0;
  do_arg.s_node_sites = 0;
  do_arg.gfix_chkb = 1;
}

inline void cps_begin(int* argc, char** argv[], const Coordinate& total_site, const bool is_initialize_rng = false)
{
  cps::Start(argc, argv);
  cps::DoArg do_arg;
  set_do_arg(do_arg, total_site);
  cps::GJP.Initialize(do_arg);
  if (is_initialize_rng) {
    cps::LRG.Initialize();
  }
  Coordinate size_node(cps::SizeX(), cps::SizeY(), cps::SizeZ(), cps::SizeT());
  begin(cps::UniqueID(), size_node);
  sync_node();
}

inline void cps_end(const bool is_preserving_cache = false)
{
  end(is_preserving_cache);
  cps::End();
  // Currently CPS does not call destroy_qmp() in End(), so call it explicitly now.
  cps::QMPSCU::destroy_qmp();
}

template <class M, class N>
inline void value_convert(M& x, const N& y)
{
  qassert(sizeof(M) == sizeof(N));
  std::memcpy(&x, &y, sizeof(M));
}

}  // namespace qlat

#endif
