#ifndef UTILS_MCG_SIGNATURE_H
#define UTILS_MCG_SIGNATURE_H
#pragma once

// Diagnostic control: qlat_mcg_signature=0 (default) disables automatic
// fingerprints; any nonzero value enables them. Read once per process.
// Enabled calls allocate host buffers, copy GPU data, evaluate CPU probes,
// perform MPI sums and print records. These observations do not modify the
// solve, but are not free even when excluded from the reported solve timer.
// Automatic callers guard all work with mcg_signature_enabled(); direct
// helper calls below explicitly request an observation. See the flag table
// at the top of utils_staggered_multishift_mcg.h for the complete reference.

#include "general_funs.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace qlat {

// A deliberately small, reusable observer for custom staggered MCG data.
//
// Public staggered vectors contain three colors on the even checkerboard:
// the logical site `s` is the even full-lattice site whose lexicographic
// index is either 2*s or 2*s+1, and `staggered_mcg_layout::public_index`
// maps (s,color) to the native public storage offset.  The custom MCG
// workspace uses `storage_index(site,color,sites)` instead; its conversion
// kernel reads public_index and writes storage_index.  The helper therefore
// never treats either buffer as lexicographic storage.  Coordinates are
// global coordinates from Geometry::g_coordinate_from_index, normalized as
// 2*x/(L-1)-1 (and zero for a length-one direction).  The five reductions are
// ordinary norm, coordinate-weighted norm, a signed phase probe, its
// coordinate-weighted form, and a signed staggered/color probe.  Each local
// five-tuple is combined across ranks with sum_all_size. The phase probes
// also use local site indices: compare them with the same MPI decomposition.
//
// The phase and staggered/color probes use deterministic coordinates, color,
// and (for gauge data) direction/matrix indices.  They detect sign and phase
// changes which an ordinary norm cannot.  Matching fingerprints are strong
// evidence that two vectors agree, but are not a proof of vector equality.
struct MCGSignature {
  double l2 = 0;
  double coord_l2 = 0;
  double probe0 = 0;
  double probe1 = 0;
  double probe2 = 0;
};

// Process control and record format
inline bool mcg_signature_enabled() {
  static const bool enabled = qlat::get_env_long_default(std::string("qlat_mcg_signature"), 0) != 0;
  return enabled;
}

inline void mcg_signature_config_once() {
  if (!mcg_signature_enabled())
    return;
  static const bool printed = [] {
    qmessage("MCG_SIG_CONFIG version=1 basis=coord-normalized probes=phase-signed global=1\n");
    return true;
  }();
  (void)printed;
}

// Coordinate weights and deterministic projections
inline double mcg_signature_normalized_coordinate(const Coordinate& x, const Coordinate& extent, Int mu) {
  return extent[mu] > 1 ? 2.0 * double(x[mu]) / double(extent[mu] - 1) - 1.0 : 0.0;
}

inline double mcg_signature_weight(const Coordinate& x, const Coordinate& extent) {
  double sum = 0;
  for (Int mu = 0; mu < 4; ++mu)
    sum += mcg_signature_normalized_coordinate(x, extent, mu);
  return 1.0 + 0.25 * sum / 4.0;
}

inline double mcg_signature_projection(const ComplexD& value, double phase) {
  return value.real() * std::cos(phase) - value.imag() * std::sin(phase);
}

inline void mcg_signature_add(MCGSignature& out, const ComplexD& value, double weight, double phase,
                              double signed_phase) {
  const double norm = value.real() * value.real() + value.imag() * value.imag();
  out.l2 += norm;
  out.coord_l2 += weight * norm;
  out.probe0 += mcg_signature_projection(value, phase);
  out.probe1 += weight * mcg_signature_projection(value, phase);
  out.probe2 += signed_phase * mcg_signature_projection(value, phase + 0.37);
}

inline MCGSignature mcg_signature_global(MCGSignature local) {
  double values[5] = {local.l2, local.coord_l2, local.probe0, local.probe1, local.probe2};
  sum_all_size(values, 5, 0);
  return {values[0], values[1], values[2], values[3], values[4]};
}

// Return the full-lattice coordinate belonging to an even checkerboard site.
// This is the same 2*s / (2*s+1) mapping used by the qlat QUDA conversion.
inline Coordinate mcg_signature_even_coordinate(const Geometry& geo, Long site) {
  const Long first = 2 * site;
  Coordinate x = geo.coordinate_from_index(first);
  if ((x[0] + x[1] + x[2] + x[3]) % 2 != 0)
    x = geo.coordinate_from_index(first + 1);
  return geo.coordinate_g_from_l(x);
}

// Explicit vector and gauge observations
template <class U, class Layout>
inline MCGSignature mcg_signature_public_vector(const Geometry& geo, const Layout& layout,
                                                const U* device_data, Long sites) {
  const Long elements = sites * 3;
  std::vector<U> host(elements);
  if (elements)
    cpy_GPU(host.data(), device_data, elements, 0, 1);
  const Coordinate extent = geo.total_site();
  MCGSignature local;
  for (Long site = 0; site < sites; ++site) {
    const Coordinate x = mcg_signature_even_coordinate(geo, site);
    const double weight = mcg_signature_weight(x, extent);
    const double phase = 0.173 * (1.0 + double(site)) + 0.271 * double(x[0] + 2 * x[1] + 3 * x[2] + 5 * x[3]);
    for (Int color = 0; color < 3; ++color) {
      const U value = host[layout.template public_index<U>(site, color, sites)];
      const ComplexD v = ComplexD(value);
      const double signed_phase = ((site + color + (x[0] + x[1] + x[2] + x[3])) & 1) ? -1.0 : 1.0;
      mcg_signature_add(local, v, weight, phase + 0.131 * double(color), signed_phase);
    }
  }
  return mcg_signature_global(local);
}

inline MCGSignature mcg_signature_gauge(const Geometry& geo, const ComplexD* device_data, Long elements) {
  std::vector<ComplexD> host(elements);
  if (elements)
    cpy_GPU(host.data(), device_data, elements, 0, 1);
  const Long volume = geo.local_volume(), half = volume / 2;
  const Coordinate extent = geo.total_site();
  MCGSignature local;
  for (Long full = 0; full < volume; ++full) {
    const Coordinate x = geo.g_coordinate_from_index(full);
    const double weight = mcg_signature_weight(x, extent);
    const Int parity = (x[0] + x[1] + x[2] + x[3]) & 1;
    const Long site = parity * half + full / 2;
    for (Int dir = 0; dir < 4; ++dir)
      for (Int row = 0; row < 3; ++row)
        for (Int col = 0; col < 3; ++col) {
          const Long index = ((site * 4 + dir) * 3 + row) * 3 + col;
          if (index >= elements)
            continue;
          const double phase =
              0.097 * (1.0 + double(full)) + 0.193 * double(dir) + 0.071 * double(row) + 0.113 * double(col);
          const double signed_phase = ((full + dir + row + 2 * col) & 1) ? -1.0 : 1.0;
          mcg_signature_add(local, host[index], weight, phase, signed_phase);
        }
  }
  return mcg_signature_global(local);
}

inline void mcg_signature_print_vector(Long call, const char* stage, Int rhs, Long shift,
                                       const MCGSignature& sig) {
  if (shift < 0)
    qmessage("  MCG_VEC_SIG call=%ld stage=%s rhs=%d l2=%.17e coord_l2=%.17e probe0=%.17e probe1=%.17e "
             "probe2=%.17e\n",
             long(call), stage, rhs, sig.l2, sig.coord_l2, sig.probe0, sig.probe1, sig.probe2);
  else
    qmessage("  MCG_VEC_SIG call=%ld stage=%s rhs=%d shift=%ld l2=%.17e coord_l2=%.17e probe0=%.17e "
             "probe1=%.17e probe2=%.17e\n",
             long(call), stage, rhs, long(shift), sig.l2, sig.coord_l2, sig.probe0, sig.probe1, sig.probe2);
}

} // namespace qlat
#endif
