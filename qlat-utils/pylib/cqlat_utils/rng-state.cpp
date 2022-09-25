#include "lib.h"

EXPORT(mk_rng, {
  using namespace qlat;
  PyObject* p_rng = NULL;
  PyObject* p_seed = NULL;
  if (!PyArg_ParseTuple(args, "|OO", &p_rng, &p_seed)) {
    return NULL;
  }
  RngState* prng_new = NULL;
  if (NULL == p_rng) {
    prng_new = new RngState();
  } else {
    RngState& rng = py_convert_type<RngState>(p_rng);
    if (p_seed == NULL) {
      prng_new = new RngState(rng);
    } else {
      const std::string seed = py_convert_data<std::string>(p_seed);
      prng_new = new RngState(rng, seed);
    }
  }
  return py_convert((void*)prng_new);
})

EXPORT(free_rng, {
  using namespace qlat;
  return free_obj<RngState>(args);
})

EXPORT(set_rng, {
  using namespace qlat;
  return set_obj<RngState>(args);
})

EXPORT(rand_gen, {
  using namespace qlat;
  PyObject* p_rng = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_rng)) {
    return NULL;
  }
  RngState& rng = py_convert_type<RngState>(p_rng);
  const uint64_t x = rand_gen(rng);
  return py_convert(x);
})

EXPORT(u_rand_gen, {
  using namespace qlat;
  PyObject* p_rng = NULL;
  double upper = 1.0;
  double lower = 0.0;
  if (!PyArg_ParseTuple(args, "O|dd", &p_rng, &upper, &lower)) {
    return NULL;
  }
  RngState& rng = py_convert_type<RngState>(p_rng);
  const double x = u_rand_gen(rng, upper, lower);
  return py_convert(x);
})

EXPORT(g_rand_gen, {
  using namespace qlat;
  PyObject* p_rng = NULL;
  double center = 0.0;
  double sigma = 1.0;
  if (!PyArg_ParseTuple(args, "O|dd", &p_rng, &center, &sigma)) {
    return NULL;
  }
  RngState& rng = py_convert_type<RngState>(p_rng);
  const double x = g_rand_gen(rng, center, sigma);
  return py_convert(x);
})
