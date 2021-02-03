#include "lib.h"

EXPORT(mk_rng, {
  using namespace qlat;
  RngState* prng = NULL;
  PyObject* p_seed = NULL;
  if (!PyArg_ParseTuple(args, "|lO", &prng, &p_seed)) {
    return NULL;
  }
  std::string seed;
  if (p_seed != NULL) {
    py_convert(seed, p_seed);
  }
  RngState* prng_new = NULL;
  if (NULL == prng) {
    prng_new = new RngState();
  } else {
    prng_new = new RngState(*prng, seed);
  }
  return py_convert((void*)prng_new);
});

EXPORT(free_rng, {
  using namespace qlat;
  RngState* prng = NULL;
  if (!PyArg_ParseTuple(args, "l", &prng)) {
    return NULL;
  }
  pqassert(prng != NULL);
  delete prng;
  Py_RETURN_NONE;
});

EXPORT(rand_gen, {
  using namespace qlat;
  RngState* prng = NULL;
  if (!PyArg_ParseTuple(args, "l", &prng)) {
    return NULL;
  }
  pqassert(prng != NULL);
  RngState& rng = *prng;
  const uint64_t x = rand_gen(rng);
  return py_convert(x);
});

EXPORT(u_rand_gen, {
  using namespace qlat;
  RngState* prng = NULL;
  double upper = 1.0;
  double lower = 0.0;
  if (!PyArg_ParseTuple(args, "l|dd", &prng, &upper, &lower)) {
    return NULL;
  }
  pqassert(prng != NULL);
  RngState& rng = *prng;
  const double x = u_rand_gen(rng, upper, lower);
  return py_convert(x);
});

EXPORT(g_rand_gen, {
  using namespace qlat;
  RngState* prng = NULL;
  double center = 0.0;
  double sigma = 1.0;
  if (!PyArg_ParseTuple(args, "l|dd", &prng, &center, &sigma)) {
    return NULL;
  }
  pqassert(prng != NULL);
  RngState& rng = *prng;
  const double x = g_rand_gen(rng, center, sigma);
  return py_convert(x);
});

